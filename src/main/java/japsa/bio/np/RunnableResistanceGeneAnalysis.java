package japsa.bio.np;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonObject;
import htsjdk.samtools.*;
import japsa.bio.alignment.ProbFSM;
import japsa.bio.detection.DetectedGeneRecord;
import japsa.seq.Alphabet;
import japsa.seq.Sequence;
import japsa.seq.SequenceOutputStream;
import japsa.util.HTSUtilities;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.*;

public class RunnableResistanceGeneAnalysis extends RunnableAnalysis {
  private static final Logger LOG = LoggerFactory.getLogger(RunnableResistanceGeneAnalysis.class);
  private final Map<String, List<Sequence>> alignmentMapSnap = new HashMap<>();
  private static final Gson gson = new GsonBuilder().serializeNulls().create();
  private String recordPrefix = "tmp";
  private final RealtimeResistanceGene resistGene;
  private SequenceOutputStream sequenceOutputStream;
  InputStream samInputStream;

  //Set of genes confirmed to have found
  private final HashSet<String> predictedGenes = new HashSet<>();

  private int runIndex = 0;

  public RunnableResistanceGeneAnalysis(RealtimeResistanceGene resistGene, String outputFile, String recordPrefix) throws IOException {
    this(resistGene, new FileOutputStream(outputFile), recordPrefix);
  }

  RunnableResistanceGeneAnalysis(RealtimeResistanceGene resistGene, OutputStream outStream, String recordPrefix) throws IOException {
    this.resistGene = resistGene;
    this.recordPrefix = recordPrefix;
    this.sequenceOutputStream = new SequenceOutputStream(outStream);
  }

  @Override
  public void analysis() {
    SamReaderFactory.setDefaultValidationStringency(ValidationStringency.SILENT);
    SamReader samReader = SamReaderFactory.makeDefault().open(SamInputResource.of(samInputStream));
    SAMRecordIterator samIter = samReader.iterator();

    String readName = "";
    //A dummy sequence
    Sequence readSequence = new Sequence(Alphabet.DNA(), 1, "");

    while (samIter.hasNext()) {
      LOG.info("next SAM record");
      SAMRecord record = samIter.next();
      String geneID = record.getReferenceName();
      int refLength = RealtimeResistanceGene.getGeneSequence(geneID).length();
      int refStart = record.getAlignmentStart();
      int refEnd = record.getAlignmentEnd();

      if (RealtimeResistanceGene.twoDOnly && !record.getReadName().contains("twodim"))
        continue;

      if (!record.getReadName().equals(readName)) {
        readName = record.getReadName();

        RealtimeResistanceGene.observedReadCount++;
        RealtimeResistanceGene.observedBaseCount += record.getReadLength();

        //Get the read
        if (!record.getReadUnmappedFlag()) {
          readSequence = new Sequence(Alphabet.DNA(), record.getReadString(), readName);
          if (record.getReadNegativeStrandFlag()) {
            readSequence = Alphabet.DNA.complement(readSequence);
            readSequence.setName(readName);
          }
        }
      }

      if (record.getReadUnmappedFlag())
        continue;
      if (RealtimeResistanceGene.getGeneSequence(geneID) == null)
        continue;
      if (refStart > 99 || refEnd < refLength - 99)
        continue;

      synchronized (resistGene) {
        if (resistGene.getAlignmentMap().get(geneID) == null)
          resistGene.getAlignmentMap().put(geneID, new ArrayList<>());

        //put the sequence into alignment list
        Sequence readSeq = HTSUtilities.readSequence(record, readSequence, 99, refLength - 99);
        resistGene.getAlignmentMap().get(geneID).add(readSeq);
      }

      try {
        if (!RealtimeResistanceGene.JSON) {
          sequenceOutputStream.print("##" + timeNow + "\t" + (this.lastTime - this.startTime) + "\t" + RealtimeResistanceGene.observedReadCount + "\n");
        } else {
          LOG.info("status update");

          JsonObject jo = new JsonObject();
          jo.addProperty("timestamp", timeNow);
          jo.addProperty("timeLast", this.lastTime);
          jo.addProperty("timeStart", this.startTime);
          jo.addProperty("timeWaited", (this.lastTime - this.startTime));
          jo.addProperty("currentReadCount", RealtimeResistanceGene.observedReadCount);
          sequenceOutputStream.print(gson.toJson(jo));
          sequenceOutputStream.println();
        }
        sequenceOutputStream.flush();

        //1. Make a snapshot of the current alignment
        synchronized (resistGene) {
          for (String gene : resistGene.getAlignmentMap().keySet()) {
            List<Sequence> readMap = resistGene.getAlignmentMap().get(gene);
            List<Sequence> sequences = new ArrayList<>();
            sequences.addAll(readMap);
            alignmentMapSnap.put(gene, sequences);
          }
        }

        runIndex++;
        //Now can make the call
        antiBioticsProfile();
      } catch (IOException | InterruptedException e) {
        e.printStackTrace();
      }
    }
  }

  private void antiBioticsProfile() throws IOException, InterruptedException {
    for (String geneID : RealtimeResistanceGene.getGeneSequences()) {
      if (predictedGenes.contains(geneID))
        continue;

      List<Sequence> alignmentList = alignmentMapSnap.get(geneID);
      Sequence consensus = ErrorCorrection.consensusSequence(alignmentList, recordPrefix + "_" + geneID + "_thread-" + runIndex, RealtimeResistanceGene.msa);

      if (consensus == null)
        continue;

      Sequence gene = RealtimeResistanceGene.getGeneSequence(geneID);
      gene = gene.subSequence(99, gene.length() - 100);

      ProbFSM.ProbOneSM tsmF = new ProbFSM.ProbOneSM(gene);
      double cost = 100000000;
      for (int c = 0; c < 10; c++) {
        tsmF.resetCount();
        ProbFSM.Emission retState = tsmF.alignGenerative(consensus);
        if (cost <= retState.myCost)
          break;//for c

        cost = retState.myCost;
        int emitCount = tsmF.updateCount(retState);
        LOG.info("Iter " + c + " : " + emitCount + " states and " + cost + " bits " + consensus.length() + "bp " + consensus.getName() + " by " + gene.getName());
        tsmF.reEstimate();
      }
      double score = (consensus.length() * 2 - cost) / gene.length();
      LOG.info("SGF: score=" + score + " geneID=" + geneID + " group=" + RealtimeResistanceGene.getGeneGroup(geneID));

      if (score >= RealtimeResistanceGene.scoreThreshold) {
        LOG.info("ADDF " + geneID);
        predictedGenes.add(geneID);
        if (!RealtimeResistanceGene.JSON) {
          sequenceOutputStream.print(timeNow +
              "\t" + (this.lastTime - this.startTime) / 1000 +
              "\t" + RealtimeResistanceGene.observedReadCount +
              "\t" + RealtimeResistanceGene.observedBaseCount +
              "\t" + geneID +
              "\t" + RealtimeResistanceGene.getGeneName(geneID) +
              "\t" + RealtimeResistanceGene.getGeneGroup(geneID) +
              "\n");
        } else {
          DetectedGeneRecord gr = RealtimeResistanceGene.createDetectedGeneRecord(geneID);
          JsonObject jo = new JsonObject();
          jo.addProperty("timeStamp", timeNow);
          jo.addProperty("timeLast", this.lastTime);
          jo.addProperty("timeStart", this.startTime);
          jo.addProperty("timeWaited", (this.lastTime - this.startTime));
          jo.addProperty("currentBaseCount", RealtimeResistanceGene.observedBaseCount);
          jo.addProperty("currentReadCount", RealtimeResistanceGene.observedReadCount);
          jo.addProperty("geneID", geneID);
          jo.addProperty("geneName", gr.geneName);
          jo.addProperty("geneGroup", gr.geneGroup);

          LOG.info(gson.toJson(jo));
          sequenceOutputStream.print(gson.toJson(jo));
          sequenceOutputStream.println();
        }
        sequenceOutputStream.flush();

      }
    }
  }

  @Override
  protected void close() {
    try {
      sequenceOutputStream.close();
    } catch (IOException e) {
      // TODO Auto-generated catch block
      e.printStackTrace();
    }
  }
  @Override
  protected int getCurrentRead() {
    return RealtimeResistanceGene.observedReadCount;

  }
}
