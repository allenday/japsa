package japsa.bio.np;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonObject;
import japsa.bio.alignment.ProbFSM;
import japsa.seq.*;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

//TODO: way to improve performance:
//1.
//3. options: gene or antibiotics class
//4.
//5. Future improve: incrementally multiple alignment

public class ResistanceGeneAnalysis extends RunnableAnalysis {
  private static final Logger LOG = LoggerFactory.getLogger(ResistanceGeneAnalysis.class);

  private final Map<String, List<Sequence>> alignmentMapSnap = new HashMap<String, List<Sequence>>();
  private final Map<String, String> gene2GeneName = new HashMap<String, String>();
  private final Map<String, String> gene2Group = new HashMap<String, String>();
  private final Map<String, Sequence> geneMap = new HashMap<String, Sequence>();
  private final List<String> geneList = new ArrayList<String>();

  public Boolean JSON = false;
  private static final Gson gson = new GsonBuilder().serializeNulls().create();

  //Set of genes confirmed to have found
  private final HashSet<String> predictedGenes = new HashSet<String>();

  private String recordPrefix = "tmp";
  private RealtimeResistanceGene resistGene;
  private SequenceOutputStream sequenceOutputStream;

  public int runIndex = 0;

  public ResistanceGeneAnalysis(RealtimeResistanceGene resistGene, String outputFile, String resDB, String recordPrefix) throws IOException {
    this(resistGene, new FileOutputStream(outputFile), resDB, recordPrefix);
  }

  public ResistanceGeneAnalysis(RealtimeResistanceGene resistGene, OutputStream outStream, String resDB, String recordPrefix) throws IOException {
    this(resistGene, outStream, new FileInputStream(new File(resDB + "/geneList")), new FileInputStream(new File(resDB + "/DB.fasta")), recordPrefix);
  }

  public ResistanceGeneAnalysis(RealtimeResistanceGene resistGene, OutputStream outStream, InputStream resDBInputStream, InputStream fastaInputStream, String recordPrefix) throws IOException {
    this.resistGene = resistGene;

    readGeneClassInformation(fastaInputStream); // step 1
    readGeneInformation(resDBInputStream);      // step 2

    LOG.info("geneList = " + geneList.size());
    LOG.info("geneMap = " + geneMap.size());
    LOG.info("gene2Group = " + gene2Group.size());
    LOG.info("gene2GeneName = " + gene2GeneName.size());

    this.recordPrefix = recordPrefix;
    this.sequenceOutputStream = new SequenceOutputStream(outStream);
  }

  private void readGeneClassInformation(InputStream fastaInputStream) throws IOException {
    List<Sequence> drGeneList = new ArrayList<Sequence>();

    SequenceReader sequenceReader = FastaReader.getReader(fastaInputStream);

    while (true) {
      Sequence fastaSequence = sequenceReader.nextSequence(Alphabet.DNA());
      if (fastaSequence == null)
        break;
      //LOG.info("fastaSequence = "+fastaSequence);
      drGeneList.add(fastaSequence);
    }
    for (Sequence seq : drGeneList) {
      geneMap.put(seq.getName(), seq);
      geneList.add(seq.getName());

      String desc = seq.getDesc();
      String[] toks = desc.split(";");
      for (String tok : toks) {
        if (tok.startsWith("dg=")) {
          addGeneInfo(gene2Group, seq.getName(), tok.substring(3));
        }
        if (tok.startsWith("geneID=")) {
          String proteinID = tok.substring(7);
          addGeneInfo(gene2GeneName, seq.getName(), proteinID);
        }
      }
    }
  }

  private void readGeneInformation(InputStream geneInfoInputStream) throws IOException {
    BufferedReader bf = SequenceReader.openInputStream(geneInfoInputStream);
    String line = "";
    while ((line = bf.readLine()) != null) {
      String[] toks = line.trim().split(" ");
      if (toks.length < 3)
        continue;

      addGeneInfo(gene2Group, toks[0], toks[2]);
      addGeneInfo(gene2GeneName, toks[0], toks[1]);

    }
    bf.close();
  }

  private void addGeneInfo(Map<String, String> map, String key, String info) {
    String s = map.get(key);
    if (s == null)
      map.put(key, info);
    else {
      if (!s.contains(info)) {
        s = s + ", " + info;
        map.put(key, s);
      }
    }
  }

  @SuppressWarnings("unchecked")
  private void antiBioticAnalysis() {
    try {
      if (!JSON) {
        sequenceOutputStream.print("##" + timeNow + "\t" + (this.lastTime - this.startTime) + "\t" + this.lastReadNumber + "\n");
      } else {
        JsonObject jo = new JsonObject();
        jo.addProperty("timestamp", timeNow);
        jo.addProperty("timeLast", this.lastTime);
        jo.addProperty("timeStart", this.startTime);
        jo.addProperty("timeWaited", (this.lastTime - this.startTime));
        jo.addProperty("lastReadNumber", this.lastReadNumber);
        sequenceOutputStream.print(gson.toJson(jo));
        sequenceOutputStream.println();
      }
      sequenceOutputStream.flush();

      //1. Make a snapshot of the current alignment
      synchronized (resistGene) {
        //lastTime = System.currentTimeMillis();
        //lastReadNumber = resistGene.currentReadCount;
        for (String gene : resistGene.getAlignmentMap().keySet()) {
          List<Sequence> readMap = resistGene.getAlignmentMap().get(gene);
          List<Sequence> sequences = new ArrayList<Sequence>();
          sequences.addAll(readMap);
          alignmentMapSnap.put(gene, sequences);
        }
      }//synchronized(resistGene)

      runIndex++;
      //Now can make the call
      antiBioticsProfile();
    } catch (IOException e) {
      e.printStackTrace();
    } catch (InterruptedException e) {
      e.printStackTrace();
    }
  }

  private void antiBioticsProfile() throws IOException, InterruptedException {
    int jobNo = 0;
    //Get list of genes from my
    ExecutorService executor = Executors.newFixedThreadPool((resistGene.numThead > 2) ? resistGene.numThead - 2 : 1);

    for (String geneID : geneList) {
      if (predictedGenes.contains(geneID))
        continue;

      List<Sequence> alignmentList = alignmentMapSnap.get(geneID);
      Sequence consensus = ErrorCorrection.consensusSequence(alignmentList, recordPrefix + "_" + geneID + "_" + runIndex, resistGene.msa);

      if (consensus == null) {
        continue;//gene
      }

      Sequence gene = geneMap.get(geneID);
      gene = gene.subSequence(99, gene.length() - 100);

      RealtimeResistanceGene.FSMThread thread = new RealtimeResistanceGene.FSMThread();
      thread.resGeneFinder = this;
      thread.consensus = consensus;
      thread.gene = gene;
      thread.geneID = geneID;

      executor.execute(thread);
      jobNo++;
    }
    executor.shutdown();
    executor.awaitTermination(3, TimeUnit.DAYS);
    LOG.info("===Found " + predictedGenes.size() + " vs " + geneMap.size() + "  " + alignmentMapSnap.size() + " with " + jobNo);
  }

  private RealtimeResistanceGene.GeneRecord createGeneRecord(String geneID) {
    RealtimeResistanceGene.GeneRecord gr = new RealtimeResistanceGene.GeneRecord();
    gr.geneID = geneID;
    gr.geneName = gene2GeneName.get(geneID);
    gr.geneGroup = gene2Group.get(geneID);
    return gr;
  }

  private void addPredictedGene(String geneID) throws IOException {
    predictedGenes.add(geneID);
    if (!JSON) {
      sequenceOutputStream.print(timeNow + "\t" + (this.lastTime - this.startTime) / 1000 + "\t" + lastReadNumber + "\t" + resistGene.getBasesObserved() + "\t" + geneID + "\t" + gene2GeneName.get(geneID) + "\t" + gene2Group.get(geneID) + "\n");
    } else {

      RealtimeResistanceGene.GeneRecord gr = createGeneRecord(geneID);

      JsonObject jo = new JsonObject();
      jo.addProperty("timeStamp", timeNow);
      jo.addProperty("timeLast", this.lastTime);
      jo.addProperty("timeStart", this.startTime);
      jo.addProperty("timeWaited", (this.lastTime - this.startTime));
      jo.addProperty("lastReadNumber", this.lastReadNumber);
      jo.addProperty("currentBaseCount", resistGene.getBasesObserved());
      jo.addProperty("currentReadCount", resistGene.getReadsObserved());
      jo.addProperty("geneID", geneID);
      jo.addProperty("geneName", gr.geneName);
      jo.addProperty("geneGroup", gr.geneGroup);

      sequenceOutputStream.print(gson.toJson(jo));
      sequenceOutputStream.println();
    }
    sequenceOutputStream.flush();
  }

  private static double fsmAlignment(Sequence consensus, Sequence gene) {
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
    return (consensus.length() * 2 - cost) / gene.length();
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
  protected void analysis() {
    antiBioticAnalysis();
  }

  @Override
  protected int getCurrentRead() {
    return resistGene.getReadsObserved();

  }
}
