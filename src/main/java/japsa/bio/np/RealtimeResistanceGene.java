/*****************************************************************************
 * Copyright (c) Minh Duc Cao, Monash Uni & UQ, All rights reserved.         *
 *                                                                           *
 * Redistribution and use in source and binary forms, with or without        *
 * modification, are permitted provided that the following conditions        *
 * are met:                                                                  * 
 *                                                                           *
 * 1. Redistributions of source code must retain the above copyright notice, *
 *    this list of conditions and the following disclaimer.                  *
 * 2. Redistributions in binary form must reproduce the above copyright      *
 *    notice, this list of conditions and the following disclaimer in the    *
 *    documentation and/or other materials provided with the distribution.   *
 * 3. Neither the names of the institutions nor the names of the contributors*
 *    may be used to endorse or promote products derived from this software  *
 *    without specific prior written permission.                             *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS   *
 * IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, *
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR    *
 * PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR         *
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,     *
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,       *
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR        *
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    *
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      *
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              *
 ****************************************************************************/

/**************************     REVISION HISTORY    **************************
 * 07/09/2014 - Minh Duc Cao: Created                                        
 *
 ****************************************************************************/

package japsa.bio.np;

import japsa.bio.detection.DetectedGeneRecord;
import japsa.seq.*;
import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import htsjdk.samtools.*;
import japsa.util.HTSUtilities;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.*;

/**
 * @author minhduc
 */

public class RealtimeResistanceGene {
  private static final Logger LOG = LoggerFactory.getLogger(RealtimeResistanceGene.class);
  private RunnableResistanceGeneAnalysis resistFinder;
  private HashMap<String, List<Sequence>> alignmentMap = new HashMap<>();
  private static Map<String, String> geneNames = new HashMap<String, String>();

  private static Map<String, String> geneGroups = new HashMap<String, String>();
  private static Map<String, Sequence> geneSequences = new HashMap<String, Sequence>();
  private Integer readsObserved = 0;
  private Long basesObserved = 0L;
  private static final Gson gson = new GsonBuilder().serializeNulls().create();
  static Boolean JSON = false;
  static Long observedBaseCount = 0L;
  static Integer observedReadCount = 0;

  public static String msa = "kalign";

  public static Double scoreThreshold = 2D;
  public static Boolean twoDOnly = false;
  public static Integer numThead = 4;

  Thread analysisThread = null;
  Long startTime = -1L;

  public RealtimeResistanceGene(Integer read, Integer time, String outputFile, String resDB, String recordPrefix) throws IOException {
    File geneListFile = new File(resDB + "/geneList");
    File fastaFile = new File(resDB + "/DB.fasta");
    OutputStream outStream = SequenceOutputStream.makeOutputStream(outputFile);
    InputStream resDBInputStream = new FileInputStream(geneListFile);
    InputStream fastaInputStream = new FileInputStream(fastaFile);

    resistFinder = new RunnableResistanceGeneAnalysis(this, outStream, recordPrefix);
    resistFinder.setReadPeriod(read);
    resistFinder.setTimePeriod(time * 1000);

    loadGeneInfo(fastaInputStream, resDBInputStream);
  }

  RealtimeResistanceGene(Integer read, Integer time, OutputStream outStream, InputStream resDBInputStream, InputStream fastaInputStream, String recordPrefix) throws IOException {
    resistFinder = new RunnableResistanceGeneAnalysis(this, outStream, recordPrefix);
    resistFinder.setReadPeriod(read);
    resistFinder.setTimePeriod(time * 1000);

    loadGeneInfo(fastaInputStream, resDBInputStream);
  }

  private void loadGeneInfo(InputStream fastaInputStream, InputStream resDBInputStream) throws IOException {
    readGeneClassInformation(fastaInputStream); // step 1
    LOG.info("geneGroups = " + getGeneGroupCount());
    LOG.info("geneNames = " + getGeneNameCount());
    LOG.info("geneSequences = " + getGeneSequenceCount());

    readGeneInformation(resDBInputStream);      // step 2
    LOG.info("geneGroups = " + getGeneGroupCount());
    LOG.info("geneNames = " + getGeneNameCount());
    LOG.info("geneSequences = " + getGeneSequenceCount());
  }

  Map<String,List<Sequence>> getAlignmentMap() {
    return alignmentMap;
  }
  static String getGeneGroup(String geneID) {
    return geneGroups.get(geneID);
  }
  static String getGeneName(String geneID) {
    return geneNames.get(geneID);
  }
  public static Set<String> getGeneSequences() {
    return geneSequences.keySet();
  }
  static Sequence getGeneSequence(String geneID) {
    return geneSequences.get(geneID);
  }
  private Integer getGeneSequenceCount() {
    return geneSequences.keySet().size();
  }
  private Integer getGeneGroupCount() {
    return geneGroups.size();
  }
  private Integer getGeneNameCount() {
    return geneNames.size();
  }
  private void putGeneName(String geneID, String geneName) {
    geneNames.put(geneID, geneName);
  }
  private void putGeneGroup(String geneID, String geneGroup) {
    geneGroups.put(geneID, geneGroup);
  }
  private void putGeneSequence(String geneID, Sequence sequence) {
    geneSequences.put(geneID, sequence);
  }

  public static DetectedGeneRecord createDetectedGeneRecord(String geneID) {
    DetectedGeneRecord gr = new DetectedGeneRecord();
    gr.geneID = geneID;
    gr.geneName = getGeneName(geneID);
    gr.geneGroup = getGeneGroup(geneID);
    return gr;
  }

  public void typing(String samFile) throws IOException, InterruptedException {
    InputStream samInputStream;

    if ("-".equals(samFile))
      samInputStream = System.in;
    else
      samInputStream = new FileInputStream(samFile);

    typing(samInputStream);
  }

  public void typing(InputStream samInputStream) throws IOException, InterruptedException {
    LOG.info("Resistance identification ready at " + new Date());
    resistFinder.samInputStream = samInputStream;
    startTime = System.currentTimeMillis();
    analysisThread = new Thread(resistFinder, "analysisThread");
    analysisThread.start();
  }

  public Set<DetectedGeneRecord> getGenes() throws IOException, InterruptedException {
    Set<DetectedGeneRecord> geneRecordSet = new HashSet<>();

    String readName = "";
    //A dummy sequence
    //Sequence readSequence = new Sequence(Alphabet.DNA(), 1, "");

    Long timeLast = System.currentTimeMillis();
    synchronized(this) {
      for (String geneID : getAlignmentMap().keySet()) {
        DetectedGeneRecord gr = new DetectedGeneRecord();
        gr.geneID = geneID;
        gr.geneGroup = getGeneGroup(geneID);
        gr.geneName = getGeneName(geneID);
        gr.timeStart = startTime;
        gr.timeNow = timeLast;
        gr.observedBaseCount = observedBaseCount;
        gr.observedReadCount = observedReadCount;
        geneRecordSet.add(gr);
      }
      getAlignmentMap().clear();
    }

    return geneRecordSet;
  }

  private void readGeneClassInformation(InputStream fastaInputStream) throws IOException {

    SequenceReader sequenceReader = FastaReader.getReader(fastaInputStream);
    assert sequenceReader != null;

    Sequence seq = null;
    while ((seq = sequenceReader.nextSequence(Alphabet.DNA())) != null) {
      putGeneSequence(seq.getName(), seq);
      String desc = seq.getDesc();
      //LOG.info("got geneDesc: " + desc);
      String[] toks = desc.split(";");
      for (String tok : toks) {
        if (tok.startsWith("dg=")) {
          String geneGroup = tok.substring(3);
          //LOG.info("* geneGroup: " + geneGroup);
          putGeneGroup(seq.getName(), geneGroup);
        }
        if (tok.startsWith("geneID=")) {
          String proteinID = tok.substring(7);
          //LOG.info("* proteinID: " + proteinID);
          putGeneName(seq.getName(), proteinID);
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

      String geneID = toks[0];
      String geneName = toks[1];
      String geneGroup = toks[2];

      //LOG.info("got geneInfo: " + line);
      //LOG.info("* geneName: " + geneName);
      //LOG.info("* geneGroup: " + geneGroup);

      if (getGeneName(geneID) == null)
        putGeneName(geneID, geneName);
      else
        putGeneName(geneID, getGeneName(geneID) + ", " + geneName);

      if (getGeneGroup(geneID) == null)
        putGeneGroup(geneID, geneGroup);
      else
        putGeneGroup(geneID, getGeneGroup(geneID) + ", " + geneGroup);
    }
    bf.close();
  }
}
