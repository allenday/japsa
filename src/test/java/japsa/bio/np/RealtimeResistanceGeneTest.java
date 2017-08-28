package japsa.bio.np;

import japsa.bio.detection.DetectedGeneRecord;
import japsa.seq.FastaReader;
import japsa.seq.SequenceOutputStream;
import org.apache.commons.io.FileUtils;
import org.junit.Test;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.*;
import java.util.Date;
import java.util.Set;

import static junit.framework.TestCase.assertTrue;

public class RealtimeResistanceGeneTest {
  private static final Logger LOG = LoggerFactory.getLogger(RealtimeResistanceGeneTest.class);

  //	public RealtimeResistanceGene(int read, int time, String output, String resDB, String tmp) throws IOException{
  @Test
  public void testTyping() throws Exception {
    int readNumber = 0;
    int timeNumber = 1;
    Double scoreThreshold = 1D;
    String resDir  = "src/test/resources/resFinder";
    String recordPrefix = "JUNIT";

    File outFile = File.createTempFile("AnalysisResult_",".json");
    LOG.info("output tmp file = "+outFile.getAbsolutePath());

    File fastaFile = new File("src/test/resources/resFinder/DB.fasta");
    assertTrue(fastaFile.exists());
    InputStream fastaInputStream0 = new FileInputStream(fastaFile);

    File resDBFile = new File("src/test/resources/resFinder/geneList");
    assertTrue(resDBFile.exists());
    InputStream resDBInputStream0 = new FileInputStream(resDBFile);

    File bamFile = new File("src/test/resources/jsa361.sam");
    assertTrue(bamFile.exists());
    InputStream bamInputStream0 = new FileInputStream(bamFile);

    BufferedReader bamReader0 = new BufferedReader(new InputStreamReader(bamInputStream0));
    assertTrue(bamReader0.ready());

    ByteArrayOutputStream outStream = new ByteArrayOutputStream();
    RealtimeResistanceGene rg = new RealtimeResistanceGene(readNumber, timeNumber, outStream, resDBInputStream0, fastaInputStream0, recordPrefix);
    RealtimeResistanceGene.JSON = true;
    RealtimeResistanceGene.scoreThreshold = scoreThreshold;

    rg.typing(bamInputStream0);
    try {
      Thread.sleep(10000);
    } catch (InterruptedException e) {
      e.printStackTrace();
    }

    String jsonLine = null;

    BufferedReader br = new BufferedReader(new InputStreamReader(new ByteArrayInputStream(outStream.toByteArray())));
    jsonLine = br.readLine();
    LOG.info("got JSON 1: "+jsonLine);
    assertTrue(jsonLine.indexOf("currentReadCount\":1}") > 0);

    jsonLine = br.readLine();
    LOG.info("got JSON 2: "+jsonLine);
    LOG.info(""+jsonLine.indexOf("currentReadCount"));
    assertTrue(jsonLine.indexOf("currentBaseCount\":842") > 0);
    assertTrue(jsonLine.indexOf("JSA_361") > 0);

    Set<DetectedGeneRecord> geneRecords = rg.getGenes();
    for (DetectedGeneRecord gr : geneRecords) {
      LOG.info("PULL1: "+gr.asJsonString());
    }
    geneRecords = rg.getGenes();
    for (DetectedGeneRecord gr : geneRecords) {
      LOG.info("PULL2: "+gr.asJsonString());
    }

  }
}