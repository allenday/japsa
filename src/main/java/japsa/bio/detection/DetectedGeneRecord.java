package japsa.bio.detection;

import com.google.gson.Gson;
import com.google.gson.GsonBuilder;
import com.google.gson.JsonObject;

public class DetectedGeneRecord {
  public String geneID;
  public String geneName;
  public String geneGroup;
  public Long timeStart = null;
  public Long timeNow = null;
  public Integer observedReadCount = 0;
  public Long observedBaseCount = 0L;
  private static final Gson gson = new GsonBuilder().serializeNulls().create();

  public DetectedGeneRecord() {
  }

  public String asJsonString() {
    JsonObject jo = asJsonObject();
    return gson.toJson(jo);
  }

  public JsonObject asJsonObject() {
    JsonObject jo = new JsonObject();
    jo.addProperty("geneID", geneID);
    jo.addProperty("geneName", geneName);
    jo.addProperty("geneGroup", geneGroup);
    jo.addProperty("observedBaseCount", observedBaseCount);
    jo.addProperty("observedReadCount", observedReadCount);
    if (timeStart != null)
      jo.addProperty("timeStart", timeStart);
    if (timeNow != null)
      jo.addProperty("timeNow", timeNow);
    return jo;
  }
}
