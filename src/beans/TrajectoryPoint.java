/**
 * Created by eric-lin on 2017/10/9.
 */

package beans;
import java.util.*;

public class TrajectoryPoint {
    private String sim = null;
    private int apc = 0;
    private String tflag = null;
    private String ostdesc = null;
    private String utc = null;
    private int dst = 0;
    private String ost = null;
    private Double longitude = 0.0;
    private Double latitude = 0.0;
    private Double distance = 0.0;
    private String suid = null;
    private String oem = null;
    private String head = null;
    private int vflag = 0;
    private Double speed = 0.0;


    public String getSim() {
        return sim;
    }

    public int getApc() {
        return apc;
    }

    public String getTflag() {
        return tflag;
    }

    public String getOstdesc() {
        return ostdesc;
    }

    public String getUtc() {
        return utc;
    }

    public int getDst() {
        return dst;
    }

    public String getOst() {
        return ost;
    }

    public Double getLongitude() {
        return longitude;
    }

    public Double getLatitude() {
        return latitude;
    }

    public Double getDistance() {
        return distance;
    }

    public String getSuid() {
        return suid;
    }

    public String getOem() {
        return oem;
    }

    public String getHead() {
        return head;
    }

    public int getVflag() {
        return vflag;
    }

    public Double getSpeed() {
        return speed;
    }

    public void setSim(String sim) {
        this.sim = sim;
    }

    public void setApc(int apc) {
        this.apc = apc;
    }

    public void setTflag(String tflag) {
        this.tflag = tflag;
    }

    public void setOstdesc(String ostdesc) {
        this.ostdesc = ostdesc;
    }

    public void setUtc(String utc) {
        this.utc = utc;
    }

    public void setDst(int dst) {
        this.dst = dst;
    }

    public void setOst(String ost) {
        this.ost = ost;
    }

    public void setLongitude(Double longitude) {
        this.longitude = longitude;
    }

    public void setLatitude(Double latitude) {
        this.latitude = latitude;
    }

    public void setDistance(Double distance) {
        this.distance = distance;
    }

    public void setSuid(String suid) {
        this.suid = suid;
    }

    public void setOem(String oem) {
        this.oem = oem;
    }

    public void setHead(String head) {
        this.head = head;
    }

    public void setVflag(int vflag) {
        this.vflag = vflag;
    }

    public void setSpeed(Double speed) {
        this.speed = speed;
    }




}
