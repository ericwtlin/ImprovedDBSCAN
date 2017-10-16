/**
 * Created by eric-lin on 2017/10/9.
 */

import beans.TrajectoryPoint;
import org.apache.commons.math3.ml.clustering.*;
import java.util.ArrayList;
import java.util.List;
import java.util.HashSet;
import java.io.*;
import com.alibaba.fastjson.JSON;

public class Main {
    public static void main(String[] args) throws FileNotFoundException, IOException {


        /*
        File[] files = getFiles("./data");

        MyDBSCANClusterer<DoublePoint> dbscan = new MyDBSCANClusterer<DoublePoint>(.05, 50);
        List<DoublePoint> points = getGPS(files);
        List<Cluster<DoublePoint>> cluster = dbscan.cluster(points);

        for(Cluster<DoublePoint> c: cluster){
            System.out.println(c.getPoints().get(0));
        }
        */
    }

    private static File[] getFiles(String args) {
        return new File(args).listFiles();
    }

    private static List<DoublePoint> getGPS(File[] files) throws FileNotFoundException, IOException {

        List<DoublePoint> points = new ArrayList<DoublePoint>();
        for (File f : files) {
            BufferedReader in = new BufferedReader(new FileReader(f));
            String line;

            while ((line = in.readLine()) != null) {
                try {
                    TrajectoryPoint p = JSON.parseObject(line, TrajectoryPoint.class);
                    double[] d = new double[2];
                    d[0] = p.getLongitude();
                    d[1] = p.getLatitude();
                    points.add(new DoublePoint(d));
                } catch (ArrayIndexOutOfBoundsException e) {
                } catch(NumberFormatException e){
                }
            }
        }
        return points;
    }


}
