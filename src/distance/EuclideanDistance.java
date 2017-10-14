package distance;

/**
 * Created by eric-lin on 17-10-14.
 */
public class EuclideanDistance implements Distance{
    public double distance(double[] point1, double[] point2){
        double dist = 0;
        for (int i = 0; i < point1.length; i ++){
            dist += Math.pow(point1[i] - point2[i], 2);
        }
        return Math.sqrt(dist);
    }
}
