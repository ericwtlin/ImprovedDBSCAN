import distance.EuclideanDistance;
import net.sf.javaml.core.kdtree.KDTree;
import org.apache.commons.math3.ml.clustering.Cluster;
import org.apache.commons.math3.ml.clustering.Clusterable;
import org.la4j.iterator.VectorIterator;
import org.la4j.matrix.sparse.CRSMatrix;

import java.io.*;
import java.util.*;
import org.apache.log4j.*;

/**
 *
 * Developed based on SNN of Cássio M. M. Pereira <cassiomartini@gmail.com>
 */
public class MyClusterer {

	/**
	 *
	 *
	 * @param X
	 *            the data matrix, one example per row
	 * @param K
	 *            number of neighbors to form the sparse similarity matrix
	 * @param SNNEps
	 *            the threshold of the number of shared K nearest neighbors,
     *            used for determining whether to merge clusters of two core points
     * @param distEps
     *            used for determining the core points which should have MinPts neighbors under the distance DistEps
	 * @param minPts
	 *            used for determining core points
	 * @return cluster labels for the points; Note: the cluster labels might not range between [0, cluster_num)
	 */
	int K;
    int SNNEps;
    double distEps;
    int minPts;
    double[][] X = null;
    int[] clusterResult = null;

    Logger logger  =  Logger.getLogger(MyClusterer.class);

	public MyClusterer(int K, int SNNEps, double distEps, int minPts, double[][] X){
        this.K = K;
        this.SNNEps = SNNEps;
        this.distEps = distEps;
        this.minPts = minPts;
        this.X = X;
    }

    public int[] cluster(){
	    return this.cluster(this.X);
    }

	public int[] cluster(double[][] X) {
		int N = X.length; // number of points
		int d = X[0].length; // dimensionality

        /*
		if (minPts >= K) {
			throw new RuntimeException(
					"MinPts has to be smaller than K. No sense in a point having more than K neighbors.");
		}
		*/

		// STEP 1 - get a similarity matrix

		// construct the kd-tree for knn queries
		KDTree kdtree = new KDTree(d);
		for (int i = 0; i < N; i++)
			kdtree.insert(X[i], i);

		logger.info(String.format("%d points inserted to the KDTree", N));

		// STEP 2 - sparsify the matrix by keeping only the k most similar neighbors, get core points

		// find the K-neighbors of each point
		HashMap<Integer, HashSet<Integer>> kns = new HashMap<Integer, HashSet<Integer>>();
		HashSet<Integer> hs;

        boolean[] cores = new boolean[N];
        ArrayList<Integer> corePts = new ArrayList<Integer>(N/3);
		for (int i = 0; i < N; i++) {
            //get the nearest neighbors. Neighbors are returned in ascending order of distance to key.
			Object[] nns = kdtree.nearest(X[i], Math.max(K, minPts) + 1);

			hs = new HashSet<Integer>();
			for (int j = 1; j <= K; j++) // start from the 2nd nn
				hs.add((Integer) nns[j]);
            kns.put(i, hs);

            /**
             * If the distance between current node and the (minPts+1)-th nearest neighbors is larger than Eps,
             * current node is not a core point.
             */
            if (EuclideanDistance(X[i], X[(int) (nns[minPts + 1])]) <= distEps){
                cores[i] = true;
                corePts.add(i);
            }
		}
		logger.info(String.format("Number of core points: %d", corePts.size()));

		logger.info("K-neighbors of each point found and core points found.");

		// STEP 3 - construct the shared nearest neighbor graph from the
		// sparsified matrix

		// The sparse matrix S holds in element (i,j) the SNN-similarity between
		// points i and j.
		CRSMatrix S = new CRSMatrix(N, N);
		int count;

		Scheduler scheduler = new Scheduler(N-1);
		int percent;
		for (int i = 0; i < (N - 1); i++) {
		    percent = scheduler.getNewSchedule();
		    if (percent != -1){
		        logger.info(String.format("Calc SNN similartiy %d%%", percent));
            }
			for (int j = i + 1; j < N; j++) {
				// create a link between i-j only if i is in j's kNN neighborhood
				// and j is in i's kNN neighborhood
				if (kns.get(i).contains(j) && kns.get(j).contains(i)) {
					count = countIntersect(kns.get(i), kns.get(j));
					S.set(i, j, count);
					S.set(j, i, count);
				}
			}
		}

		logger.info("SNN similarity figured out.");

		// STEP 4 - form clusters from the core points. If two core pts are
		// within Eps of each other, then place them in the same cluster
        int[] cluster = new int[N]; //labels[i] indicates node i belongs to cluster of labels[i]; labels[i] == -1 indicates node i is a noise;
        Arrays.fill(cluster, -1);
        scheduler = new Scheduler(corePts.size());
        for (int i = 0; i < corePts.size(); i ++){
            percent = scheduler.getNewSchedule();
            if (percent != -1) {
                logger.info(String.format("Forming clusters %d%%", percent));
            }
            findNeighborsInDistEps(corePts.get(i), X, cluster);
        }

        logger.info("Clusters formed.");

        scheduler = new Scheduler(corePts.size());
        UnionFind uf = new UnionFind(cluster);
        int countExpanded = 0;
        for(int i = 0; i < corePts.size() - 1; i ++){
            percent = scheduler.getNewSchedule();
		    if (percent != -1){
		        logger.info(String.format("Expanding clusters %d%%", percent));
            }
            for (int j = i + 1; j < corePts.size(); j ++){
                if (i != j){
                    if (S.get(corePts.get(i), corePts.get(j)) >= SNNEps){
                        uf.union(corePts.get(i), corePts.get(j));
                        countExpanded ++;
                    }
                }
            }
        }

        logger.info(String.format("Clusters expanded, totally %d pairs of clusters merged", countExpanded));
        clusterResult = uf.getGroup();
        HashSet<Integer> groupId = new HashSet<Integer>();
	    for (int idx : clusterResult){
	    	groupId.add(idx);
		}
		logger.info(String.format("Clustering end! Got %d clusters", groupId.size()));
        return uf.getGroup();
	}

	public void findNeighborsInDistEps(int point, double[][] X,
                                                       //KDTree kdtree,
                                                       //HashMap<Integer, ArrayList<Integer>> kns,
                                                       int[] cluster){
        // TO BE REFINED
        double dist = 0;
        for (int i = 0; i < X.length; i ++){
            dist = EuclideanDistance(X[point], X[i]);
            if (dist <= distEps){
                if (cluster[i] == -1 || dist < EuclideanDistance(X[i], X[cluster[i]])){
                    cluster[i] = point;
                }
            }
        }
    }

	public static int countIntersect(HashSet<Integer> h1, HashSet<Integer> h2) {
        Set<Integer> result = new HashSet<Integer>();
        result.addAll(h1);
        result.retainAll(h2);
		return result.size();
	}


    public double EuclideanDistance(double[] point1, double[] point2){
        double dist = 0;
        for (int i = 0; i < point1.length; i ++){
            dist += Math.pow(point1[i] - point2[i], 2);
        }
        return Math.sqrt(dist);
    }

    public void output(String groupPath){
        BufferedWriter out = null, outCenter = null;
        StringBuilder sb = new StringBuilder();
        try {
            out = new BufferedWriter((new FileWriter(groupPath)));

            for (int i = 0; i < clusterResult.length; i ++){
                sb = new StringBuilder();
                for (int j = 0; j < X[0].length; j ++){
                    if (j != 0){
                        sb.append(",");
                    }
                    sb.append(X[i][j]);
                }
                sb.append(clusterResult[i]);
                sb.append("\n");
                out.write(sb.toString());
            }
        }catch (IOException e){
            e.printStackTrace();
        }finally {
            try {
                out.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }


	}


    public static void main(String[] args) throws FileNotFoundException, IOException {
        PropertyConfigurator.configure( "log4j.properties" );



        BufferedReader in = new BufferedReader(new FileReader("/home/eric-lin/StateGrid/TrajectoriesMining/station_candidate_withtime_1314.txt"));
        String line;
        String[] lineSplit;
        line = in.readLine();
        lineSplit = line.split(" ");
        int N = Integer.valueOf(lineSplit[0]);
        N = 50000;
        int dim = Integer.valueOf(lineSplit[1]);
        double[][] X = new double[N][dim];
        double[] point;
        for (int i = 0; i < N; i ++){
            line = in.readLine();
            lineSplit = line.split(",");
            for(int j = 0; j < dim; j ++){
                X[i][j] = Double.valueOf(lineSplit[j]);
            }
        }
        in.close();

        System.out.println("Reading data finished!");

        //MyClusterer cls = new MyClusterer(100, 30, 20, 50, X);
        MyClusterer cls = new MyClusterer(200, 30, 100, 100, X);

        cls.cluster(X);
        cls.output("/home/eric-lin/StateGrid/TrajectoriesMining/station_candidate_withtime_1314.txt.cluster");

    }
}
