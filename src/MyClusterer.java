import distance.EuclideanDistance;
import net.sf.javaml.core.kdtree.KDTree;
import org.apache.commons.math3.ml.clustering.Cluster;
import org.apache.commons.math3.ml.clustering.Clusterable;
import org.la4j.iterator.VectorIterator;
import org.la4j.matrix.sparse.CRSMatrix;

import java.io.*;
import java.nio.Buffer;
import java.util.*;
import org.apache.log4j.*;
import utils.FileOperation;
/**
 *
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

    int corePointNum = -1;
    int mergeTimes = -1;

    static Logger logger;

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

    /**
     *
     * @param X
     * @return int[] result,  the result[i] indicates the group num, which might be sparse. Specially, -1 indicates noise.
     */
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
		HashSet<Integer>[] kns = new HashSet[N];
		HashSet<Integer> hs;

        boolean[] cores = new boolean[N];
        ArrayList<Integer> corePts = new ArrayList<Integer>(N/3);

        Scheduler scheduler = new Scheduler(N);
        int percent;
		for (int i = 0; i < N; i++) {
            //get the nearest neighbors. Neighbors are returned in ascending order of distance to key.
		    percent = scheduler.getNewSchedule();
		    if (percent != -1){
		        logger.info(String.format("Get nearest neighbors %d%%", percent));
            }

			Object[] nns = kdtree.nearest(X[i], Math.max(K, minPts) + 1);

			hs = new HashSet<Integer>();
			for (int j = 1; j <= K; j++) // start from the 2nd nn
				hs.add((Integer) nns[j]);
            kns[i] = hs;

            /**
             * If the distance between current node and the (minPts+1)-th nearest neighbors is larger than Eps,
             * current node is not a core point.
             */
            if (EuclideanDistance(X[i], X[(int) (nns[minPts])]) <= distEps){
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
        //int[][] S = new int[corePts.size()][corePts.size()];    //storing the SNN similarity between core points
		int count;

		scheduler = new Scheduler(corePts.size() - 1);
		for (int i = 0; i < (corePts.size() - 1); i++) {
		    percent = scheduler.getNewSchedule();
		    if (percent != -1){
		        logger.info(String.format("Calc SNN similartiy %d%%", percent));
            }
			for (int j = i + 1; j < corePts.size(); j++) {
				// create a link between i-j only if i is in j's kNN neighborhood
				// and j is in i's kNN neighborhood
				if (kns[corePts.get(i)].contains(corePts.get(j)) && kns[corePts.get(j)].contains(corePts.get(i))) {
					count = countIntersect(kns[i], kns[j]);
					S.set(i, j, count);
					S.set(j, i, count);
                    //S[i][j] = count;
                    //S[j][i] = count;
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
                    //if (S[i][j] >= SNNEps){
                        uf.union(corePts.get(i), corePts.get(j));
                        countExpanded ++;
                    }
                }
            }
        }

        logger.info(String.format("Number of core points: %d", corePts.size()));
        logger.info(String.format("Clusters expanded, totally %d pairs of clusters merged", countExpanded));
        this.corePointNum = corePts.size();
        this.mergeTimes = countExpanded;

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

	public static int countIntersect (HashSet<Integer> h1, HashSet<Integer> h2) {
        HashSet<Integer> a;
        HashSet<Integer> b;
        if (h1.size() <= h2.size()) {
            a = h1;
            b = h2;
        } else {
            a = h2;
            b = h1;
        }
        int count = 0;
        for (Integer e : a) {
            if (b.contains(e)) {
                count++;
            }
        }
        return count;
    }


    public double EuclideanDistance(double[] point1, double[] point2){
        double dist = 0;
        for (int i = 0; i < point1.length; i ++){
            dist += Math.pow(point1[i] - point2[i], 2);
        }
        return Math.sqrt(dist);
    }

    public void outputInList(String groupPath){
        /*
        output in a single file, in which each row indicates a point.
        The last column of a row means the group index.
         */
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

    /**
     * output in a directory which contains files correponding to each gruops
     * In the group file, each row indicates a point
     *
     * In the directory, there is another summary file that contains summary information.
     * @param groupDirPath
     * @throws Exception
     */
	public void outputInGroups(String groupDirPath) throws Exception {

        if (!groupDirPath.endsWith(File.separator)){
            groupDirPath = groupDirPath + File.separator;
        }

        FileOperation.mkDir(groupDirPath);

        HashMap<Integer, ArrayList<Integer>> groupListMap = new HashMap<>();   //key: intial sparse cluster result index, value: point indices
        HashMap<Integer, Integer> groupSizeMap = new HashMap<>(); // key: final cluster index; value: group size
        groupListMap.put(-1, new ArrayList<>());
        for (int i = 0; i < clusterResult.length; i ++){
            if (!groupListMap.containsKey(clusterResult[i])){
                groupListMap.put(clusterResult[i], new ArrayList<>());
            }
            groupListMap.get(clusterResult[i]).add(i);
        }

        StringBuilder sb = null;
        BufferedWriter out = null;
        int groupCount = -1;
        if (groupListMap.containsKey(-1))  // noise
            groupCount ++;

        int maxGroupSize = -1;
        int minGroupSize = Integer.MAX_VALUE;
        int averageGroupSize = 0;       //except noise
        try {
            for (HashMap.Entry entry: groupListMap.entrySet()){
                if ((int)entry.getKey() == -1){
                    out = new BufferedWriter(new FileWriter(new File(String.format("%snoise.txt", groupDirPath))));
                }else{
                    groupCount ++;
                    out = new BufferedWriter(new FileWriter(new File(String.format("%scluster_%d.txt", groupDirPath, groupCount))));
                }

                for (int pointIndex: groupListMap.get(entry.getKey())){
                    sb = new StringBuilder();
                    for (int j = 0; j < X[pointIndex].length; j ++){
                        if (j != 0){
                            sb.append(",");
                        }
                        sb.append(X[pointIndex][j]);
                    }
                    sb.append("\n");
                    out.write(sb.toString());
                }
                out.close();
                if ((int)entry.getKey() == -1){
                    groupSizeMap.put(0, groupListMap.get(entry.getKey()).size());
                }else {
                    groupSizeMap.put(groupCount, groupListMap.get(entry.getKey()).size());
                    averageGroupSize += groupListMap.get(entry.getKey()).size();
                    if ( groupListMap.get(entry.getKey()).size() > maxGroupSize){
                        maxGroupSize =  groupListMap.get(entry.getKey()).size();
                    }
                    if ( groupListMap.get(entry.getKey()).size() < minGroupSize){
                        minGroupSize =  groupListMap.get(entry.getKey()).size();
                    }
                }
            }
            averageGroupSize = averageGroupSize / groupCount;

            out = new BufferedWriter(new FileWriter(new File(String.format("%sa_summary.txt", groupDirPath))));
            out.write(String.format("Core point num: %d\n", this.corePointNum));
            out.write(String.format("Merged times: %d\n", this.mergeTimes));
            out.write(String.format("Final cluster num(except noise): %d\n", groupCount));
            out.write(String.format("Max cluster size: %d; min cluster size: %d; average cluster size: %d\n\n", maxGroupSize, minGroupSize, averageGroupSize));
            out.write(String.format("noise points: %d\n", groupSizeMap.get(0)));
            for (int i = 1; i <= groupCount; i ++){
                out.write(String.format("cluster_%d:%d\n", i, groupSizeMap.get(i)));
            }
            out.close();
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


    public static void main(String[] args) throws Exception {
	    PropertyConfigurator.configure( "log4j.properties" );

        if (args.length <= 0) {
          System.out.println("Command: java Clusterer" +
                          "--K=K nearest neighbors " +
                          "--SNNEps=SNNEps " +
                          "--distEps=distEps " +
                          "--minPts= minPts " +
                          "--candidate_path=candidate path "
          );
          System.exit(1);
        }

        ConfigParser config = new ConfigParser(args);
        config.checkRequiredArgs(new String[] { "K", "SNNEps", "distEps",
            "minPts", "candidate_path"});

        int K = Integer.valueOf(config.getOption("K"));
        int SNNEps = Integer.valueOf(config.getOption("SNNEps"));
        int distEps = Integer.valueOf(config.getOption("distEps"));
        int minPts = Integer.valueOf(config.getOption("minPts"));
        String candidatePath = config.getOption("candidate_path");

        /*
        int K = 500;
        int SNNEps = 50;
        int distEps = 100;
        int minPts = 2000;
        */

        logger = Logger.getLogger(MyClusterer.class);

        BufferedReader in = new BufferedReader(new FileReader(candidatePath));
        String line;
        String[] lineSplit;
        line = in.readLine();
        lineSplit = line.split(" ");
        int N = Integer.valueOf(lineSplit[0]);
        N = 150000;
        int dim = Integer.valueOf(lineSplit[1]);
        double[][] X = new double[N][dim];
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
        MyClusterer cls = new MyClusterer(K, SNNEps, distEps, minPts, X);

        cls.cluster(X);
        //cls.outputInList(String.format("/home/eric-lin/StateGrid/TrajectoriesMining/cluster_results/station_candidate_withtime_1314.N%d_K%d_SNNEps%d_distEps%d_minPts%d.cluster", N, K, SNNEps, distEps, minPts));
        cls.outputInGroups(String.format("/home/eric-lin/workspace/ImprovedDBSCAN/cluster_results/station_candidate_withtime_1314.N%d_K%d_SNNEps%d_distEps%d_minPts%d.clusters", N, K, SNNEps, distEps, minPts));


    }
}
