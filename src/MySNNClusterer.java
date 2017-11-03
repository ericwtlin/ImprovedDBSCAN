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


public class MySNNClusterer {
    static Logger logger;
    static int corePointNum = -1;
    static int mergeTimes = -1;
    public MySNNClusterer(){
        ;
	}

    public static int[] snn(double[][] X, int K, double Eps, int MinPts, double expandEps) {
		int N = X.length; // number of points
		int d = X[0].length; // dimensionality

		if (MinPts >= K) {
			throw new RuntimeException(
					"MinPts has to be smaller than K. No sense in a point having more than K neighbors.");
		}

		int[] labels = new int[N];

		// STEP 1 - get a similarity matrix

		// construct the kd-tree for knn queries
		KDTree kdtree = new KDTree(d);

		for (int i = 0; i < N; i++)
			kdtree.insert(X[i], i);

		// STEP 2 - sparsify the matrix by keeping only the k most similar
		// neighbors

		// find the K-neighbors of each point
		HashMap<Integer, HashSet<Integer>> kns = new HashMap<Integer, HashSet<Integer>>();
		HashSet<Integer> hs;

		Scheduler scheduler = new Scheduler(N);
		int percent;
		for (int i = 0; i < N; i++) {
			// we will query for K + 1 nns because the
			// first nn is always the point itself
			percent = scheduler.getNewSchedule();
		    if (percent != -1){
		        logger.info(String.format("Get K nearest neighbors %d%%", percent));
            }
			Object[] nns = kdtree.nearest(X[i], K + 1);

			hs = new HashSet<Integer>();

			for (int j = 1; j < nns.length; j++) // start from the 2nd nn
				hs.add((Integer) nns[j]);

			kns.put(i, hs);
		}

		// STEP 3 - construct the shared nearest neighbor graph from the
		// sparsified matrix

		// The sparse matrix S holds in element (i,j) the SNN-similarity between
		// points i and j.
		CRSMatrix S = new CRSMatrix(N, N);
		int count;

		scheduler = new Scheduler(N-1);
		for (int i = 0; i < (N - 1); i++) {
			percent = scheduler.getNewSchedule();
		    if (percent != -1){
		        logger.info(String.format("construct SNN graph %d%%", percent));
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

		// System.out.println(S.toCSV());

		// STEP 4 - find the SNN density of each point
		double[] snnDens = new double[N]; // should only contain ints though
		VectorIterator vi;
		double snnSim;

		scheduler = new Scheduler(N);
		for (int i = 0; i < N; i++) {
			percent = scheduler.getNewSchedule();
		    if (percent != -1){
		        logger.info(String.format("find SNN density of each point %d%%", percent));
            }
			vi = S.nonZeroIteratorOfRow(i);
			while (vi.hasNext()) {
				snnSim = vi.next();
				if (snnSim >= Eps)
					snnDens[i]++;
			}
		}

		// STEP 5 - find the core points
		// using MinPts, find all points that have SNN density greater than
		// MinPts
		ArrayList<Integer> corePts = new ArrayList<Integer>(N);
		boolean[] cores = new boolean[N]; // initialized to false by default

		for (int i = 0; i < N; i++) {
			if (snnDens[i] >= MinPts) {
				corePts.add(i);
				cores[i] = true;
			}
		}
		corePointNum = corePts.size();

		System.out.println("Core pts list:");
		System.out.println(corePts.toString());

		// System.out.println("similarities for point 0:");
		// vi = S.nonZeroIteratorOfRow(0);
		// while(vi.hasNext()) {
		// System.out.println("sim to 0: " + vi.next());
		// }

		// STEP 6 - form clusters from the core points. If two core pts are
		// within Eps of each other, then place them in the same cluster
		int C = 0;
		HashSet<Integer> visited = new HashSet<Integer>(corePts.size());

		scheduler = new Scheduler(corePts.size());
		for (int i = 0; i < corePts.size(); i++) {
			percent = scheduler.getNewSchedule();
		    if (percent != -1){
		        logger.info(String.format("form clusters from core points %d%%", percent));
            }
			int p = corePts.get(i);
			if (visited.contains(p))
				continue;
			visited.add(p);
			C++;
			labels[p] = C;
			ArrayDeque<Integer> neighCore = findCoreNeighbors(p, corePts, S, Eps);
			expandCluster(labels, neighCore, corePts, C, S, Eps, visited);
		}

		System.out.println("labels after corepts merges:");
		System.out.println(Arrays.toString(labels));

		// STEP 7 & STEP 8
		//
		// All points that are not within a radius of Eps of a core point are discarded (noise);
		//
		// Assign all non-noise, non-core points to their nearest core point

		scheduler = new Scheduler(N);
		for (int i = 0; i < N; i++) {
			percent = scheduler.getNewSchedule();
		    if (percent != -1){
		        logger.info(String.format("Assign non-noise non-core points to nearest core point %d%%", percent));
            }
			boolean notNoise = false;
			double maxSim = Double.MIN_VALUE;
			int bestCore = -1;
			double sim;

			if (cores[i]) // this is a core point
				continue;

			for (int j = 0; j < corePts.size(); j++) {
				int p = corePts.get(j);
				sim = S.get(i, p);
				if (sim >= Eps)
					notNoise = true;
				if (sim > maxSim) {
					maxSim = sim;
					bestCore = p;
				}
			}

			if (notNoise)
				labels[i] = labels[bestCore];
		}

		return labels;
	}

	private static void expandCluster(int[] labels, ArrayDeque<Integer> neighbors, ArrayList<Integer> corePts, int C,
			CRSMatrix S, double Eps, HashSet<Integer> visited) {

		while (neighbors.size() > 0) {
			int p = neighbors.poll();

			if (visited.contains(p))
				continue;

			mergeTimes ++;
			labels[p] = C;
			visited.add(p);

			ArrayDeque<Integer> neigh = findCoreNeighbors(p, corePts, S, Eps);
			neighbors.addAll(neigh);
		}

	}

	private static ArrayDeque<Integer> findCoreNeighbors(final int p, ArrayList<Integer> corePts, CRSMatrix S,
			final double Eps) {
		ArrayDeque<Integer> neighbors = new ArrayDeque<Integer>(corePts.size() / 2);
		int p2;
		for (int i = 0; i < corePts.size(); i++) {
			p2 = corePts.get(i);
			if (p != p2 && S.get(p, p2) >= Eps)
				neighbors.add(p2);
		}
		return neighbors;
	}

	public static int countIntersect(HashSet<Integer> h1, HashSet<Integer> h2) {
		int count = 0;
		for (Integer i : h1)
			if (h2.contains(i))
				count++;
		return count;
	}


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
    /**
     * output in a directory which contains files correponding to each gruops
     * In the group file, each row indicates a point
     *
     * In the directory, there is another summary file that contains summary information.
     * @param groupDirPath
     * @throws Exception
     */
	public void outputInGroups(String groupDirPath, int[] clusterResult, double[][] X) throws Exception {

        if (!groupDirPath.endsWith(File.separator)){
            groupDirPath = groupDirPath + File.separator;
        }

        //FileOperation.mkDir(groupDirPath);

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
        double averageGroupSize = 0;       //except noise
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
			out.write(String.format("merge times: %d\n", this.mergeTimes));
            out.write(String.format("Final cluster num(except noise): %d\n", groupCount));
            out.write(String.format("Max cluster size: %d; min cluster size: %d; average cluster size: %f\n", maxGroupSize, minGroupSize, averageGroupSize));
            out.write(String.format("noise points: %d\n", groupSizeMap.get(0)));

			logger.info(String.format("Core point num: %d;", this.corePointNum));
			logger.info(String.format("merge times: %d;", this.mergeTimes));
            logger.info(String.format("Final cluster num(except noise): %d;", groupCount));
            logger.info(String.format("Max cluster size: %d; min cluster size: %d; average cluster size: %f", maxGroupSize, minGroupSize, averageGroupSize));
            logger.info(String.format("noise points: %d", groupSizeMap.get(0)));

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

        if (args.length <= 0) {
          System.out.println("Command: java Clusterer" +
                          "--K=K nearest neighbors " +
                          "--Eps=SNNEps " +
				          "--expandEps=expandEps " +
                          "--minPts= minPts " +
                          "--candidate_path=candidate path " +
                          "--N=N " +
                          "--force={true,false}"
          );
          System.exit(1);
        }

        ConfigParser config = new ConfigParser(args);
        config.checkRequiredArgs(new String[] { "K", "Eps", "minPts", "candidate_path"});

        int K = Integer.valueOf(config.getOption("K"));
        int minPts = Integer.valueOf(config.getOption("minPts"));
        double Eps = Double.valueOf(config.getOption("Eps"));
		double expandEps = Double.valueOf(config.getOption("expandEps"));
        String candidatePath = config.getOption("candidate_path");
        int NSmall = Integer.valueOf(config.getOption("N", "-1"));
        boolean force = config.getOption("force", "false").toLowerCase().equals("false") ? false : true;

        BufferedReader in = new BufferedReader(new FileReader(candidatePath));
        String line;
        String[] lineSplit;
        line = in.readLine();
        lineSplit = line.split(" ");
        int N;
        if (NSmall == -1){
            N = Integer.valueOf(lineSplit[0]);
        }else{
            N = NSmall;
        }

        String resultDirName = String.format("/home/eric-lin/workspace/ImprovedDBSCAN/snn_results/station_candidate_notime.N%d_K%d_minPts%d.clusters/", N, K, minPts);
        if (!resultDirName.endsWith(File.separator)){
            resultDirName = resultDirName + File.separator;
        }
        FileOperation.mkDir(resultDirName, force);
        String logPath = resultDirName + "cluster.log";
        System.setProperty("logfile.name", logPath);
        PropertyConfigurator.configure( "log4j.properties" );
        logger = Logger.getLogger(MyClusterer.class);

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
        MySNNClusterer cls = null;
        cls = new MySNNClusterer();

        int[] clusterResult = cls.snn(X, K, Eps, minPts, expandEps);
        //cls.outputInList(String.format("/home/eric-lin/StateGrid/TrajectoriesMining/cluster_results/station_candidate_withtime_1314.N%d_K%d_SNNEps%d_distEps%d_minPts%d.cluster", N, K, SNNEps, distEps, minPts));
        cls.outputInGroups(resultDirName, clusterResult, X);

    }
}
