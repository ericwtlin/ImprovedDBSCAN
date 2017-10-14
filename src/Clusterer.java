import distance.EuclideanDistance;
import net.sf.javaml.core.kdtree.KDTree;
import org.apache.commons.math3.ml.clustering.Clusterable;
import org.la4j.iterator.VectorIterator;
import org.la4j.matrix.sparse.CRSMatrix;

import java.util.*;

/**
 *
 * Developed based on SNN of Cássio M. M. Pereira <cassiomartini@gmail.com>
 */
public class Clusterer {

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

	public Clusterer(int K, int SNNEps, double distEps, int minPts){
        this.K = K;
        this.SNNEps = SNNEps;
        this.distEps = distEps;
        this.minPts = minPts;
    }

	public int[] snn(double[][] X) {
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

		// STEP 3 - construct the shared nearest neighbor graph from the
		// sparsified matrix

		// The sparse matrix S holds in element (i,j) the SNN-similarity between
		// points i and j.
		CRSMatrix S = new CRSMatrix(N, N);
		int count;

		for (int i = 0; i < (N - 1); i++) {
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

		// STEP 4 - form clusters from the core points. If two core pts are
		// within Eps of each other, then place them in the same cluster
        int[] cluster = new int[N]; //labels[i] indicates node i belongs to cluster of labels[i]; labels[i] == -1 indicates node i is a noise;
        Arrays.fill(cluster, -1);
        for (int i = 0; i < corePts.size(); i ++){
            findNeighborsInDistEps(corePts.get(i), X, cluster);
        }

        UnionFind uf = new UnionFind(cluster);
        for(int i = 0; i < corePts.size(); i ++){
            for (int j = 0; j < corePts.size(); i ++){
                if (i != j){
                    if (S.get(corePts.get(i), corePts.get(j)) >= SNNEps){
                        uf.union(corePts.get(i), corePts.get(j));
                    }
                }
            }
        }

        return uf.getGroup();
	}

	public void findNeighborsInDistEps(int point, double[][] X,
                                                       //KDTree kdtree,
                                                       //HashMap<Integer, ArrayList<Integer>> kns,
                                                       int[] labels){
        // TO BE REFINED
        double dist = 0;
        for (int i = 0; i < X.length; i ++){
            dist = EuclideanDistance(X[point], X[i]);
            if (dist <= distEps){
                if (labels[i] == -1 || dist < EuclideanDistance(X[i], X[labels[i]])){
                    labels[i] = point;
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

}
