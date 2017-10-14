/**
 * Created by eric-lin on 17-10-14.
 */
import java.util.Arrays;
import java.util.ArrayList;
import java.util.HashMap;

public class UnionFind {
    private int[] group;
	private int[] rank;

	public UnionFind(int[] initialGroups){
		group = new int[initialGroups.length];
		rank = new int[initialGroups.length];
		Arrays.fill(rank, 1);
		for (int i = 0; i < initialGroups.length; i ++){
			group[i] = initialGroups[i];
			if (initialGroups[i] != -1)
				rank[initialGroups[i]] = 2;		//the rank of origin cluster root is 2
		}
	}

	public UnionFind(int n){
		group = new int[n];
		rank = new int[n];
        Arrays.fill(rank, 1);
		for (int i = 0; i < n; i ++){
			group[i] = i;
		}
	}

	/*
	 * find运算(加速)
	 * 从元素e相应的结点走到树根处，找出所在集合的名字
	 */
	private int find(int p){
        //路径压缩
        if (p == -1){	  // noise data
			return -1;
		}

		if (p != group[p]){
			group[p] = find(group[p]);
		}
		return group[p];

	}
	/*
	 * 合并两个集合(加速)
	 * 将表示小树的数根改为表示大树的数根的儿子结点
	 */
	public void union(int a, int b){
		int rootA = find(a);
		int rootB = find(b);

		if (rootA == rootB){
			return;
		}else{
			if (rank[rootA] <= rank[rootB]){
				group[rootA] = rootB;
				if (rank[rootA] == rank[rootB]){
					rank[rootB] = rank[rootB] + 1;
				}
			}else{
				group[rootB] = rootA;
			}
		}
	}

	public int[] getGroup(){
		return group;
	}

	public HashMap<Integer, ArrayList<Integer>> getGroups(){
        HashMap<Integer, ArrayList<Integer>> result = new HashMap<>();
        return result;
	}

	public static void main(String[] args) {
        int[] groups = {2, -1, 2, -1, 7, 5, 6, 7};
		UnionFind fuf = new UnionFind(groups);
		fuf.union(1, 3);
		fuf.union(2, 6);
		fuf.union(4, 5);
		for (int i = 0; i < groups.length; i ++){
			System.out.println(String.format("%d parent: %d", i, fuf.find(groups[i])));
		}
	}
}