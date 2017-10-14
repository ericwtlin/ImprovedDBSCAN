/**
 * Created by eric-lin on 17-10-14.
 */
import java.util.HashMap;


public class UnionFind {
    HashMap<Integer, Integer> group = new HashMap<Integer, Integer>();
	HashMap<Integer, Integer> rank = new HashMap<Integer, Integer>();

	public UnionFind(int[] initialGroups){
        for (int i = 0; i < initialGroups.length; i ++){
			group.put(i, initialGroups[i]);
            if (rank.containsKey(i) == false)
				rank.put(i, 1);
			rank.put(initialGroups[i], 2);		//the rank of origin cluster root is 2
		}
		group.put(-1, -1);		//-1 stands for noisy data
	}

	public UnionFind(int n){
		for (int i = 0; i < n; i ++){
			group.put(i, i);
			rank.put(i, 1);
		}
	}

	/*
	 * find运算(加速)
	 * 从元素e相应的结点走到树根处，找出所在集合的名字
	 */
	private int find(int p){
        //路径压缩
		if (p != group.get(p)){
			group.put(p, find(group.get(p)));
		}
		return group.get(p);

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
			if (rank.get(rootA) <= rank.get(rootB)){
				group.put(rootA, rootB);
				if (rank.get(rootA) == rank.get(rootB)){
					rank.put(rootB, rank.get(rootB) + 1);
				}
			}else{
				group.put(rootB, rootA);
			}
		}
	}

	public static void main(String[] args) {
        int[] groups = {2, 1, 2, 3, 7, 5, 6, 7};
		UnionFind fuf = new UnionFind(groups);
		fuf.union(1, 3);
		fuf.union(2, 6);
		fuf.union(3, 5);
		for (int i = 0; i < groups.length; i ++){
			System.out.println(String.format("%d parent: %d", i, fuf.find(groups[i])));
		}
	}
}