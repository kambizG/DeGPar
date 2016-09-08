import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;

public class Node {

	int dist = -1;
	boolean marked = false;

	private String ID;
	ArrayList<Edge> edges = new ArrayList<Edge>();
	double totalWeight = 0;
	public int numTry = 0;
	public int dominantColor = -1;
	Boolean internal = null;
	public double[] localBucket;
	int numClusters = -1;

	public Node(String in_ID, int in_numClusters) {
		ID = in_ID;
		numClusters = in_numClusters;
		localBucket = new double[numClusters];
	}

	public void addEdge(Edge in_edge) {
		edges.add(in_edge);
	}

	public int run(Hashtable<String, Edge> edgeList, List<String> edgeIds,
			double[] bucket, double expansionThreshold, double exchangeRatio,
			Edge[] heads, Edge[] tails, int iteration) {
		numTry = 0;
		int numSuc = 0;
		for (int i = 0; i < 1; i++) {

			if (!this.isInternal()) {
				try {
					int dominantColor = getDominantColor();
					// This always returns an edge because we only examine
					// non-internal edges.
					Edge selfEdge = selectLargestNonDominant(dominantColor);

					if (selfEdge != null
							&& bucket[dominantColor] > expansionThreshold
							&& bucket[selfEdge.color] > expansionThreshold) {
						numTry++;
						double totalWeight = 0.0;
						List<Edge> candidates = new ArrayList<Edge>();
						// ------------------------------------------------------------------------
						// Distant LinkedList Candidate Selection

						Edge e = heads[dominantColor];
						boolean cohesion = true;
						while (true) {
							if (totalWeight < selfEdge.weight) {
								boolean added = false;
								if (totalWeight + e.weight < selfEdge.weight
										* (1 + exchangeRatio)) {
									if (!e.src.isInternal()
											&& !e.dst.isInternal()) {
										if (e.gain(selfEdge.color) > 0) {
											added = true;
											candidates.add(e);
											totalWeight += e.weight;
										}
									}
								}
								if (totalWeight > 0 && cohesion && !added)
									cohesion = false;
								if (e.child != null)
									e = e.child;
								else
									break;
								;
							} else {
								break;
							}
						}
						if (totalWeight > selfEdge.weight * (1 - exchangeRatio)
								&& totalWeight < selfEdge.weight
										* (1 + exchangeRatio)) {
							numSuc++;
							if (totalWeight > selfEdge.weight) {
								bucket[dominantColor] = bucket[dominantColor]
										- Math.abs(totalWeight
												- selfEdge.weight);
								bucket[selfEdge.color] = bucket[selfEdge.color]
										+ Math.abs(totalWeight
												- selfEdge.weight);
							} else {
								bucket[dominantColor] = bucket[dominantColor]
										+ Math.abs(totalWeight
												- selfEdge.weight);
								bucket[selfEdge.color] = bucket[selfEdge.color]
										- Math.abs(totalWeight
												- selfEdge.weight);
							}
							for (Edge cand : candidates) {
								cand.changeColor(heads, tails, selfEdge.color);
								cand.src.dominantColor = -1;
								cand.dst.dominantColor = -1;
								cand.src.internal = null;
								cand.dst.internal = null;

								// update incident vertices
								cand.src.localBucket = null;
								cand.dst.localBucket = null;
							}
							selfEdge.changeColor(heads, tails, dominantColor);
							selfEdge.src.internal = null;
							selfEdge.dst.internal = null;
							selfEdge.src.localBucket = null;
							selfEdge.dst.localBucket = null;
						}
					}
				} catch (Exception ex) {
					ex.printStackTrace();
				}
			}
		}
		return numSuc;
	}

	public double getDominantRatio() {
		Hashtable<Integer, Double> colorCounts = getColorCounts();
		double maxCount = Double.MIN_VALUE;
		for (int color : colorCounts.keySet()) {
			if (colorCounts.get(color) > maxCount) {
				maxCount = colorCounts.get(color);
			}
		}
		double tw = getTotalWeight();
		double val = maxCount / tw;
		if (val > 1) {
			System.out.println(getTotalWeight()
					+ "\t\tOh My God #################################");
		}
		return val;
	}

	private Hashtable<Integer, Double> getColorCounts() {
		Hashtable<Integer, Double> colorCounts = new Hashtable<Integer, Double>();
		for (Edge e : edges) {
			if (colorCounts.containsKey(e.color)) {
				double currentVal = colorCounts.get(e.color);
				colorCounts.put(e.color, currentVal + e.weight);
			} else {
				colorCounts.put(e.color, e.weight);
			}
		}
		return colorCounts;
	}

	/**
	 * Since edges is sorted return first non-dominant edge
	 * 
	 * @return Edge
	 */
	private Edge selectLargestNonDominant(int newColor) {
		int dominant = getDominantColor();
		for (Edge e : edges) {
			if (e.color != dominant) {
				if (e.gain(newColor) > 0)
					return e;
			}
		}
		return null;
	}

	public boolean isInternal() {
		if (internal == null) {
			internal = true;
			for (Edge e : edges) {
				if (e.color != edges.get(0).color) {
					internal = false;
					break;
				}
			}
		}
		return internal;
	}

	public String getId() {
		return ID;
	}

	public String getStrippedId() {
		return ID.substring(0, ID.length() - 1);
	}

	public int getDominantColor() {
		if (dominantColor == -1) {
			Hashtable<Integer, Double> colorCounts = new Hashtable<Integer, Double>();
			dominantColor = Integer.MIN_VALUE;
			double maxCount = Double.MIN_VALUE;
			for (Edge e : edges) {
				if (colorCounts.containsKey(e.color))
					colorCounts.put(e.color, colorCounts.get(e.color)
							+ e.weight);
				else
					colorCounts.put(e.color, e.weight);
				if (colorCounts.get(e.color) > maxCount) {
					maxCount = colorCounts.get(e.color);
					dominantColor = e.color;
				}
			}
		}

		return dominantColor;
	}

	public double getTotalColorWeight(int color) {
		if (localBucket == null) {
			localBucket = new double[numClusters];
			for (Edge e : edges) {
				localBucket[e.color] += e.weight;
			}
		}
		return localBucket[color];
		/*
		 * double count = 0; for(Edge e : edges){ if(e.color == color) count +=
		 * e.weight; } return count;
		 */
	}

	public double getTotalWeight() {
		if (totalWeight == 0) {
			for (Edge e : edges) {
				totalWeight += e.weight;
			}
		}
		return totalWeight;
	}

	public Integer degree() {
		return edges.size();
	}

	public ArrayList<Integer> getUniqieCommunities() {
		ArrayList<Integer> UC = new ArrayList<Integer>();
		for (Edge e : edges) {
			if (!UC.contains(e.color))
				UC.add(e.color);
		}
		return UC;
	}

	public int getNumColors() {
		ArrayList<Integer> colors = new ArrayList<Integer>();
		for (Edge e : edges) {
			if (!colors.contains(e.color))
				colors.add(e.color);
		}
		return colors.size();
	}
}
