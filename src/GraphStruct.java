import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Random;

public class GraphStruct {

	double totalWeight = 0;
	public double[] bucket;
	ArrayList<String> edgeIds = new ArrayList<String>();
	Hashtable<String, Node> nodeList = new Hashtable<String, Node>();
	Hashtable<String, Edge> edgeList = new Hashtable<String, Edge>();
	List<Edge> sortedEdgeList = new ArrayList<Edge>();
	Hashtable<String, Double> removedEdgeList = new Hashtable<String, Double>();
	Hashtable<Integer, ArrayList<Edge>> indexed_Edges = new Hashtable<Integer, ArrayList<Edge>>();

	Edge[] heads;
	Edge[] tails;

	Random randomGenerator = new Random(555);
	Logger logger;
	double maxWeight;
	double minWeight;
	String[] colors;
	ApplicationProperties ap;

	public GraphStruct(ApplicationProperties in_ap, Logger in_logger) {
		ap = in_ap;
		logger = in_logger;
		heads = new Edge[ap.NUMBER_OF_PARTITIONS];
		tails = new Edge[ap.NUMBER_OF_PARTITIONS];
		initializeGraph();
	}

	public void start() {
		ReportIteration(0, 0, 0, 0, 0, 1 - totalDominantRatio());

		int iteration = 0;
		double EXPANSION_THRESHOLD = bucket[0]
				- (bucket[0] * ap.EXPANSION_LIMIT);
		double CUTSize = 1 - totalDominantRatio();
		int countStability = 0;
		int succAccum = 0;
		while (countStability < ap.COUNT_STABILITY && succAccum < 1) {
			iteration++;
			int countSuc = 0;
			int countTry = 0;
			long start = System.currentTimeMillis();

			for (Node n : nodeList.values()) {
				if (edgeIds.size() == 0) {
					for (String id : edgeList.keySet()) {
						edgeIds.add(id);
					}
				}
				countSuc += n.run(edgeList, edgeIds, bucket,
						EXPANSION_THRESHOLD, ap.EXCHANGE_LIMIT, heads, tails,
						iteration);
				countTry += n.numTry;
			}

			if (countSuc == 0)
				succAccum++;
			double dur = (System.currentTimeMillis() - start) / 1000.0;
			double currentCUTSize = 1 - totalDominantRatio();

			if (Math.abs(CUTSize - currentCUTSize) < ap.TERMINATION_THRESHOLD)
				countStability++;

			CUTSize = currentCUTSize;
			ReportIteration(iteration, countStability, dur, countSuc, countTry,
					currentCUTSize);
		}

		String report = "#Itr: " + iteration + " #CStb: " + countStability
				+ " #CUT: " + (1 - totalDominantRatio()) + " #STDev: "
				+ standardDeviation();
		System.out.println(report);
		logger.logEntry(report);
	}

	private void ReportIteration(int iteration, int countStability, double dur,
			int countSuc, int countTry, double currentCUTSize) {
		logger.logEntry("################");
		String report = "#Itr: " + iteration + " #CStb: " + countStability
				+ " # dur: " + dur + " #NumSuc: " + countSuc + "/" + countTry
				+ " #CUT: " + currentCUTSize;// +
		logger.logEntry(report);
	}

	public String standardDeviation() {
		double[] totalPartitionWeight = getPartitionWeights();
		double average = getTotalWeight() / ap.NUMBER_OF_PARTITIONS;

		double squaredSum = 0;
		for (double d : totalPartitionWeight)
			squaredSum += Math.pow(d - average, 2);
		double MSE = Math.sqrt(squaredSum / ap.NUMBER_OF_PARTITIONS);
		return String.format("%.3f", MSE);
	}

	public double totalDominantRatio() {
		double TDR = 0.0;
		// int count = 0;
		for (Node n : nodeList.values()) {
			// if(!n.isInternal()){
			TDR += n.getDominantRatio();
			// count++;
			// }
		}
		return TDR / nodeList.size();
	}

	public String totalDistance() {
		double[] sum_w = new double[ap.NUMBER_OF_PARTITIONS];

		for (Edge e : edgeList.values()) {
			sum_w[e.color] += e.weight;
		}

		double maxDensity = Double.MIN_VALUE;
		double minDensity = Double.MAX_VALUE;
		for (double density : sum_w) {
			if (density > maxDensity)
				maxDensity = density;
			if (density < minDensity)
				minDensity = density;
		}

		return String.format("%.5f", (maxDensity - minDensity));
	}

	public double Average_Density() {
		double[] total_cluster_color = new double[ap.NUMBER_OF_PARTITIONS];
		double[] total_cluster_edges = new double[ap.NUMBER_OF_PARTITIONS];

		for (Edge e : edgeList.values()) {
			total_cluster_color[e.color] += e.weight;
			total_cluster_edges[e.color]++;
		}

		double result = 0;
		for (int i = 0; i < ap.NUMBER_OF_PARTITIONS; i++) {
			result += total_cluster_color[i] / total_cluster_edges[i];
		}
		result = result / ap.NUMBER_OF_PARTITIONS;
		return result;
	}

	private double[] getPartitionWeights() {
		double[] partitionWeights = new double[ap.NUMBER_OF_PARTITIONS];
		for (Edge e : edgeList.values()) {
			partitionWeights[e.color] += e.weight;
		}
		return partitionWeights;
	}

	private double getTotalWeight() {
		if (totalWeight == 0) {
			for (Edge e : edgeList.values()) {
				totalWeight += e.weight;
			}
		}
		return totalWeight;
	}

	public void writeEdgeList(String file, double threshold) {
		try {
			java.io.BufferedWriter rw = new BufferedWriter(new FileWriter(
					new File(file)));
			rw.write("Source\tTarget\tColor\tweight\tLabel\n");
			// Undirected Graph
			for (String edgeID : edgeList.keySet()) {
				Edge edge = edgeList.get(edgeID);
				if (edge.weight > threshold)
					rw.write(edgeID.split("-")[0] + "\t" + edgeID.split("-")[1]
							+ "\t" + edge.color + "\t" + edge.weight + "\t"
							+ edge.src.getId() + edge.dst.getId() + "\n");
			}
			rw.flush();
			rw.close();
		} catch (Exception ex) {
			System.out.println(ex.getMessage());
		}
	}

	private void initializeGraph() {

		double total_Weights = 0.0;
		double total_Weights_Selected = 0.0;
		try {
			double minWT = findMinWeightThreshold();
			BufferedReader graphReader = new BufferedReader(new FileReader(ap.GRAPH));
			String line = "";
			int uniqueEdgeId = 0;
			while ((line = graphReader.readLine()) != null) {
				double weight = Double.parseDouble(line.split(ap.SECTION_MARKER)[1]);
				String key = line.split(ap.SECTION_MARKER)[0];
				if (!key.contains("-")) {
					System.out.println(line);
				}
				String parts[] = key.split("-");

				total_Weights += weight;
				//if (weight > ap.MINIMUM_WEIGHT_THRESHOLD) {
				if (weight > minWT) {
					Node src = null;
					if (nodeList.containsKey(parts[0])) {
						src = nodeList.get(parts[0]);
					} else {
						src = new Node(parts[0], ap.NUMBER_OF_PARTITIONS);
						nodeList.put(parts[0], src);
					}

					Node dst = null;
					if (nodeList.containsKey(parts[1])) {
						dst = nodeList.get(parts[1]);
					} else {
						dst = new Node(parts[1], ap.NUMBER_OF_PARTITIONS);
						nodeList.put(parts[1], dst);
					}

					Edge e = new Edge(src, dst, weight, uniqueEdgeId);
					uniqueEdgeId++;
					total_Weights_Selected += e.weight;
					src.addEdge(e);
					dst.addEdge(e);

					edgeList.put(key, e);
				} else {
					removedEdgeList.put(key, weight);
				}
			}

			graphReader.close();

			// Fill EdgeIds list
			if (edgeIds.size() == 0) {
				for (String id : edgeList.keySet()) {
					edgeIds.add(id);
				}
			}

			// Sort edges in all vertices
			// This is very important for the functionality of the new algorithm
			// because we select the first non dominant therefore edges must be
			// sorted to return the larges-non-domiant
			for (Node n : nodeList.values()) {
				Collections.sort(n.edges, new Comparator<Edge>() {
					@Override
					public int compare(Edge e1, Edge e2) {
						if (e1.weight < e2.weight)
							return 1;
						else if (e1.weight > e2.weight)
							return -1;
						else
							return 0;
					}
				});
			}

			if (ap.COLORING.equals("UNIFORM_RANDOM")) {
				color_init_uniform_random(total_Weights_Selected);
			} else if (ap.COLORING.equals("SIMPLE_SORTED")) {
				color_init_simple_sorted(total_Weights_Selected);
			} else if (ap.COLORING.equals("WEIGHTED_BFS")) {
				color_init_weighted_BFS(total_Weights_Selected);
			} else if (ap.COLORING.equals("HIERARCHICAL")) {
				color_init_hierarchical(total_Weights_Selected);
			}

			logger.logEntry("------------- Graph info -------------");
			logger.logEntry("Num Total Nodes: " + nodeList.size());
			logger.logEntry("Num Total Edges: " + edgeList.size());
			logger.logEntry("Num Average Deg: " + averageDegree());
			logger.logEntry("Total Weight: " + total_Weights);
			logger.logEntry("Total Selected Weight: " + total_Weights_Selected);
			logger.logEntry(String.format("%.2f", ap.MINIMUM_WEIGHT_THRESHOLD)
					+ " "
					+ String.format("%.2f", total_Weights_Selected
							/ total_Weights) + " " + edgeList.size());

			// Maximum and Minimum weights
			calcMaxMin();
			logger.logEntry("MaxWeight: " + maxWeight);
			logger.logEntry("MinWeight: " + minWeight);

		} catch (Exception ex) {
			System.out.println("#ERR - Graph - Initialization: "
					+ ex.getMessage());
			ex.printStackTrace();
		}
	}

	private double findMinWeightThreshold() {
		double minWT = 0.0;
		try {
			BufferedReader graphReader = new BufferedReader(new FileReader(ap.GRAPH));
			String line = "";
			List<Double> weights = new ArrayList<Double>();
			double totalW = 0.0;
			while ((line = graphReader.readLine()) != null) {
				double weight = Double.parseDouble(line.split(ap.SECTION_MARKER)[1]);
				weights.add(weight);
				totalW += weight;
			}
			
			Collections.sort(weights);
			Collections.reverse(weights);
			double sum = 0.0;
			for(double w : weights){
				sum+= w;
				if(sum > totalW * ap.MINIMUM_WEIGHT_THRESHOLD){
					minWT = w;
					break;
				}
			}
			graphReader.close();
			}catch(Exception ex){
				ex.printStackTrace();
			}
			return minWT;
	}

	/**
	 * Given A costume sorted set of edgeIds this method colors the edges
	 * accordingly.
	 * 
	 * @param costumOrderedIds
	 * @param total_Weights_Selected
	 */
	public void colorEdges(List<String> costumOrderedIds,
			double total_Weights_Selected) {
		double[] dist = new double[ap.NUMBER_OF_PARTITIONS];
		double each = total_Weights_Selected / ap.NUMBER_OF_PARTITIONS;
		double consumed = 0.0;
		int color = 0;
		Edge prev = null;
		for (String edgeId : costumOrderedIds) {
			if (consumed > each && color < ap.NUMBER_OF_PARTITIONS) {
				consumed = 0;
				color++;
			}
			if (prev != null && prev.color != color) {
				tails[prev.color] = prev;
			}

			Edge e = edgeList.get(edgeId);
			if (heads[color] == null) {
				heads[color] = e;
			} else {
				e.parent = prev;
				prev.child = e;
			}
			prev = e;

			consumed += e.weight;
			dist[color] += e.weight;
			e.setColor(color);
			e.src.localBucket[color] += e.weight;
			e.dst.localBucket[color] += e.weight;
		}
		tails[prev.color] = prev;

		// Just for test the correctness of Heads and Tails
		int coutFirst = 0;
		for (int i = 0; i < ap.NUMBER_OF_PARTITIONS; i++) {
			// System.out.println(coutFirst);
			Edge first = heads[i];
			coutFirst++;
			if (first == null)
				System.out.println("Error in Color Edges in Graph Struct.");
			while (first.child != null) {
				coutFirst++;
				first = first.child;
			}
		}
		System.out.println("CountFirst: " + coutFirst);

		int coutLast = 0;
		for (int i = 0; i < ap.NUMBER_OF_PARTITIONS; i++) {
			Edge last = tails[i];
			coutLast++;
			while (last.parent != null) {
				coutLast++;
				last = last.parent;
			}
		}
		System.out.println("CountLast: " + coutLast);

		bucket = dist;

		logger.logEntry("--------------------------------------");
		logger.logEntry("LOG: Colors Distribution:");
		for (int k = 0; k < dist.length; k++) {
			logger.logEntry("LOG: Color " + k + "\t=\t" + dist[k]);
		}
		logger.logEntry("LOG: TotDist: \t=\t" + totalDistance());
		logger.logEntry("--------------------------------------");
	}

	/**
	 * Just color edges by starting from the first in the list and continue down
	 * to the end of the list.
	 * 
	 * @param total_Weights_Selected
	 */
	private void color_init_uniform_random(double total_Weights_Selected) {
		// Uniform Random Color Assignment
		Collections.shuffle(edgeIds);
		colorEdges(edgeIds, total_Weights_Selected);
	}

	/**
	 * Sort edges according to weight and send to #colorEdges().
	 * 
	 * @param total_Weights_Selected
	 */
	private void color_init_simple_sorted(double total_Weights_Selected) {
		// Simple Sorted Color Assignment
		List<Edge> sortedEdgeList = new ArrayList<Edge>();
		for (Edge e : edgeList.values())
			sortedEdgeList.add(e);
		Collections.sort(sortedEdgeList, new Comparator<Edge>() {
			@Override
			public int compare(Edge e1, Edge e2) {
				if (e1.weight == e2.weight)
					return 0;
				else if (e1.weight < e2.weight)
					return 1;
				else
					return -1;
			}
		});

		List<String> sortedIds = new ArrayList<String>();
		for (Edge e : sortedEdgeList) {
			String id = e.src.getId() + "-" + e.dst.getId();
			if (!edgeIds.contains(id))
				id = e.dst.getId() + "-" + e.src.getId();
			sortedIds.add(id);
		}
		colorEdges(sortedIds, total_Weights_Selected);
	}

	/**
	 * Hierarchically color edges by finding the maximum and going downwards and
	 * repeat until the whole list is exhausted.
	 * 
	 * @param total_Weights_Selected
	 */
	private void color_init_hierarchical(double total_Weights_Selected) {

		Edge source = null;
		double maxWeight = Double.MIN_VALUE;
		for (Edge e : edgeList.values()) {
			if (e.color == -1 && !e.marked && e.weight > maxWeight) {
				maxWeight = e.weight;
				source = e;
			}
		}

		double[] dist = new double[ap.NUMBER_OF_PARTITIONS];
		double each = total_Weights_Selected / ap.NUMBER_OF_PARTITIONS;
		double consumed = 0.0;
		int color = 0;

		while (source != null && color < ap.NUMBER_OF_PARTITIONS) {

			for (Edge e : edgeList.values()) {
				e.marked = false;
				e.dist = -1;
			}
			hierarchical_BFS(this, source);

			// Coloring
			// -----------------------------------------------------------
			List<Edge> sortedEdgeList = new ArrayList<Edge>();
			for (Edge e : edgeList.values())
				if (e.marked)
					sortedEdgeList.add(e);

			Collections.sort(sortedEdgeList, new Comparator<Edge>() {
				@Override
				public int compare(Edge e1, Edge e2) {
					if (e1.dist == e2.dist)
						return 0;
					else if (e1.dist > e2.dist)
						return 1;
					else
						return -1;
				}
			});

			List<String> sortedIds = new ArrayList<String>();
			for (Edge e : sortedEdgeList) {
				if (edgeList.containsKey(e.src.getId() + "\t" + e.dst.getId())) {
					sortedIds.add(e.src.getId() + "\t" + e.dst.getId());
				} else {
					sortedIds.add(e.dst.getId() + "\t" + e.src.getId());
				}
			}
			for (String id : sortedIds) {
				if (consumed < each) {
					edgeList.get(id).color = color;
					consumed += edgeList.get(id).weight;
					dist[color] += edgeList.get(id).weight;
				} else {
					consumed = 0;
					break;
				}
			}

			// -------------------------------------------------
			color++;
			source = null;
			maxWeight = Double.MIN_VALUE;
			for (Edge e : edgeList.values()) {
				if (e.color == -1 && !e.marked && e.weight > maxWeight) {
					maxWeight = e.weight;
					source = e;
				}
			}
		}

		for (Edge e : edgeList.values()) {
			if (e.color == -1) {
				e.color = ap.NUMBER_OF_PARTITIONS - 1;
			}
		}

		dist = new double[ap.NUMBER_OF_PARTITIONS];
		for (Edge e : edgeList.values()) {
			dist[e.color] += e.weight;
		}

		bucket = dist;

		logger.logEntry("--------------------------------------");
		logger.logEntry("LOG: Colors Distribution:");
		for (int k = 0; k < dist.length; k++) {
			logger.logEntry("LOG: Color " + k + "\t=\t" + bucket[k]);
		}
		logger.logEntry("LOG: TotDist: \t=\t" + totalDistance());
		logger.logEntry("--------------------------------------");
	}

	private void hierarchical_BFS(GraphStruct G, Edge s) {
		// To Do Complete the hierarchical coloring in a way that different
		// subgraphs get colored differently.
		Queue<Edge> q = new LinkedList<Edge>();
		s.dist = 0;
		s.marked = true;
		q.add(s);

		while (!q.isEmpty()) {
			Edge v = q.poll();
			for (Edge w : v.src.edges) {
				if (w.color == -1 && !w.marked && w.weight < v.weight) {
					w.dist = v.dist + Math.abs(w.weight - v.weight);
					w.marked = true;
					q.add(w);
				}
			}
			for (Edge w : v.dst.edges) {
				if (w.color == -1 && !w.marked && w.weight < v.weight) {
					w.dist = v.dist + Math.abs(w.weight - v.weight);
					w.marked = true;
					q.add(w);
				}
			}
		}
	}

	/**
	 * Find edge with maximum weight as SOURCE. Apply Weighted-BFS to the graph
	 * and extract edge-list sorted based on their distance to the SOURCE. Send
	 * Sorted-Edge-List to color_edges for coloring.
	 * 
	 * @param total_Weights_Selected
	 */
	private void color_init_weighted_BFS(double total_Weights_Selected) {

		// Weighted BFS Color Assignment
		// Start from the Edge with the Maximum Weight
		Edge source = new Edge(null, null, 0.0, 0);
		double maxWeight = Double.MIN_VALUE;
		for (Edge e : edgeList.values()) {
			if (e.weight > maxWeight) {
				maxWeight = e.weight;
				source = e;
			}
		}

		// Extract BFS order on the graph.
		bfs(this, source);

		// Extract BFS-Sorted List
		List<Edge> sortedEdgeList = new ArrayList<Edge>();
		for (Edge e : edgeList.values())
			sortedEdgeList.add(e);
		Collections.sort(sortedEdgeList, new Comparator<Edge>() {
			@Override
			public int compare(Edge e1, Edge e2) {
				if (e1.dist == e2.dist)
					return 0;
				else if (e1.dist > e2.dist)
					return 1;
				else
					return -1;
			}
		});

		List<String> sortedIds = new ArrayList<String>();
		for (Edge e : sortedEdgeList) {
			String id = e.src.getId() + "-" + e.dst.getId();
			if (!edgeIds.contains(id))
				id = e.dst.getId() + "-" + e.src.getId();
			sortedIds.add(id);
		}

		// Send the BFS-Sorted list for coloring
		colorEdges(sortedIds, total_Weights_Selected);
	}

	private void bfs(GraphStruct G, Edge s) {
		Queue<Edge> q = new LinkedList<Edge>();
		for (Edge e : G.edgeList.values())
			e.dist = Integer.MAX_VALUE;
		s.dist = 0;
		s.marked = true;
		q.add(s);

		while (!q.isEmpty()) {
			Edge v = q.poll();
			for (Edge w : v.src.edges) {
				if (!w.marked) {
					w.dist = v.dist + Math.abs(w.weight - v.weight);
					w.marked = true;
					q.add(w);
				}
			}
			for (Edge w : v.dst.edges) {
				if (!w.marked) {
					w.dist = v.dist + Math.abs(w.weight - v.weight);
					w.marked = true;
					q.add(w);
				}
			}
		}
	}

	private void calcMaxMin() {
		maxWeight = Double.MIN_VALUE;
		minWeight = Double.MAX_VALUE;
		for (Edge e : edgeList.values())
			if (e.weight > maxWeight)
				maxWeight = e.weight;
			else if (e.weight < minWeight)
				minWeight = e.weight;
	}

	public void degreeDistribution(String outPut) {
		try {
			BufferedWriter rw = new BufferedWriter(new FileWriter(new File(
					outPut)));
			for (Node n : nodeList.values()) {
				rw.write(n.degree() + "\n");
			}
			rw.flush();
			rw.close();
		} catch (Exception ex) {
			ex.printStackTrace();
		}
	}

	public void degreeVStotalweightDistribution() {
		for (Node n : nodeList.values()) {
			System.out.println(n.degree() + "\t" + n.getTotalWeight());
		}
	}

	public double averageDegree() {
		return (2.0 * edgeList.size() / nodeList.size());
	}

	public void clusteringCoefficient1(String output) {
		/*
		 * try { BufferedWriter rw = new BufferedWriter(new FileWriter(new File(
		 * output))); for (Node n : nodeList.values()) { rw.write(n.getId() +
		 * "\t" + n.degree() + "\t" + n.clusteringCoefficient() + "\n"); }
		 * rw.flush(); rw.close(); } catch (Exception ex) {
		 * System.out.println(ex.getMessage()); }
		 */
	}

	public double[] getCountOverlap(String vector) {
		// Create mini-graph for each Vector
		ArrayList<String> docEdgeList = new ArrayList<String>();
		String[] parts = vector.split(";");
		// ignore the first element start from 1
		for (int i = 1; i < parts.length - 1; i++) {
			for (int j = i + 1; j < parts.length; j++) {
				String[] src = parts[i].split("\\+|\\-", 2);
				String[] dst = parts[j].split("\\+|\\-", 2);
				String edge = src[0] + (parts[i].contains("+") ? "+" : "-")
						+ "\t" + dst[0] + (parts[j].contains("+") ? "+" : "-");
				docEdgeList.add(edge);
			}
		}

		// Count Overlap between the document and the main-graph
		double[] doc_count_Overlaps = new double[ap.NUMBER_OF_PARTITIONS + 1];
		for (String edge : docEdgeList) {
			if (edgeList.containsKey(edge)) {
				doc_count_Overlaps[edgeList.get(edge).color] += 1.0;
			} else if (removedEdgeList.containsKey(edge)) {
				doc_count_Overlaps[ap.NUMBER_OF_PARTITIONS] += 1.0;
			} else {
				System.out.println("Error: Edge Not Found");
			}
		}
		return doc_count_Overlaps;
	}

	public double[] getUNWeightedOverlap(String[] docEdgeList) {
		// Count Overlap between the document and the main-graph
		double[] doc_weigh_Overlaps = new double[ap.NUMBER_OF_PARTITIONS + 1];
		for (String edge : docEdgeList) {
			if (edgeList.containsKey(edge)) {
				Edge e = edgeList.get(edge);
				doc_weigh_Overlaps[e.color] += 1;// += e.weight;
			} else {
				doc_weigh_Overlaps[ap.NUMBER_OF_PARTITIONS] += 1; // removedEdgeList.get(edge);
			}
		}
		return doc_weigh_Overlaps;
	}

	public double[] getWeightedOverlap(Hashtable<String, Double> docEdgeList) {
		// Count Overlap between the document and the main-graph
		double[] doc_weigh_Overlaps = new double[ap.NUMBER_OF_PARTITIONS + 1];
		for (String edge : docEdgeList.keySet()) {
			if (edgeList.containsKey(edge)) {
				Edge e = edgeList.get(edge);
				// doc_weigh_Overlaps[e.color] += e.weight;
				doc_weigh_Overlaps[e.color] += docEdgeList.get(edge);
				// doc_weigh_Overlaps[e.color] += docEdgeList.get(edge) *
				// e.weight;
			} else if (removedEdgeList.containsKey(edge)) {
				doc_weigh_Overlaps[ap.NUMBER_OF_PARTITIONS] += removedEdgeList
						.get(edge);
			} else {
				System.out.println("Edge: " + edge + "\t Does not exist!!!");
			}
		}
		return doc_weigh_Overlaps;
	}

	public void weightsDistribution(String output) {
		try {
			BufferedWriter rw = new BufferedWriter(new FileWriter(new File(
					output)));
			for (Edge e : edgeList.values()) {
				rw.write(e.weight + "\n");
			}
			rw.flush();
			rw.close();

		} catch (Exception ex) {
			System.out.println(ex.getMessage());
		}
	}

	public String getPartitionEdgeCount() {
		int[] pec = new int[ap.NUMBER_OF_PARTITIONS];
		for (Edge e : edgeList.values()) {
			pec[e.color]++;
		}

		String res = "";
		for (int c : pec)
			res += c + "-";
		return res;
	}

	public void updateEdgeColors(Hashtable<String, String> convertedKies) {
		for (Edge e : edgeList.values()) {
			e.color = Integer.parseInt(convertedKies.get(e.color + ""));
		}
	}

}
