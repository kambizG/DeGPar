public class Edge {

	Edge child = null;
	Edge parent = null;

	double dist;
	boolean marked;
	int index = -1;

	Node src;
	Node dst;
	int color = -1;
	double weight;
	int uniqueId;

	String[] colors;

	public Edge(Node in_src, Node in_dst, Double in_weight, int in_id) {
		src = in_src;
		dst = in_dst;
		weight = in_weight;
		uniqueId = in_id;
	}

	// Returns the other end of the edge that it's one end has ID = id.
	public Node getNeibor(String id) {
		if (src.getId().equals(id)) {
			return dst;
		} else {
			return src;
		}
	}

	public double gain(int new_color) {

		double new_gain = this.src.getTotalColorWeight(new_color)
				/ this.src.getTotalWeight()
				+ this.dst.getTotalColorWeight(new_color)
				/ this.dst.getTotalWeight();
		double old_gain = ((this.src.getTotalColorWeight(this.color) - this.weight) / this.src
				.getTotalWeight())
				+ ((this.dst.getTotalColorWeight(this.color) - this.weight) / this.dst
						.getTotalWeight());
		return new_gain - old_gain;

		/*
		 * double total_Old_Color = 0.0; double total_New_Color = 0.0; for(Edge
		 * e: src.edges){ if(e.color == this.color) total_Old_Color += e.weight;
		 * else if(e.color == new_color) total_New_Color += e.weight; } return
		 * (total_New_Color - (total_Old_Color - 2 * this.weight));
		 */
	}

	public Double getValue(int in_color, double in_weight) {
		if (color == in_color) {
			return ((src.getTotalColorWeight(in_color) - weight) * 1.0 / src
					.getTotalWeight())
					+ ((dst.getTotalColorWeight(in_color) - weight) * 1.0 / dst
							.getTotalWeight());
		} else {
			// return (src.getWeightedColorCount(in_color) * 1.0 /
			// src.getTotalWeight()) + (dst.getWeightedColorCount(in_color) *
			// 1.0 / dst.getTotalWeight());
			return ((src.getTotalColorWeight(in_color) + in_weight) * 1.0 / src
					.getTotalWeight())
					+ ((dst.getTotalColorWeight(in_color) + in_weight) * 1.0 / dst
							.getTotalWeight());
		}
	}

	public Double getValue(int in_color, double in_weight, String id) {
		if (id.equals(this.src.getId())) {
			if (this.color == in_color) {
				return (dst.getTotalColorWeight(in_color) - weight) * 1.0
						/ dst.getTotalWeight();
			} else {
				return (dst.getTotalColorWeight(in_color) + in_weight) * 1.0
						/ dst.getTotalWeight();
			}
		} else {
			if (this.color == in_color) {
				return (src.getTotalColorWeight(in_color) - weight) * 1.0
						/ src.getTotalWeight();
			} else {
				return (src.getTotalColorWeight(in_color) + in_weight) * 1.0
						/ src.getTotalWeight();
			}
		}
	}

	public boolean isCrossLink(int color2) {
		for (Edge e : src.edges) {
			String strippedDistId = e.getNeibor(src.getId()).getStrippedId();
			String strippedSrcId = src.getStrippedId();
			if (strippedDistId.equals(strippedSrcId)) {
				if (e.color == color2) {
					return true;
				}
			}
		}

		for (Edge e : dst.edges) {
			String strippedDistId = e.getNeibor(dst.getId()).getStrippedId();
			String strippedDstId = dst.getStrippedId();
			if (strippedDistId.equals(strippedDstId)) {
				if (e.color == color2) {
					return true;
				}
			}
		}
		return false;
	}

	public boolean isEqual(String key) {
		if ((src.getId() + "\t" + dst.getId()).equals(key))
			return true;
		if ((dst.getId() + "\t" + src.getId()).equals(key))
			return true;
		return false;
	}

	public void setColor(int in_color) {
		color = in_color;
	}

	public int getColor() {
		return color;
	}

	public void changeColor(Edge[] heads, Edge[] tails, int new_color) {
		// if head
		if (heads[color] == this) {
			heads[color] = this.child;
			this.child.parent = null;
		} else if (tails[color] == this) {
			// if tail
			this.parent.child = null;
			tails[color] = this.parent;
		} else {
			// if middle
			if (this.parent == null) {
				System.out.println(this.uniqueId + "\t" + this.color);
			}
			this.parent.child = this.child;
			this.child.parent = this.parent;
		}

		// add myself to the new head
		this.parent = tails[new_color];
		this.child = null;
		tails[new_color].child = this;
		tails[new_color] = this;
		this.color = new_color;
	}

}
