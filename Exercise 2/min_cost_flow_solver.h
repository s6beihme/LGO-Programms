#include <iostream>
#include <vector>
#include <string>
#include <fstream>
class NetworkSimplexSolver {
public:
	struct Edge {
		int _start, _end;
		double _capacity;
		double _cost;
	};

	//reads in edges and supplies.
	//adds aditional node r and edges (v,r) with capacity _b(v)+1 for all nodes with _b(v)>=0 to r
	//and (r,v) with capacity -_b(v) for all nodes v with _b(v)<0, both with high costs
	NetworkSimplexSolver(std::string filename);

	void print_instance();
	double Solve();
private:

	//initializes flow according to page 61 in Script
	void initialize_flow();

	//calculates the potential in the current tree
	void update_potential();

	void update_potential_of_node(std::vector<bool>& has_been_updated, const int node);

	//finds index of arc not in tree with negative forward or backward reduced cost
	int find_arc_to_augment();

	//augments along edge with index i and updates _potential and _depth
	void augment_along_edge(int i);



	int num_original_nodes;
	int num_original_edges;
	std::vector<Edge> _edges;
	std::vector<double> _flow;
	std::vector<bool> _in_tree;
	std::vector<double> _b;
	std::vector<double> _potential;
	std::vector<int> _parent;
	std::vector<int> _edge_to_parent;
	std::vector<int> _depth;
};

