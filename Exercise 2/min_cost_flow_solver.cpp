#pragma once
#include "min_cost_flow_solver.h"
NetworkSimplexSolver::NetworkSimplexSolver(std::string filename) {
	std::ifstream myfile;
	myfile.open(filename);

	if (!myfile.is_open()) {
		std::cout << "FILE DIDNT OPEN!\n";
		exit(0);
	}

	myfile >> num_original_nodes;

	_b = _potential = std::vector<double>(num_original_nodes+1, 0);
	_parent = _edge_to_parent = _depth = std::vector<int>(num_original_nodes+1);

	for (int i = 0; i < num_original_nodes; i++) {
		myfile >> _b[i];
	}

	myfile >> num_original_edges;
	_edges = std::vector<Edge>(num_original_edges+num_original_nodes);
	_flow = std::vector<double>(num_original_edges+num_original_nodes, 0);
	_in_tree = std::vector<bool>(num_original_edges+num_original_nodes, false);

	double max_cost = 0;

	for (int i = 0; i < num_original_edges; i++) {
		myfile >> _edges[i]._start >> _edges[i]._end >> _edges[i]._capacity >> _edges[i]._cost;
		if (_edges[i]._cost > max_cost) max_cost = _edges[i]._cost;
	}

	double artificial_cost = (1 + (num_original_nodes + 1)*max_cost);

	for (int i = 0; i < num_original_nodes; i++) {
		if (_b[i] < 0) {
			_edges[i+num_original_edges] = { num_original_nodes, i, -_b[i], artificial_cost };
		}
		else {
			_edges[i+num_original_edges] = {i, num_original_nodes, _b[i]+1, artificial_cost };
		}
	}

}

void NetworkSimplexSolver::print_instance() {
	std::cout << "Printing instance" << std::endl;
	std::cout << "number of nodes: " << num_original_nodes + 1 << std::endl;
	for (int i = 0; i < num_original_nodes + 1; i++) {
		std::cout << "b[" << i << "] = " << _b[i] << std::endl;
	}
	for (Edge e : _edges) {
		std::cout << "(" << e._start << ", " << e._end << "), u = " << e._capacity << ", c = " << e._cost << std::endl;
	}
}

void NetworkSimplexSolver::initialize_flow() {
	for (int i = 0; i < num_original_nodes; i++) {
		if (_b[i] < 0) {
			_flow[i + num_original_edges] = -_b[i];
			_in_tree[i + num_original_edges] = true;
			_parent[i] = num_original_nodes;
			_edge_to_parent[i] = i + num_original_edges;
		}
		else {
			_flow[i + num_original_edges] = _b[i];
			_in_tree[i + num_original_edges] = true;
			_parent[i] = num_original_nodes;
			_edge_to_parent[i] = i + num_original_edges;
		}
	}
}

void NetworkSimplexSolver::update_potential() {
	std::vector<bool> has_been_updated(num_original_nodes + 1, false);
	has_been_updated[num_original_nodes] = true;
	for (int i = 0; i < num_original_nodes; i++) {
		update_potential_of_node(has_been_updated, i);
	}
}

void NetworkSimplexSolver::update_potential_of_node(std::vector<bool>& has_been_updated, const int node) {
	if (has_been_updated[node] == true) return;
	else {
		//if edge is downward
		if (_edges[_edge_to_parent[node]]._end == node) {
			_potential[node] = _potential[_parent[node]] + _edges[_edge_to_parent[node]]._cost;
		}
		//if edge is updard
		else {
			_potential[node] = _potential[_parent[node]] - _edges[_edge_to_parent[node]]._cost;
		}
	}
}

int main(int argc[], char* argv[]) {
	std::string filename(argv[1]);
	NetworkSimplexSolver solver(filename);
}