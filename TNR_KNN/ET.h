#ifndef PURE_BCT_H
#define PURE_BCT_H

#include "Graph.h"
#include <iostream>
#include <vector>
#include <unordered_map>
#include <fstream>
#include <cassert>
#include <algorithm>
#include <climits>

using namespace std;

class ET : public Graph {
public:


    const int SAFE_INF = 1000000000;

    inline __attribute__((always_inline)) int safe_add(int a, int b) {
        if (a >= SAFE_INF || b >= SAFE_INF) return SAFE_INF;
        return a + b;
    }

    int node_nums;
    int edge_nums;

    vector<int> node_order;
    vector<int> order_find;


    vector<int> tree_parent;
    vector<vector<int>> tree_children;

    vector<vector<int>> bag_nodes;
    vector<vector<int>> bag_weights;

    ET() : Graph() {}

    void load_raw_graph(const string& filename);
    void build_et();
    void save_index(const string& filename);
    bool load_index(const string& filename);
private:
    void build_elimination_tree();
    void compute_bag_distances();
    void generate_level_penalized_greedy_order();
};

#endif