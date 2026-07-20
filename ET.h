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
    int MAXk;

    int core_start_rank;
    vector<int> core_tns;
    vector<bool> global_is_tn;
    vector<vector<pair<int, int>>> core_network;

    vector<int> node_order;
    vector<int> order_find;
    vector<int> affected_p_tns;

    vector<int> tree_parent;
    vector<vector<int>> tree_children;


    vector<vector<int>> bag_nodes;
    vector<vector<int>> bag_weights;

    struct Cell {
        int primary_tn = 0;
        vector<int> access_tns;
        vector<int> local_nodes;
        vector<int> tn_to_node_dists;
        bool is_valid = false;
    };

    vector<Cell> tn_cells;
    vector<int> node_to_primary_tn;
    vector<bool> global_is_primary_tn;



    vector<int> static_tree_parent;
    vector<vector<int>> static_tree_children;
    vector<int> local_down_obj_count;

    void initialize_static_tree_structures();

    ET() : Graph() {}

    void load_raw_graph(const string& filename);
    void build_et(int k);
    void generate_full_greedy_order();
    void build_tree_and_truncate_by_object();
    void compute_distances_and_core_network();

    struct EdgeSupport {
        int min_weight = 1000000000;
        vector<int> weights;

        inline bool add(int w) {
            weights.push_back(w);
            if (w < min_weight) {
                min_weight = w;
                return true;
            }
            return false;
        }

        inline bool remove(int w) {
            for (size_t i = 0; i < weights.size(); ++i) {
                if (weights[i] == w) {
                    weights[i] = weights.back();
                    weights.pop_back();
                    if (w == min_weight) {
                        min_weight = 1000000000;
                        for (int x : weights) if (x < min_weight) min_weight = x;
                        return true;
                    }
                    return false;
                }
            }
            return false;
        }
    };

    vector<int> tn_infection_count;
    vector<vector<pair<int, int>>> original_upward_edges;
    vector<unordered_map<int, EdgeSupport>> dynamic_shortcuts;

    unordered_set<int> Dynamic_Topology_Update(const vector<int>& flip_nodes);
    void Assemble_Global_Cells();

};

#endif // PURE_BCT_H