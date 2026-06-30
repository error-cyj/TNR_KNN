#ifndef TNR_2HOP_H
#define TNR_2HOP_H

#include <stdio.h>
#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <unordered_map>
#include <cmath>
#include <chrono>
#include <set>
#include "unordered_set"
#include "ET.h"
#include <filesystem>
#include "Fast2Hop.h"
#include <climits>
using namespace std::chrono;
using namespace std;

class TNR_2hop : public ET {
public:

    int up_tn = 0;
    int up_ops = 0;
    double gloab = 0;
    double local = 0;
    int dij =0;
    int fw=0;
    int edg=0;
    int nearr=0;
    int nums_local = 0;
    int MAXk;
    const int LIMIT_K = 100;
    const int LIMIT_ACCESS = 564;
    int y;
    int* global_top_k_dists = nullptr;
    int* global_top_k_ids = nullptr;
    int* global_batch_targets = nullptr;
    int* global_exact_dists = nullptr;

    struct CellConn {
        int cell_id;
        int lb;
        int col;
    };
    vector<vector<CellConn>> tn_connected_cells;

    struct HeapNode {
        int dist;
        int u;
    };
    HeapNode* global_dijkstra_heap = nullptr;
    int D_k = INT_MAX;
    int current_k_count = 0;

    vector<vector<pair<int, int>>> local_ancestors;
    struct Cell {
        int cell_id;
        vector<int> tns;
        unordered_map<int, int> tn_local_idx;
        vector<int> tn_lower_bounds;
        vector<int> local_nodes;
        vector<int> tn_to_node_dists;
        vector<int> active_object_rows;
    };
    vector<vector<int>> cell_local_cells;
    size_t max_pq_capacity = 0;
    vector<Cell> global_cells;

    vector<int> node_to_cell_id;
    vector<int> node_to_cell_row;
    vector<int> node_to_active_idx;
    vector<int> global_cell_active_counts;

    vector<vector<pair<int, int>>> bag_nodes_up;
    vector<vector<pair<int, int>>> bag_nodes_down;

    struct MergeNode {
        int value;
        int original_int;
        int list_index;
        int element_index;
        bool operator>(const MergeNode& other) const {
            return value > other.value;
        }
    };

    uint32_t current_version = 0;
    uint32_t cur_src_ver = 0;
    vector<int> dist;
    vector<uint32_t> version;

    struct NodeState {
        uint32_t settled_ver = 0;
        uint32_t batch_ver = 0;
        uint32_t push_ver = 0;
        int best_tn_dist = INT_MAX;
    };
    vector<int> top_k_pos;
    vector<NodeState> node_states;
    vector<int> current_tns;
    vector<vector<pair<int, int>>> tn_network;
    struct AccessEntry {
        int dist;
        int tn_id;
    };
    struct node_tnr {
        bool isTN = false;
        vector<AccessEntry> AccessNode;
        int down_obj_count = 0;
    };
    vector<node_tnr> tnr_index;
    Fast2Hop f2h;

    vector<int> temp_dist_acc;
    vector<int> touched_acc;
    vector<int> temp_dist_anc;
    vector<int> touched_anc;
    vector<bool> global_is_tn;

    TNR_2hop(int k , int x) : ET() {
        MAXk = k;
        y = x;
    }

    ~TNR_2hop() {
        if (global_top_k_dists) {
            delete[] global_top_k_dists;
            delete[] global_top_k_ids;
            delete[] global_batch_targets;
            delete[] global_exact_dists;
            delete[] global_dijkstra_heap;
        }
    }


    void build_TN_Network();
    void preprocess_KNN(string s , string o);
    void contraction_near_find_down_new(int src);
    void contraction_near_find_up_new(int src);
    void Add_Complete(int src);
    void Delete_Complete(int src);
    void build_super_access_cells();
    void save_tnr_index(const string& filename);
    bool load_tnr_index(const string& filename);
    void quickQueryKNN(int src, int k, vector<pair<int, int>> &ans);
};

#endif