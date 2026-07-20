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

    int current_f2h_height = 10;
    bool lazyupdated = false;
    vector<int> node_f2h_level;
    void compute_f2h_levels();

    TNR_2hop(const TNR_2hop&) = delete;
    TNR_2hop& operator=(const TNR_2hop&) = delete;

    int MAXk;
    const int LIMIT_K = 50;
    const int LIMIT_ACCESS = 564;

    struct TopKNode {
        int dist;
        int id;
    };
    TopKNode* global_top_k_heap = nullptr;

    int* global_batch_targets = nullptr;
    int* global_exact_dists = nullptr;

    struct HeapNode {
        int dist;
        int u;
    };
    HeapNode* global_dijkstra_heap = nullptr;
    int D_k = INT_MAX;
    int current_k_count = 0;

    vector<vector<pair<int, int>>> local_ancestors;
    size_t max_pq_capacity = 0;

    vector<int> node_to_cell_row;
    vector<int> node_to_active_idx;


    uint32_t current_version = 0;
    uint32_t cur_src_ver = 0;
    vector<int> dist;
    vector<uint32_t> version;
    vector<int> top_k_pos;

    vector<int> global_to_dense;
    vector<int> dense_to_global;

    struct TNDenseState {
        uint32_t settled_ver = 0;
        uint32_t push_ver = 0;
        int best_tn_dist = INT_MAX;
    };
    vector<TNDenseState> tn_states;


    Fast2Hop f2h;


    struct CellConn {
        int rt_cell_idx;
        int col;
    };


    struct DenseCell {
        int global_p_tn;
        const int* tn_to_node_dists;
        int N;
        vector<pair<int, int>> active_objects;
        vector<int> dense_access_tns;
    };
    vector<DenseCell> rt_cells;
    vector<int> down_obj_count;
    TNR_2hop(int k ) : ET() {
        MAXk = k;
    }

    ~TNR_2hop() {
        if (global_top_k_heap) {
            delete[] global_top_k_heap;
            delete[] global_batch_targets;
            delete[] global_exact_dists;
            delete[] global_dijkstra_heap;
        }
    }

    void init_cell_runtime_states();
    void build_TN_Network();
    void preprocess_KNN(string s , string o);
    void Add_Complete(int src);
    void Delete_Complete(int src);
    void quickQueryKNN(int src, int k, vector<pair<int, int>> &ans);

    void Apply_Index_Updates();
    void build_Fast2Hop();

    void Update_Dirty_Fast2Hop(const unordered_set<int>& dirty_p_tns, const vector<int>& dirty_nodes);
    void Recompute_Dirty_Cells(const unordered_set<int>& dirty_p_tns, const vector<int>& dirty_nodes);
    void Update_TN_Network(const unordered_set<int>& dirty_p_tns);
    void Update_Cell_Runtime_States();


    void Update_Dirty_Core_Network(const unordered_set<int> &dirty_p_tns);
    vector<int> Reassemble_Dirty_Cells(const unordered_set<int>& dirty_p_tns);
    int Find_Best_Split_Point(int p_tn);

    std::queue<int> recycled_dense_ids;

    std::vector<std::vector<CellConn>> tn_connected_cells;

    void Update_Cell_Runtime_States(const unordered_set<int>& dirty_p_tns);
};

#endif