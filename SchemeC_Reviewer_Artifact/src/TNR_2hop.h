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

    int current_f2h_height = LOCAL_F2H_HEIGHT;

    // ===================== Engineering knobs =====================
    // Start with roomier Cells. Change this single constant for experiments
    // such as 1.5 / 2.0 / 2.5 / 3.0.
    // Main experiment setting: 5.5K is the current low-density sweet spot.
    // Keep this as a compile-time knob so the benchmark log always records
    // exactly which Cell granularity is being tested.
    static constexpr double CELL_CAPACITY_FACTOR = 5.5;

    // Fast2Hop height cap used by the current implementation.
    // This is intentionally exposed for ablation experiments.
    static constexpr int LOCAL_F2H_HEIGHT = 10;

    // Fine-grained timers are guarded by benchmark_collection_enabled. The
    // final experiment driver keeps that flag false for all clean measurements.
    static constexpr bool ENABLE_FINE_GLOBAL_PROFILING = true;

    // Default: Apply_Index_Updates() is lazy. It records stale derived indexes
    // and structural maintenance as pending; queries remain exact through the
    // fallback path until a rebuild. Set true for synchronous rebuilds.
    bool eager_maintenance = false;
    bool maintenance_pending = false;
    bool lazyupdated = false;
    vector<int> node_f2h_level;
    void compute_f2h_levels();

    TNR_2hop(const TNR_2hop&) = delete;
    TNR_2hop& operator=(const TNR_2hop&) = delete;

    int up_tn = 0;
    int up_ops = 0;
    double gloab = 0;
    double local = 0;
    int ace=0;
    int nums_local = 0;
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
    int* global_lca_ranks = nullptr;

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

    // ================================================================
    // Cell-local LCA
    //
    // The Global ET is much larger than any source Cell.  Local Fast2Hop
    // queries only compare vertices that belong to the same Cell, therefore
    // their LCA is guaranteed to stay inside that Cell's ET subtree (including
    // the primary TN as the local root).  We build one compact Euler+RMQ index
    // per primary Cell instead of one global RMQ structure.
    // ================================================================
    struct CellLCAIndex {
        vector<int> nodes;                 // local index -> global node id
        vector<int> euler;                 // Euler sequence of local indices
        vector<int> first;                 // first occurrence in Euler
        vector<int> depth;                 // local depth in the Cell subtree
        vector<int> log2;
        vector<vector<int>> sparse;        // RMQ over Euler, storing local index
    };

    vector<CellLCAIndex> cell_lca;
    vector<int> node_to_lca_local_idx;

    void build_cell_lca_index();

    inline int cell_lca_rank(int primary_tn, int u, int v) const {
        if (primary_tn <= 0 || primary_tn >= static_cast<int>(cell_lca.size())) return 0;
        if (u <= 0 || v <= 0 || u >= static_cast<int>(node_to_lca_local_idx.size()) ||
            v >= static_cast<int>(node_to_lca_local_idx.size())) return 0;

        if ((u != primary_tn && node_to_primary_tn[u] != primary_tn) ||
            (v != primary_tn && node_to_primary_tn[v] != primary_tn)) return 0;

        const auto& idx = cell_lca[primary_tn];
        const int lu = node_to_lca_local_idx[u];
        const int lv = node_to_lca_local_idx[v];
        if (lu < 0 || lv < 0 || lu >= static_cast<int>(idx.nodes.size()) ||
            lv >= static_cast<int>(idx.nodes.size())) return 0;

        int left = idx.first[lu];
        int right = idx.first[lv];
        if (left > right) std::swap(left, right);

        const int len = right - left + 1;
        const int k = idx.log2[len];
        const int a = idx.sparse[k][left];
        const int b = idx.sparse[k][right - (1 << k) + 1];
        const int l = (idx.depth[a] <= idx.depth[b]) ? a : b;
        const int lca_node = idx.nodes[l];
        return (lca_node > 0 && lca_node < static_cast<int>(order_find.size()))
               ? order_find[lca_node] : 0;
    }


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

    // Build-time Local one-to-many Fast2Hop index statistics. Kept outside
    // QueryBenchmarkStats because benchmark resets must not erase build cost.
    unsigned long long local_f2h_index_entries_total = 0;
    long long local_f2h_index_cells_total = 0;
    double local_f2h_index_build_time_ms = 0.0;


    struct CellConn {
        int rt_cell_idx;
        int col;
    };


    struct DenseCell {
        int global_p_tn = 0;
        const int* tn_to_node_dists = nullptr;
        const int* min_access_dist = nullptr;
        int N = 0;
        vector<pair<int, int>> active_objects;
        vector<int> dense_access_tns;

        // Exact Global-phase acceleration. For each AccessTN column, object
        // rows are sorted by the immutable AccessTN->object distance. Query
        // processing can therefore stop the scan as soon as D_k is met.
        vector<uint32_t> sorted_object_offsets;
        vector<int> sorted_object_rows;
        bool sorted_objects_valid = false;

        // V8.1 Local one-to-many Fast2Hop index. Built once from the
        // immutable Fast2Hop labels and reused by queries whose source belongs
        // to this Cell. The query only intersects Source Label ranks with
        // these prebuilt Hub groups; no query-time bucket construction.
        Fast2Hop::CellOneToManyIndex f2h_one_to_many_index;
    };
    vector<DenseCell> rt_cells;

    // V7: query-time Global index is TN-centric.  Each dense TN owns one
    // distance-sorted list of Objects for which this TN is an exact AccessTN.
    // Cell/AccessTN matrices remain the build source, but Cells are no longer
    // traversed in the normal Global query path.
    struct TNObjectEntry {
        int object = 0;
        int distance = 1000000000;
    };
    vector<vector<TNObjectEntry>> tn_object_index;
    bool tn_object_index_valid = false;
    unsigned long long tn_object_index_raw_entries = 0;
    unsigned long long tn_object_index_entries = 0;
    unsigned long long tn_object_index_duplicates_removed = 0;

    vector<int> down_obj_count;
    TNR_2hop(int k ) : ET() {
        MAXk = k;
        // Propagate the experiment knob into ET so static Cell construction
        // and dynamic maintenance use the same alpha*K target.
        cell_capacity_factor = CELL_CAPACITY_FACTOR;
    }

    ~TNR_2hop() {
        if (global_top_k_heap) {
            delete[] global_top_k_heap;
            delete[] global_batch_targets;
            delete[] global_exact_dists;
            delete[] global_lca_ranks;
            delete[] global_dijkstra_heap;
        }
    }

    void init_cell_runtime_states();
    void build_sorted_cell_object_index(DenseCell& dcell);
    void build_local_f2h_cell_index(DenseCell& dcell);
    void build_TN_Object_Index();
    void build_TN_Network();
    void preprocess_KNN(string s , string o);
    void Add_Complete(int src);
    void Delete_Complete(int src);
    void quickQueryKNN(int src, int k, vector<pair<int, int>> &ans);

    void Apply_Index_Updates();
    inline bool Has_Pending_Maintenance() const { return maintenance_pending; }
    inline int Cell_Split_Threshold() const {
        return CellTargetObjectCount();
    }
    inline int Cell_Merge_Threshold() const {
        return std::max(1, static_cast<int>(std::floor(static_cast<double>(MAXk) / cell_capacity_factor)));
    }
    // Rebuild only the dynamic projection (TN/Cell/Matrix/Fast2Hop) on the
    // already-built static ET. Ordinary object updates remain O(1).
    void Rebuild_Dynamic_Index_From_Objects();
    void build_Fast2Hop();

    void Update_Dirty_Fast2Hop(const unordered_set<int>& dirty_p_tns, const vector<int>& dirty_nodes);
    void Recompute_Dirty_Cells(const unordered_set<int>& dirty_p_tns, const vector<int>& dirty_nodes);
    void Update_TN_Network(const unordered_set<int>& dirty_p_tns);


    void Update_Dirty_Core_Network(const unordered_set<int> &dirty_p_tns);
    vector<int> Reassemble_Dirty_Cells(const unordered_set<int>& dirty_p_tns);
    int Find_Best_Split_Point(int p_tn);

    std::queue<int> recycled_dense_ids;

    std::vector<std::vector<CellConn>> tn_connected_cells;

    void Update_Cell_Runtime_States(const unordered_set<int>& dirty_p_tns);

    // ===================== Internal benchmark / profiling =====================
    struct QueryBenchmarkStats {
        long long queries = 0;
        double total_query_time_us = 0.0;
        double total_local_time_us = 0.0;
        double total_global_time_us = 0.0;
        long long local_k_reached = 0;
        long long local_candidate_count = 0;
        long long local_object_count = 0;
        long long local_dk_valid_count = 0;
        long double local_dk_sum = 0.0;

        // Local 2-Hop micro breakdown.
        long long local_lca_queries = 0;
        double local_lca_time_us = 0.0;
        double local_2hop_time_us = 0.0;
        double local_object_eval_time_us = 0.0;

        // Fast2Hop V8 final-stage profiling: only metrics relevant to the
        // currently optimized one-to-many query path are retained.
        long long local_f2h_batches = 0;
        long long local_f2h_target_queries = 0;
        long long local_f2h_source_label_length_sum = 0;
        long long local_f2h_target_label_length_sum = 0;
        long long local_f2h_source_entries_scanned = 0;
        long long local_f2h_target_entries_scanned = 0;
        long long local_f2h_matched_hub_groups = 0;
        long long local_f2h_hub_object_scans = 0;
        long long local_f2h_hub_dk_breaks = 0;
        long long local_f2h_source_dk_breaks = 0;
        long long local_f2h_object_best_prunes = 0;
        long long local_f2h_distance_updates = 0;
        long long local_f2h_source_hub_sweeps = 0;
        long long local_f2h_pass1_time_ns = 0;
        long long local_f2h_pass2_time_ns = 0;
        long long local_f2h_source_cursor_advances = 0;
        long long local_f2h_hub_cursor_advances = 0;
        long long local_f2h_merge_comparisons = 0;
        long long local_f2h_physical_source_passes = 0;
        long long local_f2h_matched_hub_candidates = 0;
        long long local_f2h_lb_pruned_hub_groups = 0;
        long long local_f2h_processed_hub_groups = 0;
        long long local_f2h_target_lca_checks = 0;
        long long local_f2h_sort_comparisons = 0;
        long long local_f2h_sort_time_ns = 0;
        // V8.9 Object-Coalesced K-way profiling.  Heap metrics are now
        // object-level, not Hub-stream-level.
        long long local_f2h_kway_heap_pushes = 0;
        long long local_f2h_kway_heap_pops = 0;
        long long local_f2h_kway_object_group_creations = 0;
        long long local_f2h_kway_object_group_key_decreases = 0;
        long long local_f2h_kway_coalesced_streams = 0;
        long long local_f2h_kway_stream_entry_advances = 0;
        long long local_f2h_kway_finalized_object_skips = 0;
        long long local_f2h_kway_source_entry_skips = 0;
        long long local_f2h_kway_heap_peak_size = 0;
        long long local_f2h_kway_distinct_objects_popped = 0;

        long long access_tn = 0;
        long long tn_visited = 0;
        long long cell_visit_events = 0;
        long long cells_pruned_lb = 0;
        long long object_eval_events = 0;
        long long objects_evaluated = 0;
        long long distinct_cells = 0;
        long long distinct_objects = 0;
        long long sorted_object_scan_events = 0;
        long long sorted_object_dk_breaks = 0;

        // Work counters for explaining Global CPU time.
        long long global_tn_heap_pops = 0;
        long long global_tn_heap_pushes = 0;
        long long global_core_edge_checks = 0;
        long long global_core_edge_relax_success = 0;
        long long global_cell_conn_checks = 0;
        long long global_cell_lb_checks = 0;
        long long global_object_evaluate_calls = 0;
        long long global_tn_object_list_entries = 0;
        long long global_tn_object_scans = 0;
        long long global_tn_object_dk_breaks = 0;
        long long global_tn_object_index_hits = 0;

        // Optional phase timings; the final clean benchmark leaves collection disabled.
        double global_init_time_us = 0.0;
        double global_tn_pop_time_us = 0.0;
        double global_cell_phase_time_us = 0.0;
        double global_edge_phase_time_us = 0.0;
        double global_object_phase_time_us = 0.0;
        double global_tn_object_phase_time_us = 0.0;

        void reset() { *this = QueryBenchmarkStats(); }
    };

    QueryBenchmarkStats benchmark_stats;
    bool benchmark_collection_enabled = false;

    // Final benchmark export values. These are benchmark-only fields and do
    // not participate in the query algorithm. They allow the final main
    // program to write a clean, machine-readable CSV without parsing stdout.
    double last_clean_avg_query_us = 0.0;
    long long last_clean_query_count = 0;

    // Clean-only timing path: one warmup pass plus one timed pass over the
    // provided query file; never enables profiling or runs a profiling pass.
    // Returns average query latency (microseconds) over the measured input queries.
    double run_clean_query_benchmark(const vector<pair<int, int>>& queries,
                                     int warmup = -1);

    void print_statistics();
};

#endif
