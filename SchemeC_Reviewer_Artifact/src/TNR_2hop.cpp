#include <algorithm>
#include "TNR_2hop.h"
#include <vector>
#include <queue>
#include <iostream>
#include <utility>
#include <climits>

namespace fs = std::filesystem;



// Build one compact LCA structure for each primary Cell.
//
// The primary TN is kept as the root of the local subtree because it can be
// the true LCA of two local objects. Nested TNs are separators and therefore
// are not traversed: they belong to their own Cell.
void TNR_2hop::build_cell_lca_index() {
    cell_lca.assign(node_nums + 1, CellLCAIndex());
    node_to_lca_local_idx.assign(node_nums + 1, -1);

    struct Frame {
        int u;
        size_t next_child;
    };

    vector<Frame> stack;
    stack.reserve(64);

    for (int p_tn : core_tns) {
        if (!global_is_primary_tn[p_tn]) continue;

        CellLCAIndex& idx = cell_lca[p_tn];
        idx.nodes.clear();
        idx.euler.clear();
        idx.first.clear();
        idx.depth.clear();
        idx.log2.clear();
        idx.sparse.clear();

        // Local index 0 is always the primary TN / local tree root.
        idx.nodes.push_back(p_tn);
        idx.first.push_back(-1);
        idx.depth.push_back(0);
        node_to_lca_local_idx[p_tn] = 0;

        idx.euler.push_back(0);
        idx.first[0] = 0;

        stack.clear();
        stack.push_back({p_tn, 0});

        while (!stack.empty()) {
            Frame& f = stack.back();

            if (f.next_child >= static_tree_children[f.u].size()) {
                stack.pop_back();
                if (!stack.empty()) {
                    const int parent_local = node_to_lca_local_idx[stack.back().u];
                    idx.euler.push_back(parent_local);
                }
                continue;
            }

            const int v = static_tree_children[f.u][f.next_child++];

            // Only nodes owned by this primary Cell are included.  A nested
            // TN has node_to_primary_tn[v] == 0 and is deliberately cut here.
            if (v != p_tn && (global_is_tn[v] || node_to_primary_tn[v] != p_tn)) continue;

            const int local_id = static_cast<int>(idx.nodes.size());
            idx.nodes.push_back(v);
            idx.first.push_back(-1);
            idx.depth.push_back(idx.depth[node_to_lca_local_idx[f.u]] + 1);
            node_to_lca_local_idx[v] = local_id;

            idx.euler.push_back(local_id);
            if (idx.first[local_id] < 0) {
                idx.first[local_id] = static_cast<int>(idx.euler.size()) - 1;
            }

            stack.push_back({v, 0});
        }

        if (idx.euler.empty()) continue;

        idx.log2.assign(idx.euler.size() + 1, 0);
        for (size_t i = 2; i < idx.log2.size(); ++i) {
            idx.log2[i] = idx.log2[i >> 1] + 1;
        }

        const int levels = idx.log2.back() + 1;
        idx.sparse.assign(levels, vector<int>(idx.euler.size()));

        for (size_t i = 0; i < idx.euler.size(); ++i) {
            idx.sparse[0][i] = idx.euler[i];
        }

        for (int level = 1; level < levels; ++level) {
            const size_t span = static_cast<size_t>(1) << level;
            const size_t half = span >> 1;
            for (size_t i = 0; i + span <= idx.euler.size(); ++i) {
                const int a = idx.sparse[level - 1][i];
                const int b = idx.sparse[level - 1][i + half];
                idx.sparse[level][i] = (idx.depth[a] <= idx.depth[b]) ? a : b;
            }
        }
    }
}

// 1. 抽离出的独立 F2H 重构函数
void TNR_2hop::build_Fast2Hop() {
    cout << "--- [Phase 3.5] Building Pruned Fast2Hop Index ---" << endl;
    compute_f2h_levels();

    local_ancestors.assign(node_nums + 1, vector<pair<int, int>>());
    f2h.dynamic_labels.assign(node_nums + 1, vector<Fast2Hop::Label>());
    vector<int> temp_anc_dist(node_nums + 1, INT_MAX);
    vector<int> touched_anc;

    for (int i = 1; i <= node_nums; ++i) {
        if (!global_is_tn[i]) local_ancestors[i].push_back({0, i});
    }

    for (int rank = node_nums; rank >= 1; --rank) {
        int u = node_order[rank];

        if (global_is_tn[u] && node_f2h_level[u] > current_f2h_height) continue;

        for (size_t i = 0; i < bag_nodes[u].size(); ++i) {
            int v = bag_nodes[u][i];
            int weight_uv = bag_weights[u][i];
            if (weight_uv >= SAFE_INF) continue;

            if (!global_is_tn[v]) {
                for (auto& anc : local_ancestors[v]) {
                    int hub = anc.second;
                    int new_dist = safe_add(weight_uv, anc.first);
                    if (!global_is_tn[v] || node_f2h_level[v] <= current_f2h_height) {
                        if (temp_anc_dist[hub] == INT_MAX) touched_anc.push_back(hub);
                        if (new_dist < temp_anc_dist[hub]) {
                            temp_anc_dist[hub] = new_dist;
                        }
                    }
                }
            }
        }

        // 候选 Hub 按 Rank 降序排列（保障 Top-Down 剪枝正确性）
        std::sort(touched_anc.begin(), touched_anc.end(), [&](int a, int b) {
            return order_find[a] > order_find[b];
        });

        vector<Fast2Hop::Label> u_labels;
        for (auto& anc : local_ancestors[u]) {
            u_labels.push_back({static_cast<uint32_t>(order_find[anc.second]), anc.first});
        }

        for (int hub : touched_anc) {
            int cand_dist = temp_anc_dist[hub];
            temp_anc_dist[hub] = INT_MAX; // 恢复标记

            if (hub == u) continue;

            uint32_t hub_rank = static_cast<uint32_t>(order_find[hub]);

            // 2-Hop 路径剪枝校验 (Pruning Query)
            bool pruned = false;
            if (hub < (int)f2h.dynamic_labels.size()) {
                const auto& hub_labels = f2h.dynamic_labels[hub];
                size_t i = 0, j = 0;
                size_t u_len = u_labels.size();
                size_t hub_len = hub_labels.size();

                while (i < u_len && j < hub_len) {
                    uint32_t r1 = u_labels[i].rank;
                    uint32_t r2 = hub_labels[j].rank;

                    if (r1 == r2) {
                        int path_d = u_labels[i].dist + hub_labels[j].dist;
                        if (path_d <= cand_dist) {
                            pruned = true;
                            break;
                        }
                        i++; j++;
                    } else {
                        i += (r1 > r2);
                        j += (r1 < r2);
                    }
                }
            }

            if (!pruned) {
                u_labels.push_back({hub_rank, cand_dist});
                local_ancestors[u].push_back({cand_dist, hub});
            }
        }
        touched_anc.clear();

        // 及时更新 f2h.dynamic_labels[u]，供后续较低 Rank 节点剪枝使用
        f2h.update_local_node(u, local_ancestors[u], order_find);
    }

    // 将 Label 打平为绝对连续的平坦内存布局
    f2h.flatten();
}


void TNR_2hop::build_TN_Network() {
    const int num_tns = static_cast<int>(core_tns.size());
    global_to_dense.assign(node_nums + 1, -1);
    dense_to_global.assign(num_tns, 0);

    for (int i = 0; i < num_tns; ++i) {
        int global_u = core_tns[i];
        global_to_dense[global_u] = i;
        dense_to_global[i] = global_u;
    }

    // core_network is already the exact TN projection of the static CH/ET.
    // Do not apply an additional heuristic edge deletion here: that would
    // require a separate metric-preservation proof for every TN pair.
    size_t directed_core_edges = 0;
    size_t max_access_count = 0;
    for (int u : core_tns) {
        directed_core_edges += core_network[u].size();
        if (global_is_primary_tn[u]) {
            max_access_count = std::max(max_access_count, tn_cells[u].access_tns.size());
        }
        std::sort(core_network[u].begin(), core_network[u].end(),
                  [](const pair<int, int>& a, const pair<int, int>& b) {
                      if (a.second != b.second) return a.second < b.second;
                      return a.first < b.first;
                  });
    }

    max_pq_capacity = std::max<size_t>(1, directed_core_edges + max_access_count + 1);
    recycled_dense_ids = std::queue<int>();
    cout << " TN Exact Core Network ready. directed_edges="
         << directed_core_edges << endl;
}

void TNR_2hop::build_sorted_cell_object_index(DenseCell& dcell) {
    const int K = static_cast<int>(dcell.dense_access_tns.size());
    dcell.sorted_object_offsets.assign(static_cast<size_t>(K) + 1, 0);
    dcell.sorted_object_rows.clear();

    if (K <= 0 || dcell.active_objects.empty() || dcell.N <= 0 ||
        dcell.tn_to_node_dists == nullptr) {
        dcell.sorted_objects_valid = true;
        return;
    }

    dcell.sorted_object_rows.reserve(static_cast<size_t>(K) * dcell.active_objects.size());

    // Build one immutable sorted list per AccessTN.  The sort is intentionally
    // done during index construction/rebuild rather than in the query path.
    vector<pair<int, int>> tmp;
    tmp.reserve(dcell.active_objects.size());

    for (int col = 0; col < K; ++col) {
        tmp.clear();
        const int* base = dcell.tn_to_node_dists + static_cast<size_t>(col) * dcell.N;

        for (const auto& obj : dcell.active_objects) {
            const int row = obj.first;
            const int dist = base[row];
            if (dist < SAFE_INF) {
                tmp.push_back({dist, row});
            }
        }

        sort(tmp.begin(), tmp.end(), [](const pair<int, int>& a, const pair<int, int>& b) {
            if (a.first != b.first) return a.first < b.first;
            return a.second < b.second;
        });

        dcell.sorted_object_offsets[col] = static_cast<uint32_t>(dcell.sorted_object_rows.size());
        for (const auto& item : tmp) {
            dcell.sorted_object_rows.push_back(item.second);
        }
        dcell.sorted_object_offsets[col + 1] = static_cast<uint32_t>(dcell.sorted_object_rows.size());
    }

    dcell.sorted_objects_valid = true;
}

void TNR_2hop::build_TN_Object_Index() {
    const size_t tn_count = dense_to_global.size();
    tn_object_index.assign(tn_count, {});
    tn_object_index_valid = false;
    tn_object_index_raw_entries = 0;
    tn_object_index_entries = 0;
    tn_object_index_duplicates_removed = 0;

    if (tn_count == 0) {
        tn_object_index_valid = true;
        return;
    }

    // Build directly from the exact Cell matrices.  A TN contributes an
    // (Object, TN->Object distance) pair for every Cell in which it is an
    // AccessTN.  The same Object may appear through several Cells; the
    // duplicate is removed by taking the minimum exact distance below.
    for (int p_tn : core_tns) {
        if (!global_is_primary_tn[p_tn]) continue;

        const DenseCell& dcell = rt_cells[p_tn];
        if (dcell.dense_access_tns.empty() || dcell.active_objects.empty() ||
            dcell.tn_to_node_dists == nullptr || dcell.N <= 0) {
            continue;
        }

        for (int col = 0; col < static_cast<int>(dcell.dense_access_tns.size()); ++col) {
            const int dense_tn = dcell.dense_access_tns[col];
            if (dense_tn < 0 || dense_tn >= static_cast<int>(tn_object_index.size())) continue;

            const int* base = dcell.tn_to_node_dists + static_cast<size_t>(col) * dcell.N;
            auto& out = tn_object_index[dense_tn];
            out.reserve(out.size() + dcell.active_objects.size());

            for (const auto& act_obj : dcell.active_objects) {
                const int row = act_obj.first;
                const int object = act_obj.second;
                if (row < 0 || row >= dcell.N) continue;
                const int d = base[row];
                if (d >= SAFE_INF) continue;
                out.push_back({object, d});
                ++tn_object_index_raw_entries;
            }
        }
    }

    for (auto& entries : tn_object_index) {
        if (entries.empty()) continue;

        // First canonicalize by Object ID so the same Object reached through
        // multiple Cell/AccessTN relationships is represented exactly once.
        std::sort(entries.begin(), entries.end(), [](const TNObjectEntry& a, const TNObjectEntry& b) {
            if (a.object != b.object) return a.object < b.object;
            return a.distance < b.distance;
        });

        size_t write = 0;
        for (size_t read = 0; read < entries.size();) {
            const int object = entries[read].object;
            int best_distance = entries[read].distance;
            size_t next = read + 1;
            while (next < entries.size() && entries[next].object == object) {
                if (entries[next].distance < best_distance) best_distance = entries[next].distance;
                ++next;
            }
            entries[write++] = {object, best_distance};
            tn_object_index_duplicates_removed += static_cast<unsigned long long>(next - read - 1);
            read = next;
        }
        entries.resize(write);

        // Query order: exact TN->Object distance ascending.  Once
        // curr_dist + entry.distance >= D_k, every later Object is also unable
        // to improve the current Top-K and the scan can terminate exactly.
        std::sort(entries.begin(), entries.end(), [](const TNObjectEntry& a, const TNObjectEntry& b) {
            if (a.distance != b.distance) return a.distance < b.distance;
            return a.object < b.object;
        });

        tn_object_index_entries += entries.size();
    }

    tn_object_index_valid = true;

    cout << " [Global V7] TN->Object Index built. raw="
         << tn_object_index_raw_entries
         << ", unique=" << tn_object_index_entries
         << ", duplicates_removed=" << tn_object_index_duplicates_removed << endl;
}

void TNR_2hop::build_local_f2h_cell_index(DenseCell& dcell) {
    // Preserve DenseCell::active_objects order as the query-time target slot.
    // This removes the old object-id -> slot hash/stamp lookup completely.
    f2h.build_cell_one_to_many_index(dcell.active_objects, dcell.f2h_one_to_many_index);
}

void TNR_2hop::init_cell_runtime_states() {
    cout << "\n=== [Phase 4]  ===" << endl;
    unsigned long long local_f2h_index_entries = 0;
    long long local_f2h_index_cells = 0;

    node_to_cell_row.assign(node_nums + 1, -1);
    node_to_active_idx.assign(node_nums + 1, -1);


    rt_cells.assign(node_nums + 1, DenseCell());

    // Each Cell owns its compact LCA/RMQ state.  This is rebuilt together
    // with the dynamic Cell projection, never over the entire ET.
    build_cell_lca_index();

    tn_connected_cells.assign(dense_to_global.size() + 100, vector<CellConn>());

    const auto local_f2h_index_build_t0 = steady_clock::now();
    for (int p_tn : core_tns) {
        if (!global_is_primary_tn[p_tn]) continue;


        DenseCell& dcell = rt_cells[p_tn];
        Cell& cell = tn_cells[p_tn];

        dcell.global_p_tn = p_tn;
        dcell.N = static_cast<int>(cell.local_nodes.size());
        dcell.tn_to_node_dists = cell.tn_to_node_dists.data();
        dcell.min_access_dist = cell.min_to_access_dist.data();
        dcell.active_objects.reserve(dcell.N);

        int K = cell.access_tns.size();
        dcell.dense_access_tns.resize(K);

        for (int row = 0; row < dcell.N; ++row) {
            int v = cell.local_nodes[row];
            node_to_cell_row[v] = row;

            if (nodes[v].is_Object) {
                node_to_active_idx[v] = dcell.active_objects.size();
                dcell.active_objects.push_back({row, v});
            }
        }

        for (int col = 0; col < K; ++col) {
            int acc_tn = cell.access_tns[col];
            int dense_acc = global_to_dense[acc_tn];
            dcell.dense_access_tns[col] = dense_acc;

            tn_connected_cells[dense_acc].push_back({p_tn, col});
        }

        build_sorted_cell_object_index(dcell);
        build_local_f2h_cell_index(dcell);
        if (dcell.f2h_one_to_many_index.valid) {
            local_f2h_index_entries += dcell.f2h_one_to_many_index.entry_count;
            ++local_f2h_index_cells;
        }
    }

    const double local_f2h_index_build_ms =
        duration<double, std::milli>(steady_clock::now() - local_f2h_index_build_t0).count();
    local_f2h_index_entries_total = local_f2h_index_entries;
    local_f2h_index_cells_total = local_f2h_index_cells;
    local_f2h_index_build_time_ms = local_f2h_index_build_ms;
    cout << " [Local V8.1] Cell Hub Index built. cells=" << local_f2h_index_cells
         << ", entries=" << local_f2h_index_entries
         << ", build=" << local_f2h_index_build_ms << " ms" << endl;

    // V7 query path is TN-centric; build it only after every runtime Cell has
    // valid active-object lists and exact AccessTN matrices.
    build_TN_Object_Index();
}

void TNR_2hop::preprocess_KNN(string s, string o) {
    cout << "===  start  ===" << endl;
    auto start_time = steady_clock::now();

    load_raw_graph(s + "Ein.txt");
    addObject(o);
    build_et(MAXk);

    dist.assign(node_nums + 1, INT_MAX);
    version.assign(node_nums + 1, 0);
    tn_states.assign(core_tns.size(), TNDenseState());
    top_k_pos.assign(node_nums + 1, -1);

    build_Fast2Hop();
    build_TN_Network();
    init_cell_runtime_states();

    if (global_top_k_heap == nullptr) {
        global_top_k_heap = new TopKNode[LIMIT_K];
        global_batch_targets = new int[node_nums + 1];
        global_exact_dists = new int[node_nums + 1];
        global_lca_ranks = new int[node_nums + 1];
        global_dijkstra_heap = new HeapNode[max_pq_capacity];
    }

    auto end_time = steady_clock::now();
    cout  << duration_cast<milliseconds>(end_time - start_time).count() << " ms ===" << endl;
}

void TNR_2hop::quickQueryKNN(int src, int k ,vector<pair<int, int>> &ans) {
    // Guard the public query entry point: the fixed Top-K scratch heap supports
    // k <= LIMIT_K, and all per-node arrays are indexed by valid 1-based IDs.
    if (src <= 0 || src >= static_cast<int>(nodes.size()) ||
        src >= static_cast<int>(node_to_primary_tn.size()) ||
        k <= 0 || k > LIMIT_K || global_top_k_heap == nullptr) {
        ans.clear();
        return;
    }
    const bool fine_profiling = benchmark_collection_enabled && ENABLE_FINE_GLOBAL_PROFILING;
    const auto query_bench_start = benchmark_collection_enabled
        ? steady_clock::now() : steady_clock::time_point{};
    vector<int> query_visited_cells;
    vector<int> query_evaluated_objects;
    if (benchmark_collection_enabled) {
        query_visited_cells.reserve(64);
        query_evaluated_objects.reserve(128);
    }
    if (++cur_src_ver == 0) cur_src_ver = 1;
    if (++current_version == 0) {
        std::fill(version.begin(), version.end(), 0);
        current_version = 1;
    }
    D_k = INT_MAX;
    current_k_count = 0;
    TopKNode* top_k_heap = global_top_k_heap;

    auto evaluate_obj = [&](int obj_id, int path_d) {
        if (benchmark_collection_enabled) {
            ++benchmark_stats.object_eval_events;
            query_evaluated_objects.push_back(obj_id);
        }
        if (current_k_count == k && path_d >= D_k) return;
        if (benchmark_collection_enabled) ++benchmark_stats.objects_evaluated;
        auto heapify_down = [&](int pos) {
            int target_dist = top_k_heap[pos].dist;
            int target_id = top_k_heap[pos].id;
            int limit = current_k_count;
            while (true) {
                int left = 2 * pos + 1;
                int right = 2 * pos + 2;
                int largest = pos;
                int current_max = target_dist;

                if (left < limit && top_k_heap[left].dist > current_max) {
                    largest = left;
                    current_max = top_k_heap[left].dist;
                }
                if (right < limit && top_k_heap[right].dist > current_max) {
                    largest = right;
                }

                if (largest != pos) {
                    top_k_heap[pos] = top_k_heap[largest];
                    top_k_pos[top_k_heap[pos].id] = pos;
                    pos = largest;
                } else {
                    break;
                }
            }
            top_k_heap[pos] = {target_dist, target_id};
            top_k_pos[target_id] = pos;
        };

        auto heapify_up = [&](int pos) {
            int target_dist = top_k_heap[pos].dist;
            int target_id = top_k_heap[pos].id;
            while (pos > 0) {
                int parent = (pos - 1) / 2;
                if (target_dist > top_k_heap[parent].dist) {
                    top_k_heap[pos] = top_k_heap[parent];
                    top_k_pos[top_k_heap[pos].id] = pos;
                    pos = parent;
                } else {
                    break;
                }
            }
            top_k_heap[pos] = {target_dist, target_id};
            top_k_pos[target_id] = pos;
        };

        if (version[obj_id] == current_version) {
            if (path_d >= dist[obj_id]) return;
            dist[obj_id] = path_d;

            int pos = top_k_pos[obj_id];
            if (pos != -1) {
                top_k_heap[pos].dist = path_d;
                heapify_down(pos);
                if (current_k_count == k) D_k = top_k_heap[0].dist;
                return;
            }
        } else {
            version[obj_id] = current_version;
            dist[obj_id] = path_d;
        }

        if (current_k_count < k) {
            int idx = current_k_count++;
            top_k_heap[idx] = {path_d, obj_id};
            top_k_pos[obj_id] = idx;
            heapify_up(idx);
        } else {
            int old_root_id = top_k_heap[0].id;
            top_k_pos[old_root_id] = -1;
            top_k_heap[0] = {path_d, obj_id};
            top_k_pos[obj_id] = 0;
            heapify_down(0);
        }

        if (current_k_count == k) D_k = top_k_heap[0].dist;
    };

    auto start = benchmark_collection_enabled
        ? steady_clock::now() : steady_clock::time_point{};
    microseconds last{0};

    if (nodes[src].is_Object) evaluate_obj(src, 0);

    int p_tn_src = node_to_primary_tn[src];

    if (p_tn_src != 0) {
        int batch_count = 0;
        int* batch_targets = global_batch_targets;
        int* batch_lca_ranks = global_lca_ranks;

        DenseCell& dcell = rt_cells[p_tn_src];

        if (benchmark_collection_enabled) {
            benchmark_stats.local_object_count += static_cast<long long>(dcell.active_objects.size());
        }

        const auto local_lca_t0 = fine_profiling
            ? steady_clock::now() : steady_clock::time_point{};
        const bool use_index = dcell.f2h_one_to_many_index.valid;
        for (const auto& act_obj : dcell.active_objects) {
            const int obj_id = act_obj.second;
            batch_targets[batch_count] = obj_id;
            // V8.9: when the static one-to-many index is valid, every HubObjectEntry
            // already certifies target-label membership, and Source-side compact
            // Hub-rank lookup certifies source-label membership. No per-target LCA is needed.
            // Compute LCA ranks only for the dynamic fallback path.
            if (use_index) {
                batch_lca_ranks[batch_count] = UINT_MAX;
            } else if (obj_id == src) {
                batch_lca_ranks[batch_count] = UINT_MAX;
            } else {
                batch_lca_ranks[batch_count] = cell_lca_rank(p_tn_src, src, obj_id);
                if (benchmark_collection_enabled) ++benchmark_stats.local_lca_queries;
            }
            ++batch_count;
            nums_local++;
        }
        if (benchmark_collection_enabled && ENABLE_FINE_GLOBAL_PROFILING) {
            benchmark_stats.local_lca_time_us +=
                duration<double, std::micro>(steady_clock::now() - local_lca_t0).count();
        }

        if (batch_count > 0) {
            int* exact_dists = global_exact_dists;

            const auto local_2hop_t0 = fine_profiling
                ? steady_clock::now() : steady_clock::time_point{};
            Fast2Hop::BatchQueryStats f2h_query_stats;
            auto f2h_evaluator = [&](int obj_id, int path_d) {
                if (benchmark_collection_enabled) ++benchmark_stats.local_candidate_count;
                evaluate_obj(obj_id, path_d);
            };

            if (dcell.f2h_one_to_many_index.valid) {
                f2h.query_1_to_N_hub_grouped(
                    src, batch_targets, batch_lca_ranks, batch_count, exact_dists, D_k,
                    dcell.f2h_one_to_many_index, f2h_evaluator,
                    benchmark_collection_enabled ? &f2h_query_stats : nullptr);
            } else {
                // Dynamic object updates invalidate the Cell Hub index lazily.
                // Keep an exact Pairwise F2H fallback until the next rebuild.
                f2h.query_1_to_N_pairwise(
                    src, batch_targets, batch_lca_ranks, batch_count, exact_dists, D_k,
                    nullptr);
                for (int i = 0; i < batch_count; ++i) {
                    if (exact_dists[i] != INT_MAX) f2h_evaluator(batch_targets[i], exact_dists[i]);
                }
            }
            if (benchmark_collection_enabled) {
                benchmark_stats.local_f2h_batches += f2h_query_stats.batch_queries;
                benchmark_stats.local_f2h_target_queries += f2h_query_stats.target_queries;
                benchmark_stats.local_f2h_source_label_length_sum += f2h_query_stats.source_label_length_sum;
                benchmark_stats.local_f2h_target_label_length_sum += f2h_query_stats.target_label_length_sum;
                benchmark_stats.local_f2h_source_entries_scanned += f2h_query_stats.source_entries_scanned;
                benchmark_stats.local_f2h_target_entries_scanned += f2h_query_stats.target_entries_scanned;
                benchmark_stats.local_f2h_matched_hub_groups += f2h_query_stats.matched_hub_groups;
                benchmark_stats.local_f2h_distance_updates += f2h_query_stats.distance_updates;
                benchmark_stats.local_f2h_pass1_time_ns += f2h_query_stats.pass1_time_ns;
                benchmark_stats.local_f2h_pass2_time_ns += f2h_query_stats.pass2_time_ns;
                benchmark_stats.local_f2h_hub_object_scans += f2h_query_stats.hub_object_scans;
                benchmark_stats.local_f2h_hub_dk_breaks += f2h_query_stats.hub_dk_breaks;
                benchmark_stats.local_f2h_source_dk_breaks += f2h_query_stats.source_dk_breaks;
                benchmark_stats.local_f2h_object_best_prunes += f2h_query_stats.object_best_prunes;
                benchmark_stats.local_f2h_source_cursor_advances += f2h_query_stats.source_cursor_advances;
                benchmark_stats.local_f2h_hub_cursor_advances += f2h_query_stats.hub_cursor_advances;
                benchmark_stats.local_f2h_merge_comparisons += f2h_query_stats.merge_comparisons;
                benchmark_stats.local_f2h_physical_source_passes += f2h_query_stats.physical_source_passes;
                benchmark_stats.local_f2h_matched_hub_candidates += f2h_query_stats.matched_hub_candidates;
                benchmark_stats.local_f2h_lb_pruned_hub_groups += f2h_query_stats.lb_pruned_hub_groups;
                benchmark_stats.local_f2h_processed_hub_groups += f2h_query_stats.processed_hub_groups;
                benchmark_stats.local_f2h_target_lca_checks += f2h_query_stats.target_lca_checks;
                benchmark_stats.local_f2h_sort_comparisons += f2h_query_stats.sort_comparisons;
                benchmark_stats.local_f2h_sort_time_ns += f2h_query_stats.sort_time_ns;
                benchmark_stats.local_f2h_kway_heap_pushes += f2h_query_stats.kway_heap_pushes;
                benchmark_stats.local_f2h_kway_heap_pops += f2h_query_stats.kway_heap_pops;
                benchmark_stats.local_f2h_kway_object_group_creations += f2h_query_stats.kway_object_group_creations;
                benchmark_stats.local_f2h_kway_object_group_key_decreases += f2h_query_stats.kway_object_group_key_decreases;
                benchmark_stats.local_f2h_kway_coalesced_streams += f2h_query_stats.kway_coalesced_streams;
                benchmark_stats.local_f2h_kway_stream_entry_advances += f2h_query_stats.kway_stream_entry_advances;
                benchmark_stats.local_f2h_kway_finalized_object_skips += f2h_query_stats.kway_finalized_object_skips;
                benchmark_stats.local_f2h_kway_source_entry_skips += f2h_query_stats.kway_source_entry_skips;
                benchmark_stats.local_f2h_kway_heap_peak_size += f2h_query_stats.kway_heap_peak_size;
                benchmark_stats.local_f2h_kway_distinct_objects_popped += f2h_query_stats.kway_distinct_objects_popped;
            }
            if (benchmark_collection_enabled && ENABLE_FINE_GLOBAL_PROFILING) {
                benchmark_stats.local_2hop_time_us +=
                    duration<double, std::micro>(steady_clock::now() - local_2hop_t0).count();
            }
        }
    }

    if (benchmark_collection_enabled) {
        last = duration_cast<microseconds>(steady_clock::now() - start);
        local += last.count();
        benchmark_stats.total_local_time_us += static_cast<double>(last.count());
        if (current_k_count >= k) ++benchmark_stats.local_k_reached;
        if (D_k < INT_MAX) {
            ++benchmark_stats.local_dk_valid_count;
            benchmark_stats.local_dk_sum += static_cast<long double>(D_k);
        }
    }
    start = benchmark_collection_enabled
        ? steady_clock::now() : steady_clock::time_point{};

    HeapNode* pq_arr = global_dijkstra_heap;
    int pq_size = 0;
    auto cmp = [](const HeapNode& a, const HeapNode& b) { return a.dist > b.dist; };

    const auto global_init_t0 = fine_profiling
        ? steady_clock::now() : steady_clock::time_point{};

    if (global_is_tn[src]) {
        pq_arr[pq_size++] = {0, global_to_dense[src]};
        if (benchmark_collection_enabled) ++benchmark_stats.global_tn_heap_pushes;
    } else if (p_tn_src != 0) {
        DenseCell& dcell = rt_cells[p_tn_src];
        int row = node_to_cell_row[src];

        for (int col = 0; col < (int)dcell.dense_access_tns.size(); ++col) {
            int dense_tn = dcell.dense_access_tns[col];
            int dist_s_to_tn = dcell.tn_to_node_dists[col * dcell.N + row];
            if (benchmark_collection_enabled) {
                ++ace;
                ++benchmark_stats.access_tn;
            }
            if (dist_s_to_tn < SAFE_INF && dist_s_to_tn < D_k) {
                tn_states[dense_tn].push_ver = cur_src_ver;
                tn_states[dense_tn].best_tn_dist = dist_s_to_tn;
                pq_arr[pq_size++] = {dist_s_to_tn, dense_tn};
                if (benchmark_collection_enabled) ++benchmark_stats.global_tn_heap_pushes;
            }
        }
    }
    std::make_heap(pq_arr, pq_arr + pq_size, cmp);

    if (benchmark_collection_enabled && ENABLE_FINE_GLOBAL_PROFILING) {
        benchmark_stats.global_init_time_us +=
            duration<double, std::micro>(steady_clock::now() - global_init_t0).count();
    }

    if (current_k_count == k && pq_size > 0 && pq_arr[0].dist >= D_k) pq_size = 0;

    while (pq_size > 0  && pq_arr[0].dist < D_k) {
        const auto pop_t0 = fine_profiling
            ? steady_clock::now() : steady_clock::time_point{};
        std::pop_heap(pq_arr, pq_arr + pq_size, cmp);
        auto [curr_dist, dense_u] = pq_arr[--pq_size];
        if (benchmark_collection_enabled && ENABLE_FINE_GLOBAL_PROFILING) {
            benchmark_stats.global_tn_pop_time_us +=
                duration<double, std::micro>(steady_clock::now() - pop_t0).count();
        }
        if (benchmark_collection_enabled) ++benchmark_stats.global_tn_heap_pops;

        if (curr_dist >= D_k) break;

        if (tn_states[dense_u].settled_ver == cur_src_ver) continue;
        tn_states[dense_u].settled_ver = cur_src_ver;
        if (benchmark_collection_enabled) ++benchmark_stats.tn_visited;

        int global_u = dense_to_global[dense_u];
        if (nodes[global_u].is_Object) {
            if (benchmark_collection_enabled) ++benchmark_stats.global_object_evaluate_calls;
            evaluate_obj(global_u, curr_dist);
        }

        // ============================================================
        // V7 Global query: TN -> Object.
        //
        // Dijkstra settles TNs, so Objects are expanded under the settled TN
        // itself rather than by traversing TN -> Cell -> Object.  The index is
        // built from exact Cell matrices, but the Cell layer disappears from
        // the normal query path.  This removes repeated CellConn traversal.
        //
        // If an object update has invalidated the immutable TN->Object index,
        // fall back to the V6 exact Cell path.  This preserves cheap lazy
        // updates without making the query incorrect.
        // ============================================================
        if (tn_object_index_valid) {
            const auto object_phase_t0 = fine_profiling
                ? steady_clock::now() : steady_clock::time_point{};
            if (dense_u >= 0 && dense_u < static_cast<int>(tn_object_index.size())) {
                const auto& entries = tn_object_index[dense_u];
                if (benchmark_collection_enabled) {
                    ++benchmark_stats.global_tn_object_index_hits;
                    benchmark_stats.global_tn_object_list_entries += static_cast<long long>(entries.size());
                }

                bool broke_by_dk = false;
                for (const TNObjectEntry& entry : entries) {
                    const int total_path = safe_add(curr_dist, entry.distance);
                    if (benchmark_collection_enabled) {
                        ++benchmark_stats.global_tn_object_scans;
                    }

                    if (current_k_count == k && total_path >= D_k) {
                        broke_by_dk = true;
                        if (benchmark_collection_enabled) ++benchmark_stats.global_tn_object_dk_breaks;
                        break;
                    }

                    if (total_path < D_k) {
                        if (benchmark_collection_enabled) ++benchmark_stats.global_object_evaluate_calls;
                        evaluate_obj(entry.object, total_path);
                    }
                }
                (void)broke_by_dk;
            }

            if (benchmark_collection_enabled && ENABLE_FINE_GLOBAL_PROFILING) {
                benchmark_stats.global_tn_object_phase_time_us +=
                    duration<double, std::micro>(steady_clock::now() - object_phase_t0).count();
            }
        } else {
            // Exact V6 fallback for lazy object updates.  This branch is not
            // used in a clean static/rebuilt benchmark.
            const auto cell_phase_t0 = fine_profiling
                ? steady_clock::now() : steady_clock::time_point{};
            for (const CellConn& conn : tn_connected_cells[dense_u]) {
                if (benchmark_collection_enabled) ++benchmark_stats.global_cell_conn_checks;
                DenseCell& dcell = rt_cells[conn.rt_cell_idx];
                if (dcell.N <= 0 || dcell.tn_to_node_dists == nullptr) continue;

                if (dcell.min_access_dist == nullptr || conn.col < 0 ||
                    conn.col >= static_cast<int>(dcell.dense_access_tns.size())) continue;

                const int min_local_dist = dcell.min_access_dist[conn.col];
                if (benchmark_collection_enabled) ++benchmark_stats.global_cell_lb_checks;
                if (safe_add(curr_dist, min_local_dist) >= D_k) {
                    if (benchmark_collection_enabled) ++benchmark_stats.cells_pruned_lb;
                    continue;
                }

                if (benchmark_collection_enabled) {
                    ++benchmark_stats.cell_visit_events;
                    query_visited_cells.push_back(conn.rt_cell_idx);
                }

                const int* base_matrix = dcell.tn_to_node_dists + (conn.col * dcell.N);
                const auto object_phase_t0 = fine_profiling
                    ? steady_clock::now() : steady_clock::time_point{};
                if (dcell.sorted_objects_valid && !dcell.sorted_object_offsets.empty()) {
                    const uint32_t begin = dcell.sorted_object_offsets[conn.col];
                    const uint32_t end = dcell.sorted_object_offsets[conn.col + 1];
                    for (uint32_t pos = begin; pos < end; ++pos) {
                        const int row = dcell.sorted_object_rows[pos];
                        const int exact_dist = base_matrix[row];
                        const int total_path = safe_add(curr_dist, exact_dist);
                        if (benchmark_collection_enabled) ++benchmark_stats.sorted_object_scan_events;
                        if (current_k_count == k && total_path >= D_k) {
                            if (benchmark_collection_enabled) ++benchmark_stats.sorted_object_dk_breaks;
                            break;
                        }
                        if (total_path < D_k) {
                            if (benchmark_collection_enabled) ++benchmark_stats.global_object_evaluate_calls;
                            evaluate_obj(tn_cells[conn.rt_cell_idx].local_nodes[row], total_path);
                        }
                    }
                } else {
                    for (const auto& act_obj : dcell.active_objects) {
                        const int exact_dist = base_matrix[act_obj.first];
                        if (exact_dist >= SAFE_INF) continue;
                        const int total_path = safe_add(curr_dist, exact_dist);
                        if (total_path < D_k) {
                            if (benchmark_collection_enabled) ++benchmark_stats.global_object_evaluate_calls;
                            evaluate_obj(act_obj.second, total_path);
                        }
                    }
                }
                if (benchmark_collection_enabled && ENABLE_FINE_GLOBAL_PROFILING) {
                    benchmark_stats.global_object_phase_time_us +=
                        duration<double, std::micro>(steady_clock::now() - object_phase_t0).count();
                }
            }
            if (benchmark_collection_enabled && ENABLE_FINE_GLOBAL_PROFILING) {
                benchmark_stats.global_cell_phase_time_us +=
                    duration<double, std::micro>(steady_clock::now() - cell_phase_t0).count();
            }
        }

        const auto edge_phase_t0 = fine_profiling
            ? steady_clock::now() : steady_clock::time_point{};
        for (const auto& edge : core_network[dense_to_global[dense_u]]) {
            if (benchmark_collection_enabled) ++benchmark_stats.global_core_edge_checks;
            int w_uv = edge.second;
            int dense_v = global_to_dense[edge.first];

            int edge_bound = safe_add(curr_dist, w_uv);
            if (edge_bound >= D_k) break;

            if (tn_states[dense_v].settled_ver != cur_src_ver) {
                int next_d = edge_bound;
                if (tn_states[dense_v].push_ver != cur_src_ver) {
                    tn_states[dense_v].best_tn_dist = INT_MAX;
                }
                if (next_d < tn_states[dense_v].best_tn_dist) {
                    tn_states[dense_v].push_ver = cur_src_ver;
                    tn_states[dense_v].best_tn_dist = next_d;
                    pq_arr[pq_size++] = {next_d, dense_v};
                    std::push_heap(pq_arr, pq_arr + pq_size, cmp);
                    if (benchmark_collection_enabled) {
                        ++benchmark_stats.global_tn_heap_pushes;
                        ++benchmark_stats.global_core_edge_relax_success;
                    }
                    if (benchmark_collection_enabled) ++up_tn;
                }
            }
        }
        if (benchmark_collection_enabled && ENABLE_FINE_GLOBAL_PROFILING) {
            benchmark_stats.global_edge_phase_time_us +=
                duration<double, std::micro>(steady_clock::now() - edge_phase_t0).count();
        }
    }
    if (benchmark_collection_enabled) {
        last = duration_cast<microseconds>(steady_clock::now() - start);
        gloab += last.count();
        benchmark_stats.total_global_time_us += static_cast<double>(last.count());

        std::sort(query_visited_cells.begin(), query_visited_cells.end());
        query_visited_cells.erase(std::unique(query_visited_cells.begin(), query_visited_cells.end()), query_visited_cells.end());
        benchmark_stats.distinct_cells += static_cast<long long>(query_visited_cells.size());

        std::sort(query_evaluated_objects.begin(), query_evaluated_objects.end());
        query_evaluated_objects.erase(std::unique(query_evaluated_objects.begin(), query_evaluated_objects.end()), query_evaluated_objects.end());
        benchmark_stats.distinct_objects += static_cast<long long>(query_evaluated_objects.size());

        ++benchmark_stats.queries;
        benchmark_stats.total_query_time_us += duration<double, std::micro>(steady_clock::now() - query_bench_start).count();
    }
    // The API returns exactly min(k, |O|) answers.  Retaining entries from a
    // previous query makes an under-filled KNN result appear to contain stale
    // objects.  resize preserves capacity for the usual preallocated caller.
    ans.resize(current_k_count);
    for (int i = 0; i < current_k_count; ++i) ans[i] = {top_k_heap[i].dist, top_k_heap[i].id};

    std::sort(ans.begin(), ans.end(), [](const pair<int, int>& a, const pair<int, int>& b) {
        if (a.first != b.first) return a.first < b.first;
        return a.second < b.second;
    });
}

void TNR_2hop::Add_Complete(int src) {
    // Public update API: reject invalid IDs before indexing any per-node vector.
    if (src <= 0 || src >= static_cast<int>(nodes.size())) return;
    if (nodes[src].is_Object) return;

    int p_tn = 0;
    int c_row = -1;
    if (!global_is_tn[src]) {
        // Every non-TN vertex must be assigned to a primary Cell. Validate this
        // invariant before changing nodes[src].is_Object, otherwise an invalid
        // mapping could silently create an Object that no query path can see.
        if (src >= static_cast<int>(node_to_primary_tn.size()) ||
            src >= static_cast<int>(node_to_cell_row.size()) ||
            src >= static_cast<int>(node_to_active_idx.size())) {
            cerr << "[Update Error] Add_Complete: missing node-to-Cell mapping for vertex "
                 << src << '\n';
            return;
        }
        p_tn = node_to_primary_tn[src];
        c_row = node_to_cell_row[src];
        if (p_tn <= 0 || p_tn >= static_cast<int>(rt_cells.size()) ||
            c_row < 0 || c_row >= rt_cells[p_tn].N) {
            cerr << "[Update Error] Add_Complete: invalid Cell mapping for vertex "
                 << src << " (primary TN=" << p_tn << ", row=" << c_row << ")\n";
            return;
        }
    }

    nodes[src].is_Object = true;
    ++nums_object;

    // TN Objects are evaluated directly when the corresponding TN is settled;
    // they are not stored in a Cell's active-object list.
    if (global_is_tn[src]) return;

    DenseCell& dcell = rt_cells[p_tn];
    node_to_active_idx[src] = static_cast<int>(dcell.active_objects.size());
    dcell.active_objects.push_back({c_row, src});
    dcell.sorted_objects_valid = false;
    dcell.f2h_one_to_many_index.valid = false;
    tn_object_index_valid = false;
    maintenance_pending = true;
    lazyupdated = true;
}

void TNR_2hop::Delete_Complete(int src) {
    // Public update API: reject invalid IDs before indexing any per-node vector.
    if (src <= 0 || src >= static_cast<int>(nodes.size())) return;
    if (!nodes[src].is_Object) return;

    int p_tn = 0;
    int del_idx = -1;
    if (!global_is_tn[src]) {
        if (src >= static_cast<int>(node_to_primary_tn.size()) ||
            src >= static_cast<int>(node_to_active_idx.size())) {
            cerr << "[Update Error] Delete_Complete: missing node-to-Cell mapping for vertex "
                 << src << '\n';
            return;
        }
        p_tn = node_to_primary_tn[src];
        del_idx = node_to_active_idx[src];
        if (p_tn <= 0 || p_tn >= static_cast<int>(rt_cells.size()) ||
            del_idx < 0 || del_idx >= static_cast<int>(rt_cells[p_tn].active_objects.size()) ||
            rt_cells[p_tn].active_objects[del_idx].second != src) {
            cerr << "[Update Error] Delete_Complete: inconsistent active-object mapping for vertex "
                 << src << " (primary TN=" << p_tn << ", active index=" << del_idx << ")\n";
            return;
        }
    }

    nodes[src].is_Object = false;
    --nums_object;

    // TN Objects are not members of a Cell active-object list.
    if (global_is_tn[src]) return;

    DenseCell& dcell = rt_cells[p_tn];
    const int last_idx = static_cast<int>(dcell.active_objects.size()) - 1;
    if (del_idx != last_idx) {
        const auto last_obj = dcell.active_objects[last_idx];
        dcell.active_objects[del_idx] = last_obj;
        node_to_active_idx[last_obj.second] = del_idx;
    }
    dcell.active_objects.pop_back();
    node_to_active_idx[src] = -1;
    dcell.sorted_objects_valid = false;
    dcell.f2h_one_to_many_index.valid = false;
    tn_object_index_valid = false;
    maintenance_pending = true;
    lazyupdated = true;
}

void TNR_2hop::Apply_Index_Updates() {
    // Dynamic object changes invalidate the affected Cell Hub->Object index
    // and the global TN->Object index. Lazy mode keeps these indexes marked
    // stale; query remains exact via the pairwise and Cell-based fallbacks.
    // This call also checks whether the Cell layout crossed a split/merge threshold.
    bool need_maintenance = false;
    const int split_threshold = Cell_Split_Threshold();
    const int merge_threshold = Cell_Merge_Threshold();

    for (int p_tn : core_tns) {
        if (!global_is_primary_tn[p_tn]) continue;
        const int current_obj_count = static_cast<int>(rt_cells[p_tn].active_objects.size());
        if (current_obj_count >= split_threshold) {
            need_maintenance = true;
            break;
        }
        if (current_obj_count < merge_threshold && static_tree_parent[p_tn] != 0) {
            need_maintenance = true;
            break;
        }
    }

    if (need_maintenance) {
        maintenance_pending = true;
        lazyupdated = true;
    }

    // In lazy mode, invalidated object indexes remain stale and queries use
    // the exact fallback path until a rebuild. In eager mode, any pending
    // invalidation must rebuild even when no Cell threshold was crossed.
    if (eager_maintenance && maintenance_pending) {
        Rebuild_Dynamic_Index_From_Objects();
    }
}

void TNR_2hop::Rebuild_Dynamic_Index_From_Objects() {
    auto start = steady_clock::now();
    maintenance_pending = false;
    lazyupdated = false;

    // ------------------------------------------------------------
    // 1) Recompute object density on the immutable full ET.
    // ------------------------------------------------------------
    local_down_obj_count.assign(node_nums + 1, 0);
    for (int u = 1; u <= node_nums; ++u) {
        if (nodes[u].is_Object) local_down_obj_count[u] = 1;
    }
    for (int rank = 1; rank <= node_nums; ++rank) {
        int u = node_order[rank];
        int p = static_tree_parent[u];
        if (p != 0) local_down_obj_count[p] += local_down_obj_count[u];
    }

    // ------------------------------------------------------------
    // 2) Re-select TNs using the same broad policy as the original code:
    //    high object density, inherited boundary pressure, or ET leaf.
    //    Root TNs are preserved so Merge never destroys a forest root.
    // ------------------------------------------------------------
    vector<bool> old_is_tn = global_is_tn;
    global_is_tn.assign(node_nums + 1, false);
    tn_infection_count.assign(node_nums + 1, 0);

    for (int rank = 1; rank <= node_nums; ++rank) {
        int u = node_order[rank];
        bool is_tn = false;

        const int tn_threshold = Cell_Split_Threshold();
        if (local_down_obj_count[u] >= tn_threshold) {
            is_tn = true;
        } else if (tn_infection_count[u] > 0) {
            is_tn = true;
        } else if (bag_nodes[u].empty()) {
            is_tn = true;
        } else if (old_is_tn[u] && static_tree_parent[u] == 0) {
            is_tn = true;
        }

        global_is_tn[u] = is_tn;
        if (is_tn) {
            for (int w : bag_nodes[u]) {
                ++tn_infection_count[w];
            }
        }
    }

    // ------------------------------------------------------------
    // 3) Rebuild the Cell forest from static ET edges.  TN nodes are
    //    separators: traversal stops at the first descendant TN.
    // ------------------------------------------------------------
    Assemble_Global_Cells();

    // ------------------------------------------------------------
    // 4) Recompute only dynamic structures.  Static bag membership and
    //    bag weights remain untouched and are reused.
    // ------------------------------------------------------------
    compute_distances_and_core_network();
    build_Fast2Hop();
    build_TN_Network();

    // Recreate all query-side dense state so no old dense TN id, Cell pointer,
    // active-object index, or matrix pointer can survive the maintenance swap.
    tn_states.assign(core_tns.size(), TNDenseState());
    dist.assign(node_nums + 1, INT_MAX);
    version.assign(node_nums + 1, 0);
    top_k_pos.assign(node_nums + 1, -1);
    current_version = 0;
    cur_src_ver = 0;
    recycled_dense_ids = std::queue<int>();

    init_cell_runtime_states();

    if (global_top_k_heap == nullptr) {
        global_top_k_heap = new TopKNode[LIMIT_K];
    }
    if (global_batch_targets == nullptr) {
        global_batch_targets = new int[node_nums + 1];
    }
    if (global_exact_dists == nullptr) {
        global_exact_dists = new int[node_nums + 1];
    }
    if (global_lca_ranks == nullptr) {
        global_lca_ranks = new int[node_nums + 1];
    }
    if (global_dijkstra_heap != nullptr) {
        delete[] global_dijkstra_heap;
    }
    global_dijkstra_heap = new HeapNode[max_pq_capacity];

    cout << " [Maintenance] Dynamic TN/Cell index rebuilt safely in "
         << duration_cast<milliseconds>(steady_clock::now() - start).count()
         << " ms" << endl;
}

int TNR_2hop::Find_Best_Split_Point(int p_tn) {
    int best_node = -1;
    int min_diff = 1000000000;

    auto dfs = [&](auto& self, int u) -> int {
        int current_count = nodes[u].is_Object ? 1 : 0;

        for (int child : static_tree_children[u]) {
            if (node_to_primary_tn[child] == p_tn) {
                current_count += self(self, child);
            }
        }

        if (u != p_tn) {
            const int target_objects = CellTargetObjectCount();
            int diff = std::abs(current_count - target_objects);
            if (diff < min_diff) {
                min_diff = diff;
                best_node = u;
            }
        }
        return current_count;
    };

    dfs(dfs, p_tn);

    if (best_node == -1 && !static_tree_children[p_tn].empty()) {
        best_node = static_tree_children[p_tn].front();
    }

    return best_node;
}

vector<int> TNR_2hop::Reassemble_Dirty_Cells(const unordered_set<int>& dirty_p_tns) {
    vector<int> dirty_nodes;

    for (int p_tn : dirty_p_tns) {
        if (!global_is_primary_tn[p_tn]) continue;

        Cell& cell = tn_cells[p_tn];
        cell.primary_tn = p_tn;
        cell.is_valid = true;

        vector<int> locals;
        queue<int> q;
        for (int child : static_tree_children[p_tn]) q.push(child);


        while (!q.empty()) {
            int curr = q.front(); q.pop();
            if (!global_is_tn[curr]) {
                locals.push_back(curr);
                node_to_primary_tn[curr] = p_tn;
                for (int child : static_tree_children[curr]) q.push(child);
            } else {
                // A node promoted from this Cell becomes a TN, not a local
                // vertex.  Its former Cell ownership/runtime row must not
                // survive the split.
                node_to_primary_tn[curr] = 0;
                node_to_cell_row[curr] = -1;
                node_to_active_idx[curr] = -1;
            }
        }
        // Preserve the Cell invariant required by Recompute_Dirty_Cells:
        // every bag dependency is processed before the dependent row.
        std::sort(locals.begin(), locals.end(), [&](int a, int b) {
            return order_find[a] > order_find[b];
        });
        cell.local_nodes = locals;
        dirty_nodes.insert(dirty_nodes.end(), locals.begin(), locals.end());

        unordered_set<int> unique_tns;
        unique_tns.insert(p_tn);
        for (int w : bag_nodes[p_tn]) {
            if (global_is_tn[w]) unique_tns.insert(w);
        }
        for (int v : locals) {
            for (int w : bag_nodes[v]) {
                if (global_is_tn[w]) unique_tns.insert(w);
            }
        }
        cell.access_tns.assign(unique_tns.begin(), unique_tns.end());
        std::sort(cell.access_tns.begin(), cell.access_tns.end());
    }

    return dirty_nodes;
}

void TNR_2hop::Update_Dirty_Core_Network(const unordered_set<int>& dirty_p_tns) {
    for (int u : dirty_p_tns) {
        if (u >= core_network.size()) continue;
        for (auto& edge : core_network[u]) {
            int v = edge.first;
            if (dirty_p_tns.find(v) == dirty_p_tns.end() && v < core_network.size()) {
                auto& v_edges = core_network[v];
                v_edges.erase(std::remove_if(v_edges.begin(), v_edges.end(), [u](const pair<int, int>& e) {
                    return e.first == u;
                }), v_edges.end());
            }
        }
        core_network[u].clear();
    }

    for (int u : dirty_p_tns) {
        if (!global_is_primary_tn[u]) continue;

        unordered_map<int, int> best_edges;

        for (size_t i = 0; i < bag_nodes[u].size(); ++i) {
            int x = bag_nodes[u][i];
            int w_ux = bag_weights[u][i];

            if (w_ux >= SAFE_INF) continue;

            if (global_is_tn[x]) {
                if (best_edges.find(x) == best_edges.end() || w_ux < best_edges[x]) {
                    best_edges[x] = w_ux;
                }
            } else {
                int p_tn_x = node_to_primary_tn[x];
                if (p_tn_x != 0 && global_is_primary_tn[p_tn_x]) {
                    Cell& c = tn_cells[p_tn_x];

                    int row_x = -1;
                    for (size_t r = 0; r < c.local_nodes.size(); ++r) {
                        if (c.local_nodes[r] == x) { row_x = r; break; }
                    }

                    if (row_x != -1) {
                        for (size_t col = 0; col < c.access_tns.size(); ++col) {
                            int target_tn = c.access_tns[col];
                            int dist_x_T = c.tn_to_node_dists[col * c.local_nodes.size() + row_x];
                            if (dist_x_T < SAFE_INF) {
                                int path = safe_add(w_ux, dist_x_T);
                                if (best_edges.find(target_tn) == best_edges.end() || path < best_edges[target_tn]) {
                                    best_edges[target_tn] = path;
                                }
                            }
                        }
                    }
                }
            }
        }

        for (auto& edge : best_edges) {
            core_network[u].push_back({edge.first, edge.second});
        }
    }

    for (int u : dirty_p_tns) {
        if (!global_is_primary_tn[u]) continue;
        for (auto& edge : core_network[u]) {
            int v = edge.first;
            int w = edge.second;
            if (dirty_p_tns.find(v) == dirty_p_tns.end() && v < core_network.size()) {
                core_network[v].push_back({u, w});
            }
        }
    }
}

void TNR_2hop::Update_Dirty_Fast2Hop(const unordered_set<int>& dirty_p_tns, const vector<int>& dirty_nodes) {
    vector<int> temp_anc_dist(node_nums + 1, INT_MAX);
    vector<int> touched_anc;

    for (int u : dirty_nodes) {
        local_ancestors[u].clear();
        if (!global_is_tn[u]) {
            local_ancestors[u].push_back({0, u});
        }
    }

    for (int p_tn : dirty_p_tns) {
        if (global_is_tn[p_tn]) {
            local_ancestors[p_tn].clear();
            if (p_tn < (int)f2h.dynamic_labels.size()) {
                f2h.dynamic_labels[p_tn].clear();
            }
        }
    }

    for (int u : dirty_nodes) {
        if (global_is_tn[u]) continue;

        for (size_t i = 0; i < bag_nodes[u].size(); ++i) {
            int v = bag_nodes[u][i];
            int weight_uv = bag_weights[u][i];

            if (weight_uv >= SAFE_INF) continue;

            if (!global_is_tn[v]) {
                for (auto& anc : local_ancestors[v]) {
                    int hub = anc.second;
                    int new_dist = safe_add(weight_uv, anc.first);

                    if (temp_anc_dist[hub] == INT_MAX) touched_anc.push_back(hub);
                    if (new_dist < temp_anc_dist[hub]) {
                        temp_anc_dist[hub] = new_dist;
                    }
                }
            }
        }

        vector<Fast2Hop::Label> u_labels;
        for (auto& anc : local_ancestors[u]) {
            u_labels.push_back({static_cast<uint32_t>(order_find[anc.second]), anc.first});
        }

        std::sort(touched_anc.begin(), touched_anc.end(), [&](int a, int b) {
            return order_find[a] > order_find[b];
        });

        for (int hub : touched_anc) {
            int cand_dist = temp_anc_dist[hub];
            temp_anc_dist[hub] = INT_MAX;

            if (hub == u) continue;

            uint32_t hub_rank = static_cast<uint32_t>(order_find[hub]);

            bool pruned = false;
            if (hub < (int)f2h.dynamic_labels.size()) {
                const auto& hub_labels = f2h.dynamic_labels[hub];
                size_t i = 0, j = 0;
                size_t u_len = u_labels.size();
                size_t hub_len = hub_labels.size();

                while (i < u_len && j < hub_len) {
                    uint32_t r1 = u_labels[i].rank;
                    uint32_t r2 = hub_labels[j].rank;

                    if (r1 == r2) {
                        int path_d = u_labels[i].dist + hub_labels[j].dist;
                        if (path_d <= cand_dist) {
                            pruned = true;
                            break;
                        }
                        i++; j++;
                    } else {
                        i += (r1 > r2);
                        j += (r1 < r2);
                    }
                }
            }

            if (!pruned) {
                u_labels.push_back({hub_rank, cand_dist});
                local_ancestors[u].push_back({cand_dist, hub});
            }
        }
        touched_anc.clear();

        f2h.update_local_node(u, local_ancestors[u], order_find);
    }

    // 增量更新完成后同步更新连续内存
    f2h.flatten();
}

void TNR_2hop::Recompute_Dirty_Cells(const unordered_set<int>& dirty_p_tns, const vector<int>& dirty_nodes) {

    for (int p_tn : dirty_p_tns) {
        if (!global_is_primary_tn[p_tn]) continue;

        Cell& cell = tn_cells[p_tn];
        int N = cell.local_nodes.size();
        int K = cell.access_tns.size();


        cell.tn_to_node_dists.assign(static_cast<size_t>(K) * N, SAFE_INF);
        cell.min_to_access_dist.assign(K, SAFE_INF);


        for (int row = 0; row < N; ++row) {
            node_to_cell_row[cell.local_nodes[row]] = row;
        }


        for (int row_v = 0; row_v < N; ++row_v) {
            int v = cell.local_nodes[row_v];

            for (int col = 0; col < K; ++col) {
                int target_tn = cell.access_tns[col];
                int min_d = SAFE_INF;


                for (size_t i = 0; i < bag_nodes[v].size(); ++i) {
                    int w = bag_nodes[v][i];
                    int w_weight = bag_weights[v][i];

                    if (w_weight >= SAFE_INF) continue;

                    if (w == target_tn) {

                        if (w_weight < min_d) {
                            min_d = w_weight;
                        }
                    }
                    else if (node_to_primary_tn[w] == p_tn) {

                        int row_w = node_to_cell_row[w];
                        int path = safe_add(w_weight, cell.tn_to_node_dists[col * N + row_w]);
                        if (path < min_d) {
                            min_d = path;
                        }
                    }
                }

                cell.tn_to_node_dists[col * N + row_v] = min_d;
                if (min_d < cell.min_to_access_dist[col]) {
                    cell.min_to_access_dist[col] = min_d;
                }
            }
        }
    }
}

void TNR_2hop::Update_TN_Network(const unordered_set<int>& dirty_p_tns) {

    for (int u : dirty_p_tns) {
        if (global_is_primary_tn[u] && global_to_dense[u] == -1) {
            if (!recycled_dense_ids.empty()) {

                int reuse_id = recycled_dense_ids.front();
                recycled_dense_ids.pop();
                global_to_dense[u] = reuse_id;
                dense_to_global[reuse_id] = u;
                tn_states[reuse_id] = TNDenseState();
            } else {

                int new_id = dense_to_global.size();
                global_to_dense[u] = new_id;
                dense_to_global.push_back(u);
                tn_states.push_back(TNDenseState());
            }
        }
    }

    for (int u : dirty_p_tns) {
        if (!global_is_primary_tn[u] && global_to_dense[u] != -1) {
            int old_dense = global_to_dense[u];
            recycled_dense_ids.push(old_dense);
            global_to_dense[u] = -1;
            dense_to_global[old_dense] = 0;
            tn_states[old_dense] = TNDenseState();
        }
    }

    for (int u : dirty_p_tns) {
        if (!global_is_primary_tn[u]) continue;

        vector<pair<int, int>> valid_edges;
        unordered_map<int, int> dense_tn_weight;

        for (auto& edge : core_network[u]) {
            dense_tn_weight[edge.first] = edge.second;
        }


        for (auto& edge : core_network[u]) {
            int v = edge.first;
            int weight_uv = edge.second;
            bool is_redundant = false;

            for (auto& e_wv : core_network[v]) {
                int w = e_wv.first;
                if (w == u) continue;

                int weight_wv = e_wv.second;

                if (dense_tn_weight.find(w) != dense_tn_weight.end()) {
                    int weight_uw = dense_tn_weight[w];
                    if (weight_uw + weight_wv < weight_uv) {
                        is_redundant = true;
                        break;
                    } else if (weight_uw + weight_wv == weight_uv) {

                        if (w < u && w < v) { is_redundant = true; break; }
                    }
                }
            }

            if (!is_redundant) {
                valid_edges.push_back(edge);
            }
        }

        core_network[u] = valid_edges;

        std::sort(core_network[u].begin(), core_network[u].end(),
                  [](const pair<int, int>& a, const pair<int, int>& b) {
                      return a.second < b.second;
                  });
    }

    size_t directed_core_edges = 0;
    size_t max_access_count = 0;
    for (int u : core_tns) {
        if (!global_is_primary_tn[u]) continue;
        directed_core_edges += core_network[u].size();
        max_access_count = std::max(max_access_count, tn_cells[u].access_tns.size());
    }
    max_pq_capacity = std::max<size_t>(1, directed_core_edges + max_access_count);
}

void TNR_2hop::Update_Cell_Runtime_States(const unordered_set<int>& dirty_p_tns) {

    if (dense_to_global.size() > tn_connected_cells.size()) {
        tn_connected_cells.resize(dense_to_global.size() + 100);
    }

    for (int p_tn : dirty_p_tns) {


        if (!global_is_primary_tn[p_tn]) {
            DenseCell& old_cell = rt_cells[p_tn];


            for (int dense_acc : old_cell.dense_access_tns) {
                auto& conn_list = tn_connected_cells[dense_acc];
                conn_list.erase(std::remove_if(conn_list.begin(), conn_list.end(),
                                               [p_tn](const CellConn& c){ return c.rt_cell_idx == p_tn; }), conn_list.end());
            }


            old_cell.active_objects.clear();
            old_cell.dense_access_tns.clear();
            old_cell.global_p_tn = 0;
            old_cell.tn_to_node_dists = nullptr;
            old_cell.min_access_dist = nullptr;
            old_cell.N = 0;

            for (int v : tn_cells[p_tn].local_nodes) {
                // A merged Cell may already have reassigned this vertex to
                // its surviving parent in the same update batch.
                if (node_to_primary_tn[v] == p_tn) {
                    node_to_cell_row[v] = -1;
                    node_to_active_idx[v] = -1;
                }
            }
            continue;
        }


        DenseCell& dcell = rt_cells[p_tn];
        Cell& cell = tn_cells[p_tn];


        for (int dense_acc : dcell.dense_access_tns) {
            auto& conn_list = tn_connected_cells[dense_acc];
            conn_list.erase(std::remove_if(conn_list.begin(), conn_list.end(),
                                           [p_tn](const CellConn& c){ return c.rt_cell_idx == p_tn; }), conn_list.end());
        }


        dcell.global_p_tn = p_tn;
        dcell.N = static_cast<int>(cell.local_nodes.size());
        dcell.tn_to_node_dists = cell.tn_to_node_dists.data();
        dcell.min_access_dist = cell.min_to_access_dist.data();

        dcell.active_objects.clear();
        dcell.sorted_object_offsets.clear();
        dcell.sorted_object_rows.clear();
        dcell.sorted_objects_valid = false;
        for (int row = 0; row < dcell.N; ++row) {
            int v = cell.local_nodes[row];
            node_to_cell_row[v] = row;

            if (nodes[v].is_Object) {
                node_to_active_idx[v] = dcell.active_objects.size();
                dcell.active_objects.push_back({row, v});
            } else {
                node_to_active_idx[v] = -1;
            }
        }

        int K = cell.access_tns.size();
        dcell.dense_access_tns.resize(K);
        for (int col = 0; col < K; ++col) {
            int acc_tn = cell.access_tns[col];
            int dense_acc = global_to_dense[acc_tn];
            dcell.dense_access_tns[col] = dense_acc;

            tn_connected_cells[dense_acc].push_back({p_tn, col});
        }
        build_sorted_cell_object_index(dcell);
        build_local_f2h_cell_index(dcell);
    }

    // Dirty-cell maintenance changes the exact TN->Object projection.
    // Rebuild it at maintenance granularity, never on every O(1) object update.
    build_TN_Object_Index();
}


void TNR_2hop::compute_f2h_levels() {
    node_f2h_level.assign(node_nums + 1, 0);
    for (int rank = 1; rank <= node_nums; ++rank) {
        int u = node_order[rank];
        if (!global_is_tn[u]) {
            node_f2h_level[u] = 0;
        } else {
            int max_child_level = -1;
            for (int child : static_tree_children[u]) {
                if (node_f2h_level[child] > max_child_level) {
                    max_child_level = node_f2h_level[child];
                }
            }
            node_f2h_level[u] = max_child_level + 1;
        }
    }
}


double TNR_2hop::run_clean_query_benchmark(
        const vector<pair<int, int>>& queries, int warmup) {
    // This method intentionally never enables benchmark_collection_enabled.
    // Therefore query-level counters and fine-grained phase timers stay disabled.
    benchmark_collection_enabled = false;
    if (queries.empty()) {
        last_clean_query_count = 0;
        last_clean_avg_query_us = 0.0;
        return 0.0;
    }

    const int n = static_cast<int>(queries.size());
    if (warmup < 0) warmup = std::min(1000, n / 10);
    // Keep at least one query in the measured sample, even for small input files.
    warmup = std::max(0, std::min(warmup, n - 1));

    vector<pair<int, int>> answers;
    answers.reserve(std::max(1, MAXk));

    // Warmup uses only actual query records loaded from the user's query file.
    for (int i = 0; i < warmup; ++i) {
        quickQueryKNN(queries[static_cast<size_t>(i)].first,
                      queries[static_cast<size_t>(i)].second, answers);
    }

    double total_us = 0.0;
    long long measured = 0;
    for (int i = warmup; i < n; ++i) {
        const auto begin = steady_clock::now();
        quickQueryKNN(queries[static_cast<size_t>(i)].first,
                      queries[static_cast<size_t>(i)].second, answers);
        total_us += duration<double, std::micro>(steady_clock::now() - begin).count();
        ++measured;
    }

    last_clean_query_count = measured;
    last_clean_avg_query_us = measured > 0 ? total_us / static_cast<double>(measured) : 0.0;
    return last_clean_avg_query_us;
}

void TNR_2hop::print_statistics() {
    int total_tns = 0;
    int tns_with_nonempty_cell = 0;
    int total_primary_cells = 0;
    long long total_access_tns_count = 0;

    // 遍历所有 Core TN
    for (int u : core_tns) {
        if (global_is_tn[u]) {
            total_tns++;
        }

        if (global_is_primary_tn[u] && tn_cells[u].is_valid) {
            total_primary_cells++;

            // 1. 下方拥有非空 Local Nodes (子节点/内部节点) 的 Cell 数量
            if (!tn_cells[u].local_nodes.empty()) {
                tns_with_nonempty_cell++;
            }

            // 2. 累加每个 Cell 对应的 Access TN 数量
            total_access_tns_count += tn_cells[u].access_tns.size();
        }
    }

    // 计算平均 Access TN 数量
    double avg_access_tns = (total_primary_cells > 0)
                            ? static_cast<double>(total_access_tns_count) / total_primary_cells
                            : 0.0;

    cout << "\n================ [TNR & Cell 性能统计数据] ================" << endl;
    cout << "1. TN 总数量 (Total Core TNs): " << total_tns << endl;
    cout << "2. 有有效 Cell 的 Primary TN 数量: " << total_primary_cells << endl;
    cout << "   └─ 下方含有非空 Cell 节点 (local_nodes > 0) 的 TN 数量: " << tns_with_nonempty_cell << endl;
    cout << "3. 每个 Cell 平均对应的 Access TN 数量: " << fixed << setprecision(2) << avg_access_tns << endl;
    cout << "   └─ 所有 Cell 的 Access TN 关联总数: " << total_access_tns_count << endl;
    cout << " V7 TN->Object Index: raw=" << tn_object_index_raw_entries
         << ", unique=" << tn_object_index_entries
         << ", duplicates_removed=" << tn_object_index_duplicates_removed << endl;
    cout << "===========================================================\n" << endl;
}
