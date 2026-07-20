#include <algorithm>
#include "TNR_2hop.h"
#include <vector>
#include <queue>
#include <iostream>
#include <utility>
#include <climits>

namespace fs = std::filesystem;

// 1. 抽离出的独立 F2H 重构函数
void TNR_2hop::build_Fast2Hop() {
    cout << "--- [Phase 3.5]  ---" << endl;
    compute_f2h_levels();


    local_ancestors.assign(node_nums + 1, vector<pair<int, int>>());
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
                        temp_anc_dist[hub] = new_dist;
                    }
                }
            }
        }
        for (int hub : touched_anc) {
            local_ancestors[u].push_back({temp_anc_dist[hub], hub});
            temp_anc_dist[hub] = INT_MAX;
        }
        touched_anc.clear();
    }

    f2h.build_index(local_ancestors, order_find);

}


void TNR_2hop::build_TN_Network() {
    int num_tns = core_tns.size();
    global_to_dense.assign(node_nums + 1, -1);
    dense_to_global.assign(num_tns, 0);

    for (int i = 0; i < num_tns; ++i) {
        int global_u = core_tns[i];
        global_to_dense[global_u] = i;
        dense_to_global[i] = global_u;
    }

    vector<int> dense_tn_weight(node_nums + 1, INT_MAX);
    size_t csr_edges = 0;

    for (int u : core_tns) {
        vector<pair<int, int>> valid_edges;
        for (auto& edge : core_network[u]) dense_tn_weight[edge.first] = edge.second;

        for (auto& edge : core_network[u]) {
            int v = edge.first;
            int weight_uv = edge.second;
            bool is_redundant = false;

            for (auto& e_wv : core_network[v]) {
                int w = e_wv.first;
                if (w == u) continue;
                if (dense_tn_weight[w] != INT_MAX) {
                    int weight_uw = dense_tn_weight[w];
                    int weight_wv = e_wv.second;
                    if (weight_uw + weight_wv < weight_uv) {
                        is_redundant = true; break;
                    } else if (weight_uw + weight_wv == weight_uv) {
                        if (w < u && w < v) { is_redundant = true; break; }
                    }
                }
            }
            if (!is_redundant) {
                valid_edges.push_back(edge);
                csr_edges++;
            }
        }
        for (auto& edge : core_network[u]) dense_tn_weight[edge.first] = INT_MAX;

        core_network[u] = valid_edges;
        std::sort(core_network[u].begin(), core_network[u].end(),
                  [](const pair<int, int>& a, const pair<int, int>& b) { return a.second < b.second; });
    }

    max_pq_capacity = csr_edges + LIMIT_ACCESS + node_nums;
    cout << " TN Academic Core Network ！" << endl;
}

void TNR_2hop::init_cell_runtime_states() {
    cout << "\n=== [Phase 4]  ===" << endl;

    node_to_cell_row.assign(node_nums + 1, -1);
    node_to_active_idx.assign(node_nums + 1, -1);


    rt_cells.assign(node_nums + 1, DenseCell());


    tn_connected_cells.assign(dense_to_global.size() + 100, vector<CellConn>());

    for (int p_tn : core_tns) {
        if (!global_is_primary_tn[p_tn]) continue;


        DenseCell& dcell = rt_cells[p_tn];
        Cell& cell = tn_cells[p_tn];

        dcell.global_p_tn = p_tn;
        dcell.N = cell.local_nodes.size();
        dcell.tn_to_node_dists = cell.tn_to_node_dists.data();
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
    }
}

void TNR_2hop::preprocess_KNN(string s, string o) {
    cout << "===  start  ===" << endl;
    auto start_time = steady_clock::now();

    load_raw_graph(s + "graph.txt");
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
        global_dijkstra_heap = new HeapNode[max_pq_capacity];
    }

    auto end_time = steady_clock::now();
    cout  << duration_cast<milliseconds>(end_time - start_time).count() << " ms ===" << endl;
}

void TNR_2hop::quickQueryKNN(int src, int k ,vector<pair<int, int>> &ans) {
    if (++cur_src_ver == 0) cur_src_ver = 1;
    if (++current_version == 0) {
        std::fill(version.begin(), version.end(), 0);
        current_version = 1;
    }
    D_k = INT_MAX;
    current_k_count = 0;
    TopKNode* top_k_heap = global_top_k_heap;

    auto evaluate_obj = [&](int obj_id, int path_d) {
        if (current_k_count == k && path_d >= D_k) return;
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

    if (nodes[src].is_Object) evaluate_obj(src, 0);

    int p_tn_src = node_to_primary_tn[src];

    if (p_tn_src != 0) {
        int batch_count = 0;
        int* batch_targets = global_batch_targets;

        DenseCell& dcell = rt_cells[p_tn_src];

        for (const auto& act_obj : dcell.active_objects) {
            int obj_id = act_obj.second;
            if (obj_id != src) {
                batch_targets[batch_count++] = obj_id;
            }
        }

        if (batch_count > 0) {
            int* exact_dists = global_exact_dists;

            f2h.query_1_to_N_inplace(src, batch_targets, batch_count, exact_dists, D_k);

            for (int i = 0; i < batch_count; ++i) {
                if (exact_dists[i] != INT_MAX) evaluate_obj(batch_targets[i], exact_dists[i]);
            }
        }
    }

    HeapNode* pq_arr = global_dijkstra_heap;
    int pq_size = 0;
    auto cmp = [](const HeapNode& a, const HeapNode& b) { return a.dist > b.dist; };

    if (global_is_tn[src]) {
        pq_arr[pq_size++] = {0, global_to_dense[src]};
    } else if (p_tn_src != 0) {
        DenseCell& dcell = rt_cells[p_tn_src];
        int row = node_to_cell_row[src];

        for (int col = 0; col < (int)dcell.dense_access_tns.size(); ++col) {
            int dense_tn = dcell.dense_access_tns[col];
            int dist_s_to_tn = dcell.tn_to_node_dists[col * dcell.N + row];

            if (dist_s_to_tn < SAFE_INF && dist_s_to_tn < D_k) {
                tn_states[dense_tn].push_ver = cur_src_ver;
                tn_states[dense_tn].best_tn_dist = dist_s_to_tn;
                pq_arr[pq_size++] = {dist_s_to_tn, dense_tn};
            }
        }
    }
    std::make_heap(pq_arr, pq_arr + pq_size, cmp);

    if (current_k_count == k && pq_size > 0 && pq_arr[0].dist >= D_k) pq_size = 0;

    while (pq_size > 0  && pq_arr[0].dist < D_k) {
        std::pop_heap(pq_arr, pq_arr + pq_size, cmp);
        auto [curr_dist, dense_u] = pq_arr[--pq_size];

        if (curr_dist >= D_k) break;

        if (tn_states[dense_u].settled_ver == cur_src_ver) continue;
        tn_states[dense_u].settled_ver = cur_src_ver;

        int global_u = dense_to_global[dense_u];
        if (nodes[global_u].is_Object) evaluate_obj(global_u, curr_dist);

        for (const CellConn& conn : tn_connected_cells[dense_u]) {
            DenseCell& dcell = rt_cells[conn.rt_cell_idx];
            if (dcell.active_objects.empty()) continue;

            const int* base_matrix = dcell.tn_to_node_dists + (conn.col * dcell.N);

            for (const auto& act_obj : dcell.active_objects) {
                int exact_dist = base_matrix[act_obj.first];
                if (exact_dist >= SAFE_INF) continue;

                int total_path = curr_dist + exact_dist;
                if (total_path < D_k) {
                    evaluate_obj(act_obj.second, total_path);
                }
            }
        }

        for (const auto& edge : core_network[dense_to_global[dense_u]]) {
            int w_uv = edge.second;
            int dense_v = global_to_dense[edge.first];

            if (curr_dist + w_uv >= D_k) break;

            if (tn_states[dense_v].settled_ver != cur_src_ver) {
                int next_d = curr_dist + w_uv;
                if (tn_states[dense_v].push_ver != cur_src_ver) {
                    tn_states[dense_v].best_tn_dist = INT_MAX;
                }
                if (next_d < tn_states[dense_v].best_tn_dist) {
                    tn_states[dense_v].push_ver = cur_src_ver;
                    tn_states[dense_v].best_tn_dist = next_d;
                    pq_arr[pq_size++] = {next_d, dense_v};
                    std::push_heap(pq_arr, pq_arr + pq_size, cmp);
                }
            }
        }
    }

    for (int i = 0; i < current_k_count; ++i) ans[i] = {top_k_heap[i].dist, top_k_heap[i].id};

    std::sort(ans.begin(), ans.begin() + current_k_count, [](const pair<int, int>& a, const pair<int, int>& b) {
        if (a.first != b.first) return a.first < b.first;
        return a.second < b.second;
    });
}

void TNR_2hop::Add_Complete(int src) {
    if (nodes[src].is_Object) return;
    nodes[src].is_Object = true;
    nums_object++;

    if (global_is_tn[src]) return;
    int p_tn = node_to_primary_tn[src];
    if (p_tn == 0) return;

    int c_row = node_to_cell_row[src];

    DenseCell& dcell = rt_cells[p_tn];
    node_to_active_idx[src] = dcell.active_objects.size();
    dcell.active_objects.push_back({c_row, src});

}

void TNR_2hop::Delete_Complete(int src) {
    if (!nodes[src].is_Object) return;
    nodes[src].is_Object = false;
    nums_object--;

    if (global_is_tn[src]) return;

    int p_tn = node_to_primary_tn[src];
    if (p_tn == 0) return;

    int del_idx = node_to_active_idx[src];
    if (del_idx == -1) return;

    DenseCell& dcell = rt_cells[p_tn];
    int last_idx = dcell.active_objects.size() - 1;

    if (del_idx != last_idx) {
        auto last_obj = dcell.active_objects[last_idx];
        dcell.active_objects[del_idx] = last_obj;
        node_to_active_idx[last_obj.second] = del_idx;
    }
    dcell.active_objects.pop_back();
    node_to_active_idx[src] = -1;
}

void TNR_2hop::Apply_Index_Updates() {
    auto t_start_total = steady_clock::now();
    bool structure_changed = false;
    unordered_set<int> dirty_p_tns;


    for (int p_tn : core_tns) {
        if (!global_is_primary_tn[p_tn]) continue;

        int current_obj_count = rt_cells[p_tn].active_objects.size();


        if (current_obj_count >= MAXk * 2) {
            structure_changed = true;
            int new_tn = Find_Best_Split_Point(p_tn);

            global_is_tn[new_tn] = true;
            global_is_primary_tn[new_tn] = true;
            core_tns.push_back(new_tn);

            dirty_p_tns.insert(p_tn);
            dirty_p_tns.insert(new_tn);
        }

        else if (current_obj_count < MAXk / 2 && static_tree_parent[p_tn] != 0) {
            structure_changed = true;
            global_is_tn[p_tn] = false;
            global_is_primary_tn[p_tn] = false;
            tn_cells[p_tn].is_valid = false;


            auto it = std::find(core_tns.begin(), core_tns.end(), p_tn);
            if (it != core_tns.end()) core_tns.erase(it);


            int parent_tn = static_tree_parent[p_tn];
            while (parent_tn != 0 && !global_is_primary_tn[parent_tn]) {
                parent_tn = static_tree_parent[parent_tn];
            }
            if (parent_tn != 0) dirty_p_tns.insert(parent_tn);
        }
    }


    if (!structure_changed) return;




    vector<int> dirty_nodes = Reassemble_Dirty_Cells(dirty_p_tns);

    std::sort(dirty_nodes.begin(), dirty_nodes.end(), [&](int a, int b) {
        return order_find[a] > order_find[b];
    });

    Recompute_Dirty_Cells(dirty_p_tns, dirty_nodes);

    Update_Dirty_Core_Network(dirty_p_tns);

    lazyupdated = false;

    // 检查脏区内是否存在突破当前 2-hop 缓冲层高的节点
    for (int p_tn : dirty_p_tns) {
        // 防止新生成的 TN ID 越界，增加一次大小校验
        if (p_tn < node_f2h_level.size() && node_f2h_level[p_tn] >= current_f2h_height) {
            lazyupdated = true;
            break;
        }
    }

    if (lazyupdated) {

        current_f2h_height += 10;


        Update_Dirty_Fast2Hop(dirty_p_tns, dirty_nodes);

    } else {

    }



    Update_TN_Network(dirty_p_tns);


    if (global_dijkstra_heap != nullptr) delete[] global_dijkstra_heap;
    global_dijkstra_heap = new HeapNode[max_pq_capacity];
    Update_Cell_Runtime_States(dirty_p_tns);

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
            int diff = std::abs(current_count - MAXk);
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
            }
        }
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

void TNR_2hop::Recompute_Dirty_Cells(const unordered_set<int>& dirty_p_tns, const vector<int>& dirty_nodes) {

    for (int p_tn : dirty_p_tns) {
        if (!global_is_primary_tn[p_tn]) continue;

        Cell& cell = tn_cells[p_tn];
        int N = cell.local_nodes.size();
        int K = cell.access_tns.size();


        cell.tn_to_node_dists.assign(K * N, SAFE_INF);


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
            recycled_dense_ids.push(global_to_dense[u]);
            global_to_dense[u] = -1;

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

                    if (new_dist < temp_anc_dist[hub]) {
                        if (temp_anc_dist[hub] == INT_MAX) touched_anc.push_back(hub);
                        temp_anc_dist[hub] = new_dist;
                    }
                }
            }
        }

        for (int hub : touched_anc) {
            local_ancestors[u].push_back({temp_anc_dist[hub], hub});
            temp_anc_dist[hub] = INT_MAX;
        }
        touched_anc.clear();
    }

    for (int u : dirty_nodes) {
        if (!global_is_tn[u]) {
            f2h.update_local_node(u, local_ancestors[u], order_find);
        } else {
            f2h.dynamic_labels[u].clear();
        }
    }


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

            for (int v : tn_cells[p_tn].local_nodes) {
                node_to_cell_row[v] = -1;
                node_to_active_idx[v] = -1;
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
        dcell.N = cell.local_nodes.size();
        dcell.tn_to_node_dists = cell.tn_to_node_dists.data();

        dcell.active_objects.clear();
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
    }
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