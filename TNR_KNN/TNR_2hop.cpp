#include <algorithm>
#include <fstream>
#include <mutex>
#include "TNR_2hop.h"
#include "unordered_set"
#include <vector>
#include <queue>
#include <iostream>
#include <utility>
#include <climits>

namespace fs = std::filesystem;

void TNR_2hop::contraction_near_find_down_new(int src) {
    int obj_count = 0;
    bool has_tn_child = false;

    if (nodes[src].is_Object) {
        obj_count = 1;
    }


    for (int v : tree_children[src]) {
        if (tnr_index[v].isTN) {
            has_tn_child = true;
        } else {
            obj_count += tnr_index[v].down_obj_count;
        }
    }

    if (has_tn_child || obj_count >= MAXk * y) {
        tnr_index[src].isTN = true;
        tnr_index[src].down_obj_count = 0;
    } else {
        tnr_index[src].down_obj_count = obj_count;
    }
}

void TNR_2hop::contraction_near_find_up_new(int src) {
    auto& acc_list = tnr_index[src].AccessNode;

    if (global_is_tn[src]) {
        acc_list.push_back({0, src});
        local_ancestors[src].push_back({0, src});
        return;
    }

    local_ancestors[src].push_back({0, src});

    for (auto& edge : bag_nodes_up[src]) {
        int v = edge.first;
        int dist_uv = edge.second;

        for (auto& acc : tnr_index[v].AccessNode) {
            int new_dist = safe_add(dist_uv, acc.dist);
            if (new_dist < temp_dist_acc[acc.tn_id]) {
                if (temp_dist_acc[acc.tn_id] == INT_MAX) {
                    touched_acc.push_back(acc.tn_id);
                }
                temp_dist_acc[acc.tn_id] = new_dist;
            }
        }

        if (!global_is_tn[v]) {
            for (auto& anc : local_ancestors[v]) {
                int new_dist = safe_add(dist_uv, anc.first);
                int hub = anc.second;
                if (new_dist < temp_dist_anc[hub]) {
                    if (temp_dist_anc[hub] == INT_MAX) {
                        touched_anc.push_back(hub);
                    }
                    temp_dist_anc[hub] = new_dist;
                }
            }
        }
    }

    acc_list.reserve(touched_acc.size());
    for (int tn_id : touched_acc) {
        acc_list.push_back({temp_dist_acc[tn_id], tn_id});
        temp_dist_acc[tn_id] = INT_MAX;
    }
    touched_acc.clear();

    for (int hub : touched_anc) {
        local_ancestors[src].push_back({temp_dist_anc[hub], hub});
        temp_dist_anc[hub] = INT_MAX;
    }
    touched_anc.clear();

    sort(acc_list.begin(), acc_list.end(), [](const AccessEntry& a, const AccessEntry& b) {
        return a.tn_id < b.tn_id;
    });
    sort(local_ancestors[src].begin(), local_ancestors[src].end(), [](const pair<int,int>& a, const pair<int,int>& b) {
        return a.second < b.second;
    });

    acc_list.shrink_to_fit();
    local_ancestors[src].shrink_to_fit();
}

void TNR_2hop::build_TN_Network() {
    tn_network.assign(node_nums + 1, vector<pair<int, int>>());
    vector<vector<pair<int, int>>> temp_adj(node_nums + 1);

    for (int u : current_tns) {
        for (auto& edge : bag_nodes_up[u]) {
            int v = edge.first;
            if (tnr_index[v].isTN) {
                temp_adj[u].push_back({v, edge.second});
                temp_adj[v].push_back({u, edge.second});
            }
        }
    }

    vector<int> dense_tn_weight(node_nums + 1, INT_MAX);
    int removed_edges = 0;
    int total_edges = 0;

    for (int u : current_tns) {
        for (auto& edge : temp_adj[u]) {
            dense_tn_weight[edge.first] = edge.second;
        }

        for (auto& edge : temp_adj[u]) {
            int v = edge.first;
            int weight_uv = edge.second;

            if (u >= v) continue;
            total_edges++;

            bool is_redundant = false;

            for (auto& e_wv : temp_adj[v]) {
                int w = e_wv.first;
                if (w == u) continue;

                int weight_wv = e_wv.second;
                int weight_uw = dense_tn_weight[w]; 

                if (weight_uw != INT_MAX && weight_uw + weight_wv <= weight_uv) {
                    is_redundant = true;
                    break;
                }
            }

            if (!is_redundant) {
                tn_network[u].push_back({v, weight_uv});
                tn_network[v].push_back({u, weight_uv});
            } else {
                removed_edges++;
            }
        }

        for (auto& edge : temp_adj[u]) {
            dense_tn_weight[edge.first] = INT_MAX;
        }
    }
    size_t tn_edge_count = 0;
    for (int u : current_tns) {
        tn_edge_count += tn_network[u].size();
    }
    max_pq_capacity = tn_edge_count + LIMIT_ACCESS + node_nums;
}



struct VectorHash {
    size_t operator()(const std::vector<int>& v) const {
        size_t seed = v.size();
        for (auto& i : v) {
            seed ^= i + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

void TNR_2hop::build_super_access_cells() {
    node_to_cell_id.assign(node_nums + 1, -1);
    node_to_cell_row.assign(node_nums + 1, -1);
    tn_connected_cells.resize(node_nums + 1);
    node_to_active_idx.assign(node_nums + 1, -1);
    unordered_map<vector<int>, vector<int>, VectorHash> access_set_map;

    for (int i = 1; i <= node_nums; ++i) {
        if (tnr_index[i].isTN) continue;

        vector<int> tn_ids;
        for (const auto& acc : tnr_index[i].AccessNode) {
            tn_ids.push_back(acc.tn_id);
        }
        sort(tn_ids.begin(), tn_ids.end());
        if (!tn_ids.empty()) {
            access_set_map[tn_ids].push_back(i);
        }
    }

    int cell_counter = 0;
    global_cells.clear();
    global_cells.reserve(access_set_map.size());

    for (const auto& pair : access_set_map) {
        const vector<int>& shared_tns = pair.first;
        const vector<int>& nodes_in_cell = pair.second;

        Cell cell;
        cell.cell_id = cell_counter++;
        cell.tns = shared_tns;

        int K = shared_tns.size();
        int N = nodes_in_cell.size();

        cell.tn_lower_bounds.assign(K, INT_MAX);
        cell.local_nodes = nodes_in_cell;
        cell.tn_to_node_dists.assign(K * N, INT_MAX);

        for (int col = 0; col < K; ++col) {
            cell.tn_local_idx[shared_tns[col]] = col;
        }

        for (int row = 0; row < N; ++row) {
            cell.active_object_rows.reserve(N);
            int u = nodes_in_cell[row];
            node_to_cell_id[u] = cell.cell_id;
            node_to_cell_row[u] = row;

            for (const auto& acc : tnr_index[u].AccessNode) {
                int col = cell.tn_local_idx[acc.tn_id];
                int dist = acc.dist;

                cell.tn_to_node_dists[col * N + row] = dist;

                if (dist < cell.tn_lower_bounds[col]) {
                    cell.tn_lower_bounds[col] = dist;
                }
            }

            if (nodes[u].is_Object) {
                node_to_active_idx[u] = cell.active_object_rows.size();
                cell.active_object_rows.push_back(row);
            }
        }
        for (int col = 0; col < K; ++col) {
            int tn_id = shared_tns[col];
            tn_connected_cells[tn_id].push_back({cell.cell_id, cell.tn_lower_bounds[col], col});
        }
        global_cells.push_back(std::move(cell));
    }

    long long freed_access_count = 0;
    for (int i = 1; i <= node_nums; ++i) {
        freed_access_count += tnr_index[i].AccessNode.size();
        tnr_index[i].AccessNode.clear();
        tnr_index[i].AccessNode.shrink_to_fit();
    }


    global_cell_active_counts.assign(global_cells.size(), 0);
    int initial_active_cells = 0;
    for (size_t i = 0; i < global_cells.size(); ++i) {
        global_cell_active_counts[i] = global_cells[i].active_object_rows.size();
        if (global_cell_active_counts[i] > 0) {
            initial_active_cells++;
        }
    }
    cell_local_cells.assign(global_cells.size(), vector<int>());
    vector<vector<int>> hub_to_cells(node_nums + 1);

    for (int i = 1; i <= node_nums; ++i) {
        int c_id = node_to_cell_id[i];
        if (c_id == -1) continue;
        for (auto& anc : local_ancestors[i]) {
            hub_to_cells[anc.second].push_back(c_id);
        }
    }

    for (int hub = 1; hub <= node_nums; ++hub) {
        auto& cells = hub_to_cells[hub];
        if (cells.empty()) continue;

        sort(cells.begin(), cells.end());
        cells.erase(unique(cells.begin(), cells.end()), cells.end());

        for (int c1 : cells) {
            for (int c2 : cells) {
                cell_local_cells[c1].push_back(c2);
            }
        }
    }

    for (int i = 0; i < global_cells.size(); ++i) {
        auto& neighbors = cell_local_cells[i];
        sort(neighbors.begin(), neighbors.end());
        neighbors.erase(unique(neighbors.begin(), neighbors.end()), neighbors.end());
    }
}


void TNR_2hop::preprocess_KNN(string s , string o) {

    load_raw_graph(s+"Ein.txt");
    string index_file = s + "BCT_Index.bin";

    if (!load_index(index_file)) {
        build_et();
        save_index(index_file);
    }

    tnr_index.resize(node_nums + 1);
    dist.resize(node_nums + 1, INT_MAX);
    version.resize(node_nums + 1, 0);
    node_states.resize(node_nums + 1);
    local_ancestors.resize(node_nums + 1);
    bag_nodes_up.assign(node_nums + 1, vector<pair<int, int>>());
    bag_nodes_down.assign(node_nums + 1, vector<pair<int, int>>());
    top_k_pos.resize(node_nums + 1 , -1);
    temp_dist_acc.assign(node_nums + 1, INT_MAX);
    temp_dist_anc.assign(node_nums + 1, INT_MAX);
    global_is_tn.assign( node_nums+1 , 0);
    addObject(o);



    for (int u = 1; u <= node_nums; ++u) {
        for (size_t i = 0; i < bag_nodes[u].size(); ++i) {
            int w = bag_nodes[u][i];
            int weight = bag_weights[u][i];
            bag_nodes_up[u].push_back({w, weight});
            bag_nodes_down[w].push_back({u, weight});
        }
    }


    current_tns.clear();
    MAXk = min(MAXk ,node_nums);


    current_tns.clear();

    MAXk = min(MAXk, node_nums);

    for (int i = 1; i <= node_nums; ++i) {
        int u = node_order[i];

        if (i == node_nums) {
            tnr_index[u].isTN = true;
        }

        contraction_near_find_down_new(u);

        if (tnr_index[u].isTN) {
            current_tns.push_back(u);
            global_is_tn[u] = true;
        }
    }

    bag_nodes_down.clear();
    bag_nodes_down.shrink_to_fit();

    local_ancestors.assign(node_nums + 1, vector<pair<int, int>>());

    for (int i = node_nums ; i > 0; --i) {
        int u = node_order[i];
        contraction_near_find_up_new(u);
    }

    f2h.build_index(local_ancestors , order_find);

    build_TN_Network();
    bag_nodes_up.clear();
    bag_nodes_up.shrink_to_fit();
    build_super_access_cells();
    if (global_top_k_dists == nullptr) {
        global_top_k_dists = new int[LIMIT_K];
        global_top_k_ids = new int[LIMIT_K];
        global_batch_targets = new int[node_nums + 1];
        global_exact_dists = new int[node_nums + 1];
        global_dijkstra_heap = new HeapNode[max_pq_capacity];
    }

    bag_weights.clear(); bag_weights.shrink_to_fit();
    bag_nodes.clear(); bag_nodes.shrink_to_fit();
    tree_parent.clear(); tree_parent.shrink_to_fit();
    tree_children.clear(); tree_children.shrink_to_fit();
    node_order.clear(); node_order.shrink_to_fit();
    order_find.clear(); order_find.shrink_to_fit();
    graph.clear(); graph.shrink_to_fit();
    local_ancestors.clear(); local_ancestors.shrink_to_fit();

    temp_dist_acc.clear(); temp_dist_acc.shrink_to_fit();
    temp_dist_anc.clear(); temp_dist_anc.shrink_to_fit();
    touched_acc.clear(); touched_acc.shrink_to_fit();
    touched_anc.clear(); touched_anc.shrink_to_fit();


    global_is_tn.assign(node_nums + 1, false);
    for (int i = 1; i <= node_nums; ++i) {
        if (tnr_index[i].isTN) {
            global_is_tn[i] = true;
        }
    }
    tnr_index.clear();
    tnr_index.shrink_to_fit();

}

 void TNR_2hop::quickQueryKNN(int src, int k ,vector<pair<int, int>> &ans) {
    if (++cur_src_ver == 0) {
        node_states.assign(node_nums + 1, NodeState());
        cur_src_ver = 1;
    }
    if (++current_version == 0) {
        std::fill(version.begin(), version.end(), 0);
        current_version = 1;
    }
    D_k = INT_MAX;
    current_k_count = 0;
    int* top_k_dists = global_top_k_dists;
    int* top_k_ids = global_top_k_ids;
    auto evaluate_obj = [&](int obj_id, int path_d) {
        if (current_k_count == k && path_d >= D_k) return;
        if (version[obj_id] == current_version) {
            if (path_d >= dist[obj_id]) return;
            dist[obj_id] = path_d;
            int pos = top_k_pos[obj_id];
            if (pos != -1) {
                top_k_dists[pos] = path_d;
                while (true) {
                    int left = 2 * pos + 1;
                    int right = 2 * pos + 2;
                    int largest = pos;
                    if (left < current_k_count && top_k_dists[left] > top_k_dists[largest]) largest = left;
                    if (right < current_k_count && top_k_dists[right] > top_k_dists[largest]) largest = right;
                    if (largest != pos) {
                        swap(top_k_dists[pos], top_k_dists[largest]);
                        swap(top_k_ids[pos], top_k_ids[largest]);
                        top_k_pos[top_k_ids[pos]] = pos;
                        top_k_pos[top_k_ids[largest]] = largest;
                        pos = largest;
                    } else {
                        break;
                    }
                }
                if (current_k_count == k) D_k = top_k_dists[0];
                return;
            }
        } else {
            version[obj_id] = current_version;
            dist[obj_id] = path_d;
            top_k_pos[obj_id] = -1;
        }
        if (current_k_count < k) {
            int idx = current_k_count++;
            top_k_dists[idx] = path_d;
            top_k_ids[idx] = obj_id;
            top_k_pos[obj_id] = idx;
            while (idx > 0) {
                int parent = (idx - 1) / 2;
                if (top_k_dists[idx] > top_k_dists[parent]) {
                    swap(top_k_dists[idx], top_k_dists[parent]);
                    swap(top_k_ids[idx], top_k_ids[parent]);
                    top_k_pos[top_k_ids[idx]] = idx;
                    top_k_pos[top_k_ids[parent]] = parent;
                    idx = parent;
                } else break;
            }
        } else {
            int old_root_id = top_k_ids[0];
            top_k_pos[old_root_id] = -1;
            top_k_dists[0] = path_d;
            top_k_ids[0] = obj_id;
            top_k_pos[obj_id] = 0;
            int idx = 0;
            while (true) {
                int left = 2 * idx + 1;
                int right = 2 * idx + 2;
                int largest = idx;
                if (left < k && top_k_dists[left] > top_k_dists[largest]) largest = left;
                if (right < k && top_k_dists[right] > top_k_dists[largest]) largest = right;
                if (largest != idx) {
                    swap(top_k_dists[idx], top_k_dists[largest]);
                    swap(top_k_ids[idx], top_k_ids[largest]);
                    top_k_pos[top_k_ids[idx]] = idx;
                    top_k_pos[top_k_ids[largest]] = largest;
                    idx = largest;
                } else break;
            }
        }
        if (current_k_count == k) D_k = top_k_dists[0];
    };
    if (nodes[src].is_Object) evaluate_obj(src, 0);
    int src_cell_id = node_to_cell_id[src];
    if (src_cell_id != -1) {
        int batch_count = 0;
        int* batch_targets = global_batch_targets;
        for (int neighbor_cell_id : cell_local_cells[src_cell_id]) {
            Cell& n_cell = global_cells[neighbor_cell_id];
            for (int row : n_cell.active_object_rows) {
                int obj_id = n_cell.local_nodes[row];
                if (obj_id != src && node_states[obj_id].batch_ver != cur_src_ver) {
                    node_states[obj_id].batch_ver = cur_src_ver;
                    batch_targets[batch_count++] = obj_id;
                }
            }
        }
        if (batch_count > 0) {
            int* exact_dists = global_exact_dists;
            f2h.query_1_to_N_inplace(src, batch_targets, batch_count, exact_dists, D_k);
            for (int i = 0; i < batch_count; ++i) {
                if (exact_dists[i] != INT_MAX) {
                    evaluate_obj(batch_targets[i], exact_dists[i]);
                }
            }
        }
    }
    HeapNode* pq_arr = global_dijkstra_heap;
    int pq_size = 0;
    auto cmp = [](const HeapNode& a, const HeapNode& b) { return a.dist > b.dist; };

    if (global_is_tn[src]) {
        pq_arr[pq_size++] = {0, src};
    } else if (src_cell_id != -1) {
        Cell& cell = global_cells[src_cell_id];
        int row = node_to_cell_row[src];
        int N = cell.local_nodes.size();

        for (int col = 0; col < (int)cell.tns.size(); ++col) {
            int tn_id = cell.tns[col];
            int dist_s_to_tn = cell.tn_to_node_dists[col * N + row];
            if (dist_s_to_tn < D_k) {
                node_states[tn_id].push_ver = cur_src_ver;
                node_states[tn_id].best_tn_dist = dist_s_to_tn;
                pq_arr[pq_size++] = {dist_s_to_tn, tn_id};
            }
        }
    }
    std::make_heap(pq_arr, pq_arr + pq_size, cmp);

    if (current_k_count == k && pq_size > 0 && pq_arr[0].dist >= D_k) {
        pq_size = 0;
    }

    while (pq_size > 0) {
        std::pop_heap(pq_arr, pq_arr + pq_size, cmp);
        auto [curr_dist, u] = pq_arr[--pq_size];

        if (curr_dist >= D_k) break;
        if (node_states[u].settled_ver == cur_src_ver) continue;
        node_states[u].settled_ver = cur_src_ver;

        if (nodes[u].is_Object) evaluate_obj(u, curr_dist);

        for (const auto& conn : tn_connected_cells[u]) {
            int cell_id = conn.cell_id;
            if (curr_dist + conn.lb >= D_k) continue;
            if (global_cell_active_counts[cell_id] == 0) continue;

            Cell& cell = global_cells[cell_id];
            int N = cell.local_nodes.size();
            const int* base_matrix = cell.tn_to_node_dists.data() + (conn.col * N);

            for (int row : cell.active_object_rows) {
                int exact_dist = base_matrix[row];
                int total_path = curr_dist + exact_dist;
                if (total_path >= D_k) continue;
                int obj_id = cell.local_nodes[row];

                evaluate_obj(obj_id, total_path);
            }
        }
        for (auto& edge : tn_network[u]) {
            int v = edge.first;
            int w_uv = edge.second;
            if (node_states[v].settled_ver != cur_src_ver) {
                int next_d = curr_dist + w_uv;
                if (next_d < D_k) {
                    if (node_states[v].push_ver != cur_src_ver || next_d < node_states[v].best_tn_dist) {
                        node_states[v].push_ver = cur_src_ver;
                        node_states[v].best_tn_dist = next_d;
                        pq_arr[pq_size++] = {next_d, v};
                        std::push_heap(pq_arr, pq_arr + pq_size, cmp);
                    }
                }
            }
        }
    }
    for (int i = 0; i < current_k_count; ++i) {
        ans[i] = {top_k_dists[i], top_k_ids[i]};
    }
}

void TNR_2hop::Add_Complete(int src) {
    if (nodes[src].is_Object) return;
    nodes[src].is_Object = true;
    nums_object++;
    if (global_is_tn[src]) return;
    int c_id = node_to_cell_id[src];
    if (c_id == -1) return;
    int c_row = node_to_cell_row[src];
    Cell& cell = global_cells[c_id];
    node_to_active_idx[src] = cell.active_object_rows.size();
    cell.active_object_rows.push_back(c_row);
    global_cell_active_counts[c_id]++;
}

void TNR_2hop::Delete_Complete(int src) {
    if (!nodes[src].is_Object) return;
    nodes[src].is_Object = false;
    nums_object--;

    if (global_is_tn[src]) return;

    int c_id = node_to_cell_id[src];
    if (c_id == -1) return;

    Cell& cell = global_cells[c_id];

    int del_idx = node_to_active_idx[src];
    if (del_idx == -1) return;

    int last_idx = cell.active_object_rows.size() - 1;

    if (del_idx == last_idx) {
        cell.active_object_rows.pop_back();
        node_to_active_idx[src] = -1;
        global_cell_active_counts[c_id]--;
        return;
    }

    int last_row = cell.active_object_rows[last_idx];
    int last_src = cell.local_nodes[last_row];

    cell.active_object_rows[del_idx] = last_row;
    node_to_active_idx[last_src] = del_idx;

    cell.active_object_rows.pop_back();
    node_to_active_idx[src] = -1;
    global_cell_active_counts[c_id]--;
}




void TNR_2hop::save_tnr_index(const string& filename) {
    ofstream out(filename, ios::binary);
    if (!out) {
        cerr << "Error: Cannot open file for writing TNR index: " << filename << endl;
        return;
    }

    out.write(reinterpret_cast<const char*>(&node_nums), sizeof(node_nums));
    out.write(reinterpret_cast<const char*>(&nums_object), sizeof(nums_object));
    out.write(reinterpret_cast<const char*>(&max_pq_capacity), sizeof(max_pq_capacity));

    auto write_vector = [&](const auto& vec) {
        size_t size = vec.size();
        out.write(reinterpret_cast<const char*>(&size), sizeof(size));
        if (size > 0) {
            out.write(reinterpret_cast<const char*>(vec.data()), size * sizeof(vec[0]));
        }
    };

    size_t sz = global_is_tn.size();
    out.write(reinterpret_cast<const char*>(&sz), sizeof(sz));
    for (bool b : global_is_tn) {
        char c = b ? 1 : 0;
        out.write(&c, 1);
    }

    sz = nodes.size();
    out.write(reinterpret_cast<const char*>(&sz), sizeof(sz));
    for (const auto& n : nodes) {
        char c = n.is_Object ? 1 : 0;
        out.write(&c, 1);
    }

    auto write_2d_vector = [&](const auto& vec2d) {
        size_t outer_size = vec2d.size();
        out.write(reinterpret_cast<const char*>(&outer_size), sizeof(outer_size));
        for (const auto& inner_vec : vec2d) {
            write_vector(inner_vec);
        }
    };

    write_2d_vector(tn_network);
    write_2d_vector(tn_connected_cells);
    write_2d_vector(cell_local_cells);

    sz = global_cells.size();
    out.write(reinterpret_cast<const char*>(&sz), sizeof(sz));
    for (const auto& c : global_cells) {
        out.write(reinterpret_cast<const char*>(&c.cell_id), sizeof(c.cell_id));
        write_vector(c.tns);
        write_vector(c.tn_lower_bounds);
        write_vector(c.local_nodes);
        write_vector(c.tn_to_node_dists);
        write_vector(c.active_object_rows);
    }
    write_vector(node_to_cell_id);
    write_vector(node_to_cell_row);
    write_vector(node_to_active_idx);
    write_vector(global_cell_active_counts);

    write_vector(f2h.flat_labels);
    write_vector(f2h.offsets);

    out.close();
    cout << " TNR Index Saved Successfully to: " << filename << " !\n";
}


bool TNR_2hop::load_tnr_index(const string& filename) {
    ifstream in(filename, ios::binary);
    if (!in) {
        return false;
    }

    in.read(reinterpret_cast<char*>(&node_nums), sizeof(node_nums));
    in.read(reinterpret_cast<char*>(&nums_object), sizeof(nums_object));
    in.read(reinterpret_cast<char*>(&max_pq_capacity), sizeof(max_pq_capacity));

    auto read_vector = [&](auto& vec) {
        size_t size;
        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        vec.resize(size);
        if (size > 0) {
            in.read(reinterpret_cast<char*>(vec.data()), size * sizeof(vec[0]));
        }
    };

    size_t sz;
    in.read(reinterpret_cast<char*>(&sz), sizeof(sz));
    global_is_tn.resize(sz);
    for (size_t i = 0; i < sz; ++i) {
        char c;
        in.read(&c, 1);
        global_is_tn[i] = (c != 0);
    }

    in.read(reinterpret_cast<char*>(&sz), sizeof(sz));
    nodes.resize(sz);
    for (size_t i = 0; i < sz; ++i) {
        char c;
        in.read(&c, 1);
        nodes[i].is_Object = (c != 0);
    }

    auto read_2d_vector = [&](auto& vec2d) {
        size_t outer_size;
        in.read(reinterpret_cast<char*>(&outer_size), sizeof(outer_size));
        vec2d.resize(outer_size);
        for (auto& inner_vec : vec2d) {
            read_vector(inner_vec);
        }
    };

    read_2d_vector(tn_network);
    read_2d_vector(tn_connected_cells);
    read_2d_vector(cell_local_cells);

    in.read(reinterpret_cast<char*>(&sz), sizeof(sz));
    global_cells.resize(sz);
    for (auto& c : global_cells) {
        in.read(reinterpret_cast<char*>(&c.cell_id), sizeof(c.cell_id));
        read_vector(c.tns);
        read_vector(c.tn_lower_bounds);
        read_vector(c.local_nodes);
        read_vector(c.tn_to_node_dists);
        read_vector(c.active_object_rows);
    }

    read_vector(node_to_cell_id);
    read_vector(node_to_cell_row);
    read_vector(node_to_active_idx);
    read_vector(global_cell_active_counts);

    read_vector(f2h.flat_labels);
    read_vector(f2h.offsets);

    in.close();

    if (global_top_k_dists == nullptr) {
        global_top_k_dists = new int[LIMIT_K];
        global_top_k_ids = new int[LIMIT_K];
        global_batch_targets = new int[node_nums + 1];
        global_exact_dists = new int[node_nums + 1];
        global_dijkstra_heap = new HeapNode[max_pq_capacity];
    }

    dist.resize(node_nums + 1, INT_MAX);
    version.resize(node_nums + 1, 0);
    node_states.assign(node_nums + 1, NodeState());
    top_k_pos.resize(node_nums + 1, -1);
    cur_src_ver = 0;
    current_version = 0;

    cout << " TNR Index Loaded Successfully from: " << filename << " !\n";
    return true;
}