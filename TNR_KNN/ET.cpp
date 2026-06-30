#include "ET.h"

void ET::build_elimination_tree() {
    cout << "--- [Phase 2] Constructing Explicit Elimination Tree ---\n";
    for (int u = 1; u <= node_nums; ++u) {
        for (auto& edge : graph[u]) {
            int v = edge.first;
            if (order_find[v] > order_find[u]) {
                bag_nodes[u].push_back(v);
            }
        }
    }
    vector<int> edge_mark(node_nums + 1, 0);
    for (int rank = 1; rank <= node_nums; ++rank) {
        int u = node_order[rank];
        if (bag_nodes[u].empty()) continue;

        int parent = -1;
        int min_rank = SAFE_INF;
        for (int ancestor : bag_nodes[u]) {
            if (order_find[ancestor] < min_rank) {
                min_rank = order_find[ancestor];
                parent = ancestor;
            }
        }

        tree_parent[u] = parent;
        tree_children[parent].push_back(u);

        for (int w : bag_nodes[parent]) {
            edge_mark[w] = u;
        }

        for (int w : bag_nodes[u]) {
            if (w != parent && edge_mark[w] != u) {
                bag_nodes[parent].push_back(w);
                edge_mark[w] = u;
            }
        }
    }
    for (int u = 1; u <= node_nums; ++u) {
        sort(bag_nodes[u].begin(), bag_nodes[u].end(), [&](int a, int b) {
            return order_find[a] < order_find[b];
        });
        bag_nodes[u].erase(unique(bag_nodes[u].begin(), bag_nodes[u].end()), bag_nodes[u].end());
    }
    cout << " Explicit Tree Built! Parent/Children arrays are ready.\n";
}

void ET::compute_bag_distances() {
    cout << "--- [Phase 3] Computing Tree Decomposition Bag Distances ---\n";
    bag_weights.assign(node_nums + 1, vector<int>());
    for (int u = 1; u <= node_nums; ++u) {
        bag_weights[u].assign(bag_nodes[u].size(), SAFE_INF);
        for (size_t i = 0; i < bag_nodes[u].size(); ++i) {
            int w = bag_nodes[u][i];
            if (graph[u].count(w)) {
                bag_weights[u][i] = graph[u][w];
            }
        }
    }

    vector<int> dense_temp_dist(node_nums + 1, SAFE_INF);

    for (int rank = 1; rank <= node_nums; ++rank) {
        int u = node_order[rank];
        const vector<int>& bag = bag_nodes[u];
        const vector<int>& dist_from_u = bag_weights[u]; 

        for (size_t i = 0; i < bag.size(); ++i) {
            int v = bag[i];
            int dist_uv = dist_from_u[i];


            for (size_t k = 0; k < bag_nodes[v].size(); ++k) {
                int w = bag_nodes[v][k];
                dense_temp_dist[w] = bag_weights[v][k];
            }

            for (size_t j = i + 1; j < bag.size(); ++j) {
                int w = bag[j];
                int dist_uw = dist_from_u[j];
                int shortcut_dist = safe_add(dist_uv, dist_uw);

                if (shortcut_dist < dense_temp_dist[w]) {
                    dense_temp_dist[w] = shortcut_dist;
                }
            }


            for (size_t k = 0; k < bag_nodes[v].size(); ++k) {
                int w = bag_nodes[v][k];
                bag_weights[v][k] = dense_temp_dist[w];
                dense_temp_dist[w] = SAFE_INF;
            }
        }
    }
    cout << " Bag Distances Computed!\n";
}

void ET::load_raw_graph(const string& filename) {
    ifstream in(filename);
    assert(in.is_open());
    in >> node_nums >> edge_nums;

    this->nodes.resize(node_nums + 1);
    this->graph.resize(node_nums + 1);

    node_order.resize(node_nums + 1);
    order_find.resize(node_nums + 1);

    tree_parent.assign(node_nums + 1, 0);
    tree_children.resize(node_nums + 1);

    bag_nodes.resize(node_nums + 1);
    bag_weights.resize(node_nums + 1);

    int u, v, w;
    while (in >> u >> v >> w) {
        graph[u][v] = w;
        graph[v][u] = w;
    }
    in.close();
    cout << "Loaded graph: " << node_nums << " nodes.\n";
}

void ET::save_index(const string& filename) {
    cout << "--- Saving BCT Index to Binary File: " << filename << " ---\n";
    ofstream out(filename, ios::binary);
    if (!out) return;

    out.write(reinterpret_cast<const char*>(&node_nums), sizeof(node_nums));
    out.write(reinterpret_cast<const char*>(&edge_nums), sizeof(edge_nums));

    auto write_vector = [&](const vector<int>& vec) {
        size_t size = vec.size();
        out.write(reinterpret_cast<const char*>(&size), sizeof(size));
        if (size > 0) out.write(reinterpret_cast<const char*>(vec.data()), size * sizeof(int));
    };
    auto write_2d_vector = [&](const vector<vector<int>>& vec2d) {
        size_t outer_size = vec2d.size();
        out.write(reinterpret_cast<const char*>(&outer_size), sizeof(outer_size));
        for (const auto& inner_vec : vec2d) write_vector(inner_vec);
    };

    write_vector(node_order);
    write_vector(order_find);
    write_vector(tree_parent);
    write_2d_vector(tree_children);
    write_2d_vector(bag_nodes);
    write_2d_vector(bag_weights);
    out.close();
}

bool ET::load_index(const string& filename) {
    cout << "--- Loading BCT Index from Binary File: " << filename << " ---\n";
    ifstream in(filename, ios::binary);
    if (!in) return false;

    in.read(reinterpret_cast<char*>(&node_nums), sizeof(node_nums));
    in.read(reinterpret_cast<char*>(&edge_nums), sizeof(edge_nums));

    auto read_vector = [&](vector<int>& vec) {
        size_t size;
        in.read(reinterpret_cast<char*>(&size), sizeof(size));
        vec.resize(size);
        if (size > 0) in.read(reinterpret_cast<char*>(vec.data()), size * sizeof(int));
    };
    auto read_2d_vector = [&](vector<vector<int>>& vec2d) {
        size_t outer_size;
        in.read(reinterpret_cast<char*>(&outer_size), sizeof(outer_size));
        vec2d.resize(outer_size);
        for (auto& inner_vec : vec2d) read_vector(inner_vec);
    };

    read_vector(node_order);
    read_vector(order_find);
    read_vector(tree_parent);
    read_2d_vector(tree_children);
    read_2d_vector(bag_nodes);
    read_2d_vector(bag_weights);

    in.close();
    return true;
}


void ET::build_et() {
    auto start = this->get_time();
    generate_level_penalized_greedy_order();
    build_elimination_tree();
    compute_bag_distances();
    cout << " Balanced Cut Tree (BCT) Build Complete! Total time: " << this->get_time() - start << " s\n";
}


void ET::generate_level_penalized_greedy_order() {
    cout << "--- [Phase 1] Building Order (Level-Penalized Greedy) ---\n";

    vector<unordered_set<int>> dynamic_graph(node_nums + 1);
    for (int u = 1; u <= node_nums; ++u) {
        for (auto& edge : graph[u]) {
            if (edge.first != u) {
                dynamic_graph[u].insert(edge.first);
            }
        }
    }

    vector<int> level(node_nums + 1, 0);
    vector<bool> contracted(node_nums + 1, false);

    const int W_EDGE_DIFF = 1000;
    const int W_DEGREE = 100;
    const int W_LEVEL = 500;

    auto calc_score = [&](int u) {
        int deg = dynamic_graph[u].size();
        if (deg == 0) return level[u] * W_LEVEL;

        int new_edges = 0;
        vector<int> nbs(dynamic_graph[u].begin(), dynamic_graph[u].end());

        for (size_t i = 0; i < nbs.size(); ++i) {
            for (size_t j = i + 1; j < nbs.size(); ++j) {
                if (dynamic_graph[nbs[i]].find(nbs[j]) == dynamic_graph[nbs[i]].end()) {
                    new_edges++;
                }
            }
        }

        int edge_diff = new_edges - deg;
        return edge_diff * W_EDGE_DIFF + deg * W_DEGREE + level[u] * W_LEVEL;
    };

    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    for (int u = 1; u <= node_nums; ++u) {
        pq.push({calc_score(u), u});
    }

    int current_rank = 1;

    while (!pq.empty()) {
        auto [score, u] = pq.top();
        pq.pop();

        if (contracted[u]) continue;

        int real_score = calc_score(u);
        if (real_score > score) {
            pq.push({real_score, u});
            continue;
        }

        contracted[u] = true;
        order_find[u] = current_rank;
        node_order[current_rank] = u;
        current_rank++;

        vector<int> neighbors(dynamic_graph[u].begin(), dynamic_graph[u].end());
        for (size_t i = 0; i < neighbors.size(); ++i) {
            int v = neighbors[i];
            dynamic_graph[v].erase(u);
            if (level[v] <= level[u]) {
                level[v] = level[u] + 1;
            }
            for (size_t j = i + 1; j < neighbors.size(); ++j) {
                int w = neighbors[j];
                dynamic_graph[v].insert(w);
                dynamic_graph[w].insert(v);
            }
        }
        dynamic_graph[u].clear();
    }
    cout << " Level-Penalized Greedy Order Generated Successfully!\n";
}