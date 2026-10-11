#include "ET.h"
#include <unordered_set>
#include <algorithm>
#include <queue>

namespace {

// Matrix DP evaluates a node through higher-rank bag vertices.  Keep this
// order as a Cell invariant when a Cell is assembled, so dynamic matrix
// recomputation can use one linear forward scan.
inline void sort_cell_local_nodes_by_descending_rank(vector<int>& local_nodes,
                                                      const vector<int>& order_find) {
    std::sort(local_nodes.begin(), local_nodes.end(),
              [&order_find](int lhs, int rhs) {
                  return order_find[lhs] > order_find[rhs];
              });
}

} // namespace

void ET::initialize_static_tree_structures() {
    static_tree_parent.assign(node_nums + 1, 0);
    static_tree_children.assign(node_nums + 1, vector<int>());
    local_down_obj_count.assign(node_nums + 1, 0);
}


void ET::build_et(int k) {
    this->MAXk = k;
    static_distance_ready = false;
    auto start = this->get_time();

    generate_full_greedy_order();
    build_tree_and_truncate_by_object();
    compute_distances_and_core_network();

    cout << " [ET Phase] Partial ET & Core Network Build Complete! Total time: " << this->get_time() - start << " s\n";
}

void ET::generate_full_greedy_order() {
    cout << "--- [Phase 1] Building Full Order (Greedy Contraction) ---\n";
    vector<unordered_set<int>> dynamic_graph(node_nums + 1);
    for (int u = 1; u <= node_nums; ++u) {
        for (auto& edge : graph[u]) {
            if (edge.first != u) dynamic_graph[u].insert(edge.first);
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
                if (dynamic_graph[nbs[i]].find(nbs[j]) == dynamic_graph[nbs[i]].end()) new_edges++;
            }
        }
        int edge_diff = new_edges - deg;
        return edge_diff * W_EDGE_DIFF + deg * W_DEGREE + level[u] * W_LEVEL;
    };

    priority_queue<pair<int, int>, vector<pair<int, int>>, greater<pair<int, int>>> pq;
    for (int u = 1; u <= node_nums; ++u) pq.push({calc_score(u), u});

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
            if (level[v] <= level[u]) level[v] = level[u] + 1;
            for (size_t j = i + 1; j < neighbors.size(); ++j) {
                int w = neighbors[j];
                dynamic_graph[v].insert(w);
                dynamic_graph[w].insert(v);
            }
        }
        dynamic_graph[u].clear();
    }
}

void ET::compute_distances_and_core_network() {
    cout << "--- [Phase 3] Computing Static CH Distances & Assembling Core Network ---\n";

    // The elimination order and bags are static.  Their edge weights must not
    // depend on the current TN selection; otherwise a later Split/Merge would
    // silently change the meaning of the ET.  Build the bag distances once.
    if (!static_distance_ready) {
        bag_weights.assign(node_nums + 1, vector<int>());

        for (int rank = 1; rank <= node_nums; ++rank) {
            int u = node_order[rank];
            if (bag_nodes[u].empty()) continue;

            bag_weights[u].assign(bag_nodes[u].size(), SAFE_INF);
            for (size_t i = 0; i < bag_nodes[u].size(); ++i) {
                int w = bag_nodes[u][i];
                auto it = graph[u].find(w);
                if (it != graph[u].end()) {
                    bag_weights[u][i] = it->second;
                }
            }
        }

        // Top-down CH closure: when u is contracted, create/update the
        // shortest known distance between each pair of higher-ranked bag
        // vertices.  Importantly, this pass is independent of TN status.
        vector<int> temp_dist(node_nums + 1, SAFE_INF);
        for (int rank = 1; rank <= node_nums; ++rank) {
            int u = node_order[rank];
            const auto& bag = bag_nodes[u];
            if (bag.empty()) continue;

            const auto& dist_from_u = bag_weights[u];
            for (size_t i = 0; i < bag.size(); ++i) {
                int v = bag[i];
                int dist_uv = dist_from_u[i];
                if (dist_uv >= SAFE_INF) continue;

                // Start from v's current upward distances.
                for (size_t q = 0; q < bag_nodes[v].size(); ++q) {
                    temp_dist[bag_nodes[v][q]] = bag_weights[v][q];
                }

                // Relax v -> w through u.  bag[] is sorted by increasing
                // contraction rank, so j > i implies v is lower-ranked than w.
                for (size_t j = i + 1; j < bag.size(); ++j) {
                    int w = bag[j];
                    int shortcut_dist = safe_add(dist_uv, dist_from_u[j]);
                    if (shortcut_dist < temp_dist[w]) {
                        temp_dist[w] = shortcut_dist;
                    }
                }

                for (size_t q = 0; q < bag_nodes[v].size(); ++q) {
                    bag_weights[v][q] = temp_dist[bag_nodes[v][q]];
                    temp_dist[bag_nodes[v][q]] = SAFE_INF;
                }
                for (size_t j = i + 1; j < bag.size(); ++j) {
                    temp_dist[bag[j]] = SAFE_INF;
                }
            }
        }

        static_distance_ready = true;
    }

    // Rebuild only the dynamic projection of the static CH structure onto
    // the current set of Transit Nodes.  This makes Split/Merge safe: no old
    // TN edge or distance survives from a previous TN configuration.
    core_network.assign(node_nums + 1, vector<pair<int, int>>());
    vector<unordered_map<int, int>> tn_edges(node_nums + 1);

    for (int u = 1; u <= node_nums; ++u) {
        if (!global_is_tn[u]) continue;
        for (size_t i = 0; i < bag_nodes[u].size(); ++i) {
            int v = bag_nodes[u][i];
            if (!global_is_tn[v]) continue;
            int w = bag_weights[u][i];
            if (w >= SAFE_INF) continue;
            auto it = tn_edges[u].find(v);
            if (it == tn_edges[u].end() || w < it->second) {
                tn_edges[u][v] = w;
                tn_edges[v][u] = w;
            }
        }
    }

    for (int u : core_tns) {
        for (const auto& edge : tn_edges[u]) {
            core_network[u].push_back(edge);
        }
        std::sort(core_network[u].begin(), core_network[u].end(),
                  [](const pair<int, int>& a, const pair<int, int>& b) {
                      if (a.second != b.second) return a.second < b.second;
                      return a.first < b.first;
                  });
    }

    cout << "--- [Phase 3.5] Extracting Localized Cell Distance Matrix (Top-Down DP) ---\n";
    vector<int> local_node_to_row(node_nums + 1, -1);

    for (int u : core_tns) {
        if (!global_is_primary_tn[u]) continue;
        Cell& cell = tn_cells[u];
        int N = static_cast<int>(cell.local_nodes.size());
        int K = static_cast<int>(cell.access_tns.size());
        cell.tn_to_node_dists.assign(static_cast<size_t>(K) * N, SAFE_INF);
        cell.min_to_access_dist.assign(K, SAFE_INF);
        for (int row = 0; row < N; ++row) {
            local_node_to_row[cell.local_nodes[row]] = row;
        }
    }

    // local_nodes are kept in descending rank, so all local bag dependencies
    // that can be used by v have already been computed when v is processed.
    for (int rank = node_nums; rank >= 1; --rank) {
        int v = node_order[rank];
        if (global_is_tn[v]) continue;

        int p_tn = node_to_primary_tn[v];
        if (p_tn == 0) continue;

        Cell& cell = tn_cells[p_tn];
        int row = local_node_to_row[v];
        if (row < 0) continue;
        int K = static_cast<int>(cell.access_tns.size());
        int N = static_cast<int>(cell.local_nodes.size());

        for (int col = 0; col < K; ++col) {
            int target_tn = cell.access_tns[col];
            int min_d = SAFE_INF;

            for (size_t i = 0; i < bag_nodes[v].size(); ++i) {
                int w = bag_nodes[v][i];
                int w_weight = bag_weights[v][i];
                if (w_weight >= SAFE_INF) continue;

                int d_w_to_T = SAFE_INF;
                if (w == target_tn) {
                    d_w_to_T = 0;
                } else if (!global_is_tn[w]) {
                    int w_p_tn = node_to_primary_tn[w];
                    if (w_p_tn == p_tn) {
                        int w_row = local_node_to_row[w];
                        if (w_row >= 0) {
                            d_w_to_T = cell.tn_to_node_dists[col * N + w_row];
                        }
                    }
                }

                if (d_w_to_T < SAFE_INF) {
                    int path = safe_add(w_weight, d_w_to_T);
                    if (path < min_d) min_d = path;
                }
            }
            cell.tn_to_node_dists[col * N + row] = min_d;
            if (min_d < cell.min_to_access_dist[col]) {
                cell.min_to_access_dist[col] = min_d;
            }
        }
    }
}

void ET::load_raw_graph(const string& filename) {
    ifstream in(filename);
    assert(in.is_open());
    in >> node_nums >> edge_nums;

    this->nodes.resize(node_nums + 1);
    this->graph.resize(node_nums + 1);
    node_order.resize(node_nums + 1);
    order_find.resize(node_nums + 1);

    int u, v, w;
    while (in >> u >> v >> w) {
        if (graph[u].count(v)) {
            graph[u][v] = min(graph[u][v], w);
            graph[v][u] = min(graph[v][u], w);
        } else {
            graph[u][v] = w;
            graph[v][u] = w;
        }
    }
    in.close();
    cout << "Loaded graph: " << node_nums << " nodes.\n";
}

void ET::build_tree_and_truncate_by_object() {
    cout << "--- [Phase 2] Constructing Complete Static Tree & Initializing Dynamic CH ---\n";
    initialize_static_tree_structures();
    bag_nodes.assign(node_nums + 1, vector<int>());
    bag_weights.assign(node_nums + 1, vector<int>());

    tn_infection_count.assign(node_nums + 1, 0);
    dynamic_shortcuts.assign(node_nums + 1, unordered_map<int, EdgeSupport>());
    original_upward_edges.assign(node_nums + 1, vector<pair<int, int>>());
    global_is_tn.assign(node_nums + 1, false);

    for (int rank = 1; rank <= node_nums; ++rank) {
        int u = node_order[rank];
        for (auto& edge : graph[u]) {
            int v = edge.first;
            if (order_find[v] > order_find[u]) {
                original_upward_edges[u].push_back({v, edge.second});
            }
        }
    }

    this->local_down_obj_count.assign(node_nums + 1, 0);
    for (int i = 1; i <= node_nums; ++i) {
        if (nodes[i].is_Object) this->local_down_obj_count[i] = 1;
    }

    for (int rank = 1; rank <= node_nums; ++rank) {
        int u = node_order[rank];

        unordered_map<int, int> cur_edges;
        for(auto& e : original_upward_edges[u]) cur_edges[e.first] = e.second;
        for(auto& pair : dynamic_shortcuts[u]) {
            if(!pair.second.weights.empty()) {
                int w = pair.second.min_weight;
                if(cur_edges.find(pair.first) == cur_edges.end() || w < cur_edges[pair.first]) {
                    cur_edges[pair.first] = w;
                }
            }
        }

        vector<pair<int, int>> temp_bag;
        for(auto& p : cur_edges) temp_bag.push_back({p.first, p.second});
        sort(temp_bag.begin(), temp_bag.end(), [&](pair<int,int>& a, pair<int,int>& b){
            return order_find[a.first] < order_find[b.first];
        });

        for(auto& p : temp_bag) {
            bag_nodes[u].push_back(p.first);
            bag_weights[u].push_back(p.second);
        }

        int parent = -1;
        if (!bag_nodes[u].empty()) {
            parent = bag_nodes[u][0];
            static_tree_parent[u] = parent;
            static_tree_children[parent].push_back(u);
            this->local_down_obj_count[parent] += this->local_down_obj_count[u];
        }

        // Cell granularity is controlled by alpha*K.  The old implementation
        // used MAXk here, which made the static partition independent of the
        // CELL_CAPACITY_FACTOR experiment.
        bool is_tn = false;
        const int tn_threshold = CellTargetObjectCount();
        if (local_down_obj_count[u] >= tn_threshold) is_tn = true;
        else if (tn_infection_count[u] > 0) is_tn = true;
        else if (bag_nodes[u].empty()) is_tn = true;

        global_is_tn[u] = is_tn;

        if (is_tn) {
            for(int w : bag_nodes[u]) tn_infection_count[w]++;
        }

        for (size_t i = 0; i < bag_nodes[u].size(); ++i) {
            for (size_t j = i + 1; j < bag_nodes[u].size(); ++j) {
                int v = bag_nodes[u][i], w = bag_nodes[u][j];
                int dist = safe_add(bag_weights[u][i], bag_weights[u][j]);
                int tv = order_find[v] < order_find[w] ? v : w;
                int tw = order_find[v] < order_find[w] ? w : v;
                dynamic_shortcuts[tv][tw].add(dist);
            }
        }
    }

    Assemble_Global_Cells();
}

void ET::Assemble_Global_Cells() {
    // Dynamic Cell ownership is always derived from the immutable full ET.
    // A descendant TN is a separator, so traversal stops at that node.
    tree_parent.assign(node_nums + 1, 0);
    tree_children.assign(node_nums + 1, vector<int>());
    tn_cells.assign(node_nums + 1, Cell());
    node_to_primary_tn.assign(node_nums + 1, 0);
    global_is_primary_tn.assign(node_nums + 1, false);

    for (int u = 1; u <= node_nums; ++u) {
        if (global_is_tn[u]) {
            tree_parent[u] = 0;
        } else {
            int p = static_tree_parent[u];
            tree_parent[u] = p;
            if (p != 0 && !global_is_tn[p]) {
                tree_children[p].push_back(u);
            } else if (p != 0 && global_is_tn[p]) {
                tree_children[p].push_back(u);
            }
        }
    }

    core_tns.clear();
    for (int u = 1; u <= node_nums; ++u) {
        if (!global_is_tn[u]) continue;

        core_tns.push_back(u);
        global_is_primary_tn[u] = true;
        Cell& cell = tn_cells[u];
        cell.primary_tn = u;
        cell.is_valid = true;

        // Every TN directly present in the primary bag is an access boundary.
        // Additional TNs appearing in descendants' bags are added below.
        vector<int> access;
        access.push_back(u);
        for (int w : bag_nodes[u]) {
            if (global_is_tn[w]) access.push_back(w);
        }

        vector<int> locals;
        queue<int> q;
        for (int child : static_tree_children[u]) q.push(child);

        while (!q.empty()) {
            int curr = q.front();
            q.pop();

            if (global_is_tn[curr]) {
                // Nested TN: it owns its own cell.  Do not cross the boundary.
                node_to_primary_tn[curr] = 0;
                continue;
            }

            locals.push_back(curr);
            node_to_primary_tn[curr] = u;
            for (int w : bag_nodes[curr]) {
                if (global_is_tn[w]) access.push_back(w);
            }
            for (int child : static_tree_children[curr]) q.push(child);
        }

        sort_cell_local_nodes_by_descending_rank(locals, order_find);
        cell.local_nodes = std::move(locals);

        // Deduplicate and make access order deterministic.  Column identity is
        // part of the Matrix layout, so keeping it stable simplifies rebuilds.
        std::sort(access.begin(), access.end());
        access.erase(std::unique(access.begin(), access.end()), access.end());
        cell.access_tns = std::move(access);
    }
}

unordered_set<int> ET::Dynamic_Topology_Update(const vector<int>& flip_nodes) {
    auto cmp = [&](int a, int b) { return order_find[a] > order_find[b]; };
    priority_queue<int, vector<int>, decltype(cmp)> pq(cmp);
    vector<bool> in_pq(node_nums + 1, false);

    for (int u : flip_nodes) {
        if (!in_pq[u]) { pq.push(u); in_pq[u] = true; }
    }

    unordered_set<int> dirty_cells;
    vector<int> nodes_to_trace;

    while (!pq.empty()) {
        int u = pq.top(); pq.pop();
        in_pq[u] = false;

        bool old_is_tn = global_is_tn[u];
        vector<int> old_bag = bag_nodes[u];
        vector<int> old_weights = bag_weights[u];


        unordered_map<int, int> cur_edges;
        for(auto& e : original_upward_edges[u]) cur_edges[e.first] = e.second;
        for(auto& pair : dynamic_shortcuts[u]) {
            if(!pair.second.weights.empty()) {
                int w = pair.second.min_weight;
                if(cur_edges.find(pair.first) == cur_edges.end() || w < cur_edges[pair.first]) {
                    cur_edges[pair.first] = w;
                }
            }
        }
        vector<pair<int, int>> temp_bag;
        for(auto& p : cur_edges) temp_bag.push_back({p.first, p.second});
        sort(temp_bag.begin(), temp_bag.end(), [&](pair<int,int>& a, pair<int,int>& b){
            return order_find[a.first] < order_find[b.first];
        });

        bag_nodes[u].clear(); bag_weights[u].clear();
        for(auto& p : temp_bag) {
            bag_nodes[u].push_back(p.first); bag_weights[u].push_back(p.second);
        }


        bool new_is_tn = old_is_tn;
        if (old_is_tn) {
            if (local_down_obj_count[u] < CellTargetObjectCount() && tn_infection_count[u] == 0 && !bag_nodes[u].empty()) {
                new_is_tn = false;
            }
        } else {
            if (local_down_obj_count[u] >= CellTargetObjectCount() * 2 || tn_infection_count[u] > 0 || bag_nodes[u].empty()) {
                new_is_tn = true;
            }
        }
        global_is_tn[u] = new_is_tn;

        bool bag_changed = (old_is_tn != new_is_tn) || (old_bag.size() != bag_nodes[u].size());
        if (!bag_changed) {
            for (size_t i = 0; i < old_bag.size(); ++i) {
                if (old_bag[i] != bag_nodes[u][i] || old_weights[i] != bag_weights[u][i]) {
                    bag_changed = true; break;
                }
            }
        }
        if (!bag_changed) continue;

        int old_p = old_is_tn ? 0 : (old_bag.empty() ? 0 : old_bag[0]);
        int new_p = new_is_tn ? 0 : (bag_nodes[u].empty() ? 0 : bag_nodes[u][0]);

        if (old_p != new_p) {
            auto& siblings = tree_children[old_p];
            siblings.erase(std::remove(siblings.begin(), siblings.end(), u), siblings.end());
            tree_children[new_p].push_back(u);
            tree_parent[u] = new_p;


            if (node_to_primary_tn[u] != 0) dirty_cells.insert(node_to_primary_tn[u]);
        }

        if (old_is_tn != new_is_tn) {
            if (new_is_tn) {
                core_tns.push_back(u);
                global_is_primary_tn[u] = true;
                dirty_cells.insert(u);
            } else {
                global_is_primary_tn[u] = false;
                auto it = std::find(core_tns.begin(), core_tns.end(), u);
                if (it != core_tns.end()) core_tns.erase(it);
                tn_cells[u].is_valid = false;
                tn_cells[u].local_nodes.clear();
                tn_cells[u].access_tns.clear();
            }
        }


        nodes_to_trace.push_back(u);


        if (old_is_tn) {
            for (int w : old_bag) {
                tn_infection_count[w]--;

                if (tn_infection_count[w] == 0 && local_down_obj_count[w] < CellTargetObjectCount() && !bag_nodes[w].empty()) {
                    if (!in_pq[w]) { pq.push(w); in_pq[w] = true; }
                }
            }
        } else {
            for (size_t i = 0; i < old_bag.size(); ++i) {
                for (size_t j = i + 1; j < old_bag.size(); ++j) {
                    int v = old_bag[i], w = old_bag[j];
                    int dist = safe_add(old_weights[i], old_weights[j]);
                    int tv = order_find[v] < order_find[w] ? v : w;
                    int tw = order_find[v] < order_find[w] ? w : v;
                    if (dynamic_shortcuts[tv][tw].remove(dist)) {
                        if (!in_pq[tv]) { pq.push(tv); in_pq[tv] = true; }
                    }
                }
            }
        }

        if (new_is_tn) {
            for (int w : bag_nodes[u]) {
                tn_infection_count[w]++;

                if (!global_is_tn[w]) {
                    if (!in_pq[w]) { pq.push(w); in_pq[w] = true; }
                }
            }
        } else {
            for (size_t i = 0; i < bag_nodes[u].size(); ++i) {
                for (size_t j = i + 1; j < bag_nodes[u].size(); ++j) {
                    int v = bag_nodes[u][i], w = bag_nodes[u][j];
                    int dist = safe_add(bag_weights[u][i], bag_weights[u][j]);
                    int tv = order_find[v] < order_find[w] ? v : w;
                    int tw = order_find[v] < order_find[w] ? w : v;
                    if (dynamic_shortcuts[tv][tw].add(dist)) {
                        if (!in_pq[tv]) { pq.push(tv); in_pq[tv] = true; }
                    }
                }
            }
        }
    }


    for (int u : nodes_to_trace) {
        int curr = u;
        while (curr != 0 && !global_is_primary_tn[curr]) {
            curr = tree_parent[curr];
        }
        if (curr != 0) dirty_cells.insert(curr);
    }

    for (int u : dirty_cells) {
        if (!global_is_primary_tn[u]) continue;

        Cell& cell = tn_cells[u];
        cell.primary_tn = u;
        cell.is_valid = true;

        cell.access_tns.clear();
        cell.access_tns.push_back(u);
        for (int w : bag_nodes[u]) if (w != u) cell.access_tns.push_back(w);

        vector<int> locals;
        queue<int> q;
        for (int child : tree_children[u]) q.push(child);

        while (!q.empty()) {
            int curr = q.front(); q.pop();
            locals.push_back(curr);
            node_to_primary_tn[curr] = u;
            for (int child : tree_children[curr]) q.push(child);
        }
        sort_cell_local_nodes_by_descending_rank(locals, order_find);
        cell.local_nodes = locals;
    }

    return dirty_cells;
}
