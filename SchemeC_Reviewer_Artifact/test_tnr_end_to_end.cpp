#include "TNR_2hop.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <queue>
#include <string>
#include <tuple>
#include <vector>
#include <algorithm>
#include <limits>

namespace fs = std::filesystem;
using Result = std::vector<std::pair<int, int>>; // (distance, object ID)

static Result dijkstra_topk(const TNR_2hop& tnr, int source, int k) {
    const int n = static_cast<int>(tnr.nodes.size()) - 1;
    const int inf = std::numeric_limits<int>::max() / 4;
    std::vector<int> dist(n + 1, inf);
    using Item = std::pair<int, int>;
    std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;
    dist[source] = 0;
    pq.push({0, source});
    while (!pq.empty()) {
        const auto [du, u] = pq.top();
        pq.pop();
        if (du != dist[u]) continue;
        for (const auto& [v, w] : tnr.graph[u]) {
            if (v <= 0 || v > n || w < 0) continue;
            if (du <= inf - w && du + w < dist[v]) {
                dist[v] = du + w;
                pq.push({dist[v], v});
            }
        }
    }

    Result all;
    for (int v = 1; v <= n; ++v) {
        if (tnr.nodes[v].is_Object && dist[v] < inf) all.push_back({dist[v], v});
    }
    std::sort(all.begin(), all.end(), [](const auto& a, const auto& b) {
        return a.first != b.first ? a.first < b.first : a.second < b.second;
    });
    if (static_cast<int>(all.size()) > k) all.resize(static_cast<size_t>(k));
    return all;
}

static bool verify_all_queries(TNR_2hop& tnr, const std::string& phase) {
    const int n = static_cast<int>(tnr.nodes.size()) - 1;
    const int inf = std::numeric_limits<int>::max() / 4;
    for (int source = 1; source <= n; ++source) {
        // Compute all exact source distances once; equal-distance objects at
        // the K boundary are interchangeable in a distance-only TOP-K query.
        std::vector<int> all_dist(n + 1, inf);
        using Item = std::pair<int, int>;
        std::priority_queue<Item, std::vector<Item>, std::greater<Item>> pq;
        all_dist[source] = 0;
        pq.push({0, source});
        while (!pq.empty()) {
            const auto [du, u] = pq.top();
            pq.pop();
            if (du != all_dist[u]) continue;
            for (const auto& [v, w] : tnr.graph[u]) {
                if (v <= 0 || v > n || w < 0) continue;
                if (du <= inf - w && du + w < all_dist[v]) {
                    all_dist[v] = du + w;
                    pq.push({all_dist[v], v});
                }
            }
        }
        std::vector<int> expected_distances;
        for (int v = 1; v <= n; ++v) {
            if (tnr.nodes[v].is_Object && all_dist[v] < inf) expected_distances.push_back(all_dist[v]);
        }
        std::sort(expected_distances.begin(), expected_distances.end());

        for (int k = 1; k <= 3; ++k) {
            Result got;
            tnr.quickQueryKNN(source, k, got);
            const size_t expected_count = std::min(static_cast<size_t>(k), expected_distances.size());
            std::vector<int> expected_k_distances(expected_distances.begin(),
                expected_distances.begin() + static_cast<std::ptrdiff_t>(expected_count));
            std::vector<int> got_distances;
            std::vector<int> got_ids;
            bool valid = got.size() == expected_count;
            for (const auto& [d, v] : got) {
                if (v <= 0 || v > n || !tnr.nodes[v].is_Object || all_dist[v] >= inf || all_dist[v] != d) {
                    valid = false;
                }
                got_distances.push_back(d);
                got_ids.push_back(v);
            }
            std::sort(got_distances.begin(), got_distances.end());
            std::sort(got_ids.begin(), got_ids.end());
            if (std::adjacent_find(got_ids.begin(), got_ids.end()) != got_ids.end()) valid = false;
            if (got_distances != expected_k_distances) valid = false;
            if (!valid) {
                std::cerr << "[FAIL] " << phase << " src=" << source << " K=" << k << "\n  got:";
                for (const auto& [d, v] : got) std::cerr << " (" << d << ',' << v << ')';
                std::cerr << "\n  expected distances:";
                for (int d : expected_k_distances) std::cerr << ' ' << d;
                std::cerr << '\n';
                return false;
            }
        }
    }
    std::cout << "[PASS] " << phase << " all-source K=1..3 vs Dijkstra (distance-exact, ties allowed)\n";
    return true;
}

int main() {
    const fs::path root = fs::temp_directory_path() / "v89_scheme_c_e2e_test";
    std::error_code ec;
    fs::remove_all(root, ec);
    fs::create_directories(root / "data");
    const fs::path graph_dir = root / "data";
    const fs::path graph_file = graph_dir / "Ein.txt";
    const fs::path object_file = root / "objects.txt";

    // Connected, positively weighted graph with cycles and long chords. Twenty
    // active objects exceed the 5.5*K cell target for K=3 and exercise multiple
    // Cells/TNs rather than only a single-Cell toy path.
    constexpr int graph_n = 40;
    std::vector<std::tuple<int, int, int>> edges;
    for (int u = 1; u < graph_n; ++u) edges.emplace_back(u, u + 1, 1 + (u % 5));
    edges.emplace_back(1, graph_n, 7);
    for (int u = 1; u + 3 <= graph_n; ++u) edges.emplace_back(u, u + 3, 3 + (u % 4));
    for (int u = 1; u + 7 <= graph_n; u += 2) edges.emplace_back(u, u + 7, 5 + (u % 6));
    {
        std::ofstream out(graph_file);
        out << graph_n << ' ' << edges.size() << '\n';
        for (const auto& [u, v, w] : edges) out << u << ' ' << v << ' ' << w << '\n';
    }
    {
        std::ofstream out(object_file);
        out << graph_n / 2 << '\n';
        for (int v = 2; v <= graph_n; v += 2) out << v << '\n';
    }

    TNR_2hop tnr(3);
    const std::string graph_prefix = graph_dir.string() + "/";
    tnr.preprocess_KNN(graph_prefix, object_file.string());

    if (!verify_all_queries(tnr, "static build")) return 1;

    // Lazy maintenance path: mutate object membership, call Apply exactly once,
    // then verify queries while the derived indexes remain invalidated.
    int lazy_insert = -1;
    int lazy_delete = -1;
    for (int v = 1; v < static_cast<int>(tnr.nodes.size()); ++v) {
        if (!tnr.nodes[v].is_Object && !tnr.global_is_tn[v] &&
            v < static_cast<int>(tnr.node_to_primary_tn.size()) &&
            tnr.node_to_primary_tn[v] > 0 && lazy_insert < 0) lazy_insert = v;
        if (tnr.nodes[v].is_Object && !tnr.global_is_tn[v] && lazy_delete < 0) lazy_delete = v;
    }
    if (lazy_insert < 0 || lazy_delete < 0) {
        std::cerr << "[FAIL] synthetic graph did not produce suitable non-TN insert/delete candidates\n";
        return 2;
    }
    tnr.eager_maintenance = false;
    tnr.Add_Complete(lazy_insert);
    tnr.Delete_Complete(lazy_delete);
    tnr.Apply_Index_Updates();
    if (!verify_all_queries(tnr, "lazy updates + one Apply")) return 3;

    // Eager maintenance path: make a second pair of non-TN changes. The rebuild
    // should run once at Apply and the post-rebuild answers must still be exact.
    int eager_insert = -1;
    int eager_delete = -1;
    for (int v = 1; v < static_cast<int>(tnr.nodes.size()); ++v) {
        if (!tnr.nodes[v].is_Object && !tnr.global_is_tn[v] &&
            v < static_cast<int>(tnr.node_to_primary_tn.size()) &&
            tnr.node_to_primary_tn[v] > 0 && eager_insert < 0) eager_insert = v;
        if (tnr.nodes[v].is_Object && !tnr.global_is_tn[v] && eager_delete < 0) eager_delete = v;
    }
    if (eager_insert < 0 || eager_delete < 0) {
        std::cerr << "[FAIL] synthetic graph did not produce suitable eager insert/delete candidates\n";
        return 4;
    }
    tnr.eager_maintenance = true;
    tnr.Add_Complete(eager_insert);
    tnr.Delete_Complete(eager_delete);
    tnr.Apply_Index_Updates();
    if (tnr.Has_Pending_Maintenance()) {
        std::cerr << "[FAIL] pending maintenance flag was not cleared after eager rebuild\n";
        return 5;
    }
    if (!verify_all_queries(tnr, "eager updates + rebuilt indexes")) return 6;

    fs::remove_all(root, ec);
    std::cout << "TNR_END_TO_END_DYNAMIC_OK\n";
    return 0;
}
