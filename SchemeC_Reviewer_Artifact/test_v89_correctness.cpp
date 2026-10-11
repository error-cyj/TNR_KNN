#include "Fast2Hop.h"
#include <algorithm>
#include <cassert>
#include <climits>
#include <iostream>
#include <random>
#include <set>
#include <vector>

using std::cout;
using std::endl;
using std::vector;

static vector<int> dense_rank_map;
static Fast2Hop::CellOneToManyIndex test_index_holder;

static vector<int> exhaustive_topk(
    const Fast2Hop& f2h,
    int src,
    const vector<int>& targets,
    int K) {
    vector<int> best(targets.size(), INT_MAX);
    if (src <= 0 || src >= static_cast<int>(f2h.sizes.size())) return {};

    const auto& src_labels = f2h.flat_labels;
    const uint32_t soff = f2h.offsets[src];
    const uint32_t slen = f2h.sizes[src];

    for (uint32_t i = 0; i < slen; ++i) {
        const auto& sl = src_labels[soff + i];
        const uint32_t r = sl.rank;
        // Ranks are global contraction ranks and can be much larger than the
        // number of nodes in this small synthetic label graph. Do not assume
        // a small rank universe in the reference implementation.
        for (size_t slot = 0; slot < targets.size(); ++slot) {
            const int obj = targets[slot];
            if (obj <= 0 || obj >= static_cast<int>(f2h.sizes.size())) continue;
            const uint32_t off = f2h.offsets[obj];
            const uint32_t len = f2h.sizes[obj];
            for (uint32_t j = 0; j < len; ++j) {
                const auto& tl = src_labels[off + j];
                if (tl.rank == r) {
                    if (obj == src) break;
                    const long long x = static_cast<long long>(sl.dist) + tl.dist;
                    const int d = x >= INT_MAX ? INT_MAX : static_cast<int>(x);
                    if (d < best[slot]) best[slot] = d;
                    break;
                }
            }
        }
    }

    vector<int> distances;
    for (size_t i = 0; i < best.size(); ++i) {
        if (targets[i] == src) continue;
        if (best[i] < INT_MAX) distances.push_back(best[i]);
    }
    std::sort(distances.begin(), distances.end());
    if (static_cast<int>(distances.size()) > K) distances.resize(K);
    return distances;
}

static vector<int> run_v89_and_collect(
    Fast2Hop& f2h,
    int src,
    const vector<int>& targets,
    int K,
    Fast2Hop::BatchQueryStats& stats_out) {
    vector<int> out(targets.size(), INT_MAX);
    vector<int> lca(targets.size(), UINT_MAX);
    int current_D_k = INT_MAX;
    vector<int> topk;

    auto evaluator = [&](int object, int distance) {
        (void)object;
        topk.push_back(distance);
        std::sort(topk.begin(), topk.end());
        if (static_cast<int>(topk.size()) > K) topk.resize(K);
        current_D_k = (static_cast<int>(topk.size()) == K) ? topk.back() : INT_MAX;
    };

    f2h.query_1_to_N_hub_grouped(
        src, targets.data(), lca.data(), static_cast<int>(targets.size()),
        out.data(), current_D_k,
        test_index_holder,
        evaluator, &stats_out);

    std::sort(topk.begin(), topk.end());
    return topk;
}


static void build_random_case(
    Fast2Hop& f2h,
    Fast2Hop::CellOneToManyIndex& index,
    int node_count,
    int hub_count,
    int target_count,
    std::mt19937& rng) {
    f2h.dynamic_labels.assign(node_count + 1, {});
    std::uniform_int_distribution<int> dist01(1, 100);
    std::uniform_int_distribution<int> hub_dist(1, hub_count);
    std::uniform_int_distribution<int> hub_len_dist(1, std::min(hub_count, 8));

    // Source gets every hub plus a few ranks that no target uses. This stresses
    // the direct rank->group lookup and non-matching Source labels.
    vector<pair<int,int>> src_anc;
    for (int h = 1; h <= hub_count; ++h) {
        src_anc.push_back({dist01(rng), h});
    }
    for (int extra = 0; extra < hub_count / 3; ++extra) {
        const int fake_hub = hub_count + extra + 1;
        src_anc.push_back({dist01(rng), fake_hub});
    }
    f2h.update_local_node(1, src_anc, dense_rank_map);

    vector<pair<int,int>> active;
    active.reserve(target_count);
    active.push_back({0, 1}); // source itself; V8.9 must skip it.

    for (int slot = 1; slot < target_count; ++slot) {
        const int obj = slot + 1;
        active.push_back({slot, obj});
        vector<pair<int,int>> anc;
        const int len = hub_len_dist(rng);
        std::set<int> chosen;
        while (static_cast<int>(chosen.size()) < len) chosen.insert(hub_dist(rng));
        for (int h : chosen) anc.push_back({dist01(rng), h});
        // Ensure at least one target is shared by multiple Hubs; this makes the
        // object-coalescing path unavoidable in randomized cases.
        if (slot == 1 && hub_count >= 3) {
            anc.clear();
            anc.push_back({5, 1});
            anc.push_back({6, 2});
            anc.push_back({7, 3});
        }
        f2h.update_local_node(obj, anc, dense_rank_map);
    }
    // flatten all nodes, but preserve labels for node IDs used by active targets.
    f2h.flatten();
    f2h.build_cell_one_to_many_index(active, index);
}


static vector<int> pairwise_topk(
    const Fast2Hop& f2h,
    int src,
    const vector<int>& targets,
    int K) {
    vector<int> out(targets.size(), INT_MAX);
    vector<int> lca(targets.size(), 1); // All test ranks are >= 1.
    f2h.query_1_to_N_pairwise(src, targets.data(), lca.data(),
                              static_cast<int>(targets.size()), out.data(),
                              INT_MAX, nullptr);
    vector<int> distances;
    for (size_t i = 0; i < out.size(); ++i) {
        if (targets[i] == src) continue;
        if (out[i] < INT_MAX) distances.push_back(out[i]);
    }
    std::sort(distances.begin(), distances.end());
    if (static_cast<int>(distances.size()) > K) distances.resize(K);
    return distances;
}

static vector<int> ref_for_case(
    Fast2Hop& f2h,
    int src,
    const vector<int>& targets,
    int K) {
    return exhaustive_topk(f2h, src, targets, K);
}

int main() {
    // Dense rank map is only needed by update_local_node; values are identity.
    dense_rank_map.assign(64, 0);
    for (size_t i = 1; i < dense_rank_map.size(); ++i) dense_rank_map[i] = static_cast<int>(i);

    // Deterministic coalescing case.
    {
        Fast2Hop f2h;
        f2h.dynamic_labels.assign(16, {});
        f2h.update_local_node(1, {{2,1}, {3,2}, {4,3}, {100,10}}, dense_rank_map);
        f2h.update_local_node(2, {{5,1}, {50,2}, {60,3}}, dense_rank_map);
        f2h.update_local_node(3, {{6,1}, {7,2}}, dense_rank_map);
        f2h.update_local_node(4, {{8,2}, {9,3}}, dense_rank_map);
        f2h.update_local_node(5, {{100,10}}, dense_rank_map); // Hub 10 has only this target.
        f2h.flatten();

        vector<pair<int,int>> active = {{0,1},{1,2},{2,3},{3,4},{4,5}};
        Fast2Hop::CellOneToManyIndex idx;
        f2h.build_cell_one_to_many_index(active, idx);
        test_index_holder = idx;

        vector<int> targets = {1,2,3,4,5};
        const int K = 3;
        Fast2Hop::BatchQueryStats st;
        const vector<int> got = run_v89_and_collect(f2h, 1, targets, K, st);
        const vector<int> ref = ref_for_case(f2h, 1, targets, K);
        const vector<int> pair = pairwise_topk(f2h, 1, targets, K);
        assert(got == ref);
        assert(got == pair);
        assert(st.kway_heap_pops <= K);
        assert(st.kway_heap_pops <= K);
        assert(st.kway_coalesced_streams >= 1);
        assert(st.kway_distinct_objects_popped == st.kway_heap_pops);
        cout << "V8_9_DETERMINISTIC_OK\n";
        cout << "groups_created=" << st.kway_object_group_creations
             << " coalesces=" << st.kway_coalesced_streams
             << " heap_pops=" << st.kway_heap_pops
             << " updates=" << st.distance_updates << '\n';
    }

    // Regression for global ranks much larger than the toy graph size. This
    // ensures compact binary lookup does not rely on rank < |V_local|.
    {
        Fast2Hop f2h;
        f2h.dynamic_labels.assign(16, {});
        vector<int> high_rank_map(64, 0);
        for (size_t i = 1; i < high_rank_map.size(); ++i)
            high_rank_map[i] = static_cast<int>(i);
        high_rank_map[1] = 3'000'000;
        high_rank_map[2] = 7'000'000;
        high_rank_map[3] = 5'000'000;
        high_rank_map[10] = 9'000'000;

        f2h.update_local_node(1, {{2,1}, {3,2}, {4,3}, {100,10}}, high_rank_map);
        f2h.update_local_node(2, {{5,1}, {50,2}, {60,3}}, high_rank_map);
        f2h.update_local_node(3, {{6,1}, {7,2}}, high_rank_map);
        f2h.update_local_node(4, {{8,2}, {9,3}}, high_rank_map);
        f2h.update_local_node(5, {{100,10}}, high_rank_map);
        f2h.flatten();

        vector<pair<int,int>> active = {{0,1},{1,2},{2,3},{3,4},{4,5}};
        Fast2Hop::CellOneToManyIndex idx;
        f2h.build_cell_one_to_many_index(active, idx);
        assert(idx.valid);
        assert(!idx.hub_ranks.empty());
        assert(idx.hub_ranks.front() == 9'000'000U);
        for (size_t i = 1; i < idx.hub_ranks.size(); ++i)
            assert(idx.hub_ranks[i - 1] > idx.hub_ranks[i]);
        test_index_holder = idx;

        vector<int> targets = {1,2,3,4,5};
        const int K = 3;
        Fast2Hop::BatchQueryStats st;
        const vector<int> got = run_v89_and_collect(f2h, 1, targets, K, st);
        const vector<int> ref = ref_for_case(f2h, 1, targets, K);
        const vector<int> pair = pairwise_topk(f2h, 1, targets, K);
        assert(got == ref);
        assert(got == pair);
        cout << "V8_9_HIGH_GLOBAL_RANK_OK\n";
    }

    // Randomized exactness regression. The reference uses exhaustive
    // Source-label x Target-label intersections over the exact same static
    // label/index relation, so this isolates V8.9 merge correctness.
    {
        std::mt19937 rng(0x89BEEFU);
        for (int tc = 0; tc < 500; ++tc) {
            Fast2Hop f2h;
            Fast2Hop::CellOneToManyIndex idx;
            const int hub_count = 5 + (rng() % 10);
            const int target_count = 6 + (rng() % 15);
            build_random_case(f2h, idx, target_count + 2, hub_count,
                              target_count, rng);
            test_index_holder = idx;

            vector<int> targets(target_count);
            targets[0] = 1;
            for (int i = 1; i < target_count; ++i) targets[i] = i + 1;

            const int K = 1 + (rng() % std::min(5, target_count - 1));
            Fast2Hop::BatchQueryStats st;
            const vector<int> got = run_v89_and_collect(f2h, 1, targets, K, st);
            const vector<int> ref = ref_for_case(f2h, 1, targets, K);
            const vector<int> pair = pairwise_topk(f2h, 1, targets, K);

            if (got != ref || got != pair) {
                cerr << "RANDOM_MISMATCH tc=" << tc
                     << " K=" << K
                     << " got=";
                for (int x : got) cerr << x << ',';
                cerr << " ref=";
                for (int x : ref) cerr << x << ',';
                cerr << " pair=";
                for (int x : pair) cerr << x << ',';
                cerr << '\n';
                return 2;
            }
            assert(st.kway_heap_pops <= K);
            assert(st.kway_heap_pops <= K);
            assert(st.kway_distinct_objects_popped <= K);
        }
        cout << "V8_9_RANDOMIZED_500_OK\n";
    }

    return 0;
}
