#ifndef CH_TNR_FAST2HOP_H
#define CH_TNR_FAST2HOP_H

#include <vector>
#include <limits>
#include <cstdint>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <climits>
#include <chrono>

using namespace std;

class Fast2Hop {
public:
    struct QueryStats {
        long long target_queries = 0;
        long long source_label_length_sum = 0;
        long long target_label_length_sum = 0;
        long long source_entries_scanned = 0;
        long long target_entries_scanned = 0;
        long long merge_comparisons = 0;
        long long hub_matches = 0;
        long long lca_early_stops = 0;
        void reset() { *this = QueryStats(); }
    };

    struct Label {
        uint32_t rank; // order_find[hub]，严格按 Rank 降序排列
        int32_t dist;  // 到 Hub 的距离
    };

    // 动态节点 Label 数组（用于构建与增量更新阶段）
    vector<vector<Label>> dynamic_labels;

    // 平坦连续内存布局（用于极致查询性能，CSR 风格）
    vector<Label> flat_labels;
    vector<uint32_t> offsets;
    vector<uint32_t> sizes;

    const vector<int>* rank_map = nullptr;

    Fast2Hop() = default;

    // 将 dynamic_labels 打平压缩到单块连续内存 flat_labels 中
    void flatten() {
        uint32_t num_nodes = dynamic_labels.size();
        offsets.assign(num_nodes, 0);
        sizes.assign(num_nodes, 0);

        size_t total_labels = 0;
        for (uint32_t i = 1; i < num_nodes; ++i) {
            sizes[i] = static_cast<uint32_t>(dynamic_labels[i].size());
            total_labels += sizes[i];
        }

        flat_labels.resize(total_labels);

        uint32_t current_offset = 0;
        for (uint32_t i = 1; i < num_nodes; ++i) {
            offsets[i] = current_offset;
            uint32_t sz = sizes[i];
            if (sz > 0) {
                std::memcpy(flat_labels.data() + current_offset,
                            dynamic_labels[i].data(),
                            sz * sizeof(Label));
                current_offset += sz;
            }
        }
    }

    // 更新单个节点的局部祖先 Hub 集合
    void update_local_node(int u, const vector<pair<int, int>>& new_ancestors, const vector<int>& order_find) {
        if (u >= (int)dynamic_labels.size()) {
            dynamic_labels.resize(u + 1);
        }
        dynamic_labels[u].clear();
        dynamic_labels[u].reserve(new_ancestors.size());

        for (auto& anc : new_ancestors) {
            int dist = anc.first;
            int hub = anc.second;
            uint32_t r = static_cast<uint32_t>(order_find[hub]);
            dynamic_labels[u].push_back({r, dist});
        }

        // 严格按 Rank 降序排序
        std::sort(dynamic_labels[u].begin(), dynamic_labels[u].end(),
                  [](const Label& a, const Label& b) { return a.rank > b.rank; });
    }

    // 构建 Fast2Hop 索引
    void build_index(vector<vector<pair<int, int>>>& local_ancestors, const vector<int>& order_find) {
        this->rank_map = &order_find;
        int num_nodes = local_ancestors.size() - 1;
        dynamic_labels.assign(num_nodes + 1, vector<Label>());

        for (int i = 1; i <= num_nodes; ++i) {
            update_local_node(i, local_ancestors[i], order_find);
        }

        flatten();
        cout << " [Local Phase] Academic Pruned Flat Fast2Hop Built Successfully!" << endl;
    }

    // 1-to-N 批量 2hop 距离查询。
    // lca_ranks[t] 是 ET-LCA(src, targets[t]) 的 contraction rank。
    // 由于所有合法共同祖先的 rank >= rank(LCA)，因此可以在 merge 前后
    // 安全地截断低于 LCA 的 label suffix，减少无效遍历。
    struct BatchQueryStats {
        long long batch_queries = 0;
        long long target_queries = 0;
        long long source_label_length_sum = 0;
        long long target_label_length_sum = 0;
        long long source_entries_scanned = 0;
        long long target_entries_scanned = 0;
        long long matched_hub_groups = 0;
        long long hub_object_scans = 0;
        long long hub_dk_breaks = 0;
        long long source_dk_breaks = 0;
        long long object_best_prunes = 0;
        long long distance_updates = 0;
        long long pass1_time_ns = 0;
        long long pass2_time_ns = 0;
        long long source_cursor_advances = 0;
        long long hub_cursor_advances = 0;
        long long merge_comparisons = 0;
        long long physical_source_passes = 0;
        long long matched_hub_candidates = 0;
        long long lb_pruned_hub_groups = 0;
        long long processed_hub_groups = 0;
        long long target_lca_checks = 0;
        long long sort_comparisons = 0;
        long long sort_time_ns = 0;
        // V8.9 Object-Coalesced K-way counters.
        long long kway_heap_pushes = 0;              // object-group inserts
        long long kway_heap_pops = 0;                // object-group pops
        long long kway_object_group_creations = 0;
        long long kway_object_group_key_decreases = 0;
        long long kway_coalesced_streams = 0;
        long long kway_stream_entry_advances = 0;
        long long kway_finalized_object_skips = 0;
        long long kway_source_entry_skips = 0;
        long long kway_heap_peak_size = 0;
        long long kway_distinct_objects_popped = 0;
        void reset() { *this = BatchQueryStats(); }
    };

    // Query-hot-path entry. Everything is contiguous and already sorted at build time.
    struct HubObjectEntry {
        int object;
        uint32_t target_slot;
        int32_t hub_to_object;
    };

    struct CellOneToManyIndex {
        vector<uint32_t> hub_ranks;
        vector<uint32_t> hub_offsets;
        vector<HubObjectEntry> entries;

        // Compact Hub ranks are stored in strictly descending rank order.
        // The production query merges this list with the source label using
        // two pointers; no dense rank->group table or per-Hub binary search.

        bool valid = false;
        unsigned long long entry_count = 0;
    };

private:
    struct MatchedHub {
        uint32_t group;
        int32_t source_to_hub;
        uint32_t lower_bound;
        int32_t object_prev = -1;
        int32_t object_next = -1;
    };
    mutable vector<MatchedHub> matched_hubs;

    // V8.9 query workspace.  One logical heap node represents one live
    // target Object, while matched_hubs still stores the individual Hub
    // streams that currently point to that Object.  The linked lists are
    // maintained with O(1) insert/delete per stream.
    mutable vector<int32_t> object_group_head;
    mutable vector<int32_t> object_heap_pos;
    mutable vector<int32_t> object_heap_slots;
    mutable vector<int64_t> object_heap_dist;

    inline void ensure_batch_workspace(int target_count) const {
        const size_t need = static_cast<size_t>(std::max(0, target_count));
        if (matched_hubs.capacity() < need * 2 + 32)
            matched_hubs.reserve(need * 2 + 32);
        if (object_group_head.size() < need) {
            object_group_head.resize(need);
            object_heap_pos.resize(need);
        }
        object_heap_slots.reserve(need);
        object_heap_dist.reserve(need);
    }

public:
    // Build the Cell inverse index once. The preferred overload preserves the
    // Cell's active-object slot, eliminating query-time object->slot lookup.
    inline void build_cell_one_to_many_index(
        const vector<pair<int, int>>& active_objects,
        CellOneToManyIndex& out) const {
        out.hub_ranks.clear();
        out.hub_offsets.clear();
        out.entries.clear();
        out.valid = false;
        out.entry_count = 0;

        if (active_objects.empty() || flat_labels.empty()) {
            out.hub_offsets.push_back(0);
            out.valid = true;
            return;
        }

        struct TempEntry {
            uint32_t rank;
            int object;
            uint32_t slot;
            int32_t dist;
        };

        size_t reserve_sz = 0;
        for (const auto& p : active_objects) {
            const int obj = p.second;
            if (obj > 0 && obj < static_cast<int>(sizes.size()))
                reserve_sz += sizes[obj];
        }

        vector<TempEntry> tmp;
        tmp.reserve(reserve_sz);
        for (uint32_t slot = 0; slot < active_objects.size(); ++slot) {
            const int obj = active_objects[slot].second;
            if (obj <= 0 || obj >= static_cast<int>(sizes.size())) continue;
            const uint32_t off = offsets[obj], len = sizes[obj];
            const Label* p = flat_labels.data() + off;
            for (uint32_t i = 0; i < len; ++i)
                tmp.push_back({p[i].rank, obj, slot, p[i].dist});
        }

        sort(tmp.begin(), tmp.end(), [](const TempEntry& a, const TempEntry& b) {
            if (a.rank != b.rank) return a.rank > b.rank;
            if (a.dist != b.dist) return a.dist < b.dist;
            if (a.object != b.object) return a.object < b.object;
            return a.slot < b.slot;
        });

        uint32_t last_rank = 0;
        bool first = true;
        for (const auto& e : tmp) {
            if (first || e.rank != last_rank) {
                last_rank = e.rank;
                first = false;
                out.hub_ranks.push_back(e.rank);
                out.hub_offsets.push_back(static_cast<uint32_t>(out.entries.size()));
            }
            out.entries.push_back({e.object, e.slot, e.dist});
        }
        out.hub_offsets.push_back(static_cast<uint32_t>(out.entries.size()));

        // hub_ranks follows the descending-rank order used by tmp and each
        // index group is exactly its position in hub_ranks. Do not create a
        // rank_to_group array spanning the global rank universe: a Cell only
        // needs storage proportional to its actual Hub count.
        out.entry_count = out.entries.size();
        out.valid = true;
    }

    // Backward-compatible overload for the small standalone tests.
    inline void build_cell_one_to_many_index(
        const vector<int>& object_ids, CellOneToManyIndex& out) const {
        vector<pair<int, int>> active_objects;
        active_objects.reserve(object_ids.size());
        for (uint32_t i = 0; i < object_ids.size(); ++i)
            active_objects.push_back({static_cast<int>(i), object_ids[i]});
        build_cell_one_to_many_index(active_objects, out);
    }

    // V7 pairwise baseline retained for correctness and adaptive fallback.
    inline void query_1_to_N_pairwise(
        int src, const int* __restrict targets, const int* __restrict lca_ranks,
        int target_count, int* __restrict out_dists, int current_D_k,
        QueryStats* stats = nullptr) const {
        if (src >= static_cast<int>(sizes.size())) {
            for (int i=0;i<target_count;++i) out_dists[i]=INT_MAX;
            return;
        }
        const uint32_t src_len=sizes[src];
        if(src_len==0){for(int i=0;i<target_count;++i) out_dists[i]=INT_MAX; return;}
        const Label* src_ptr=flat_labels.data()+offsets[src];
        const Label* all_labels=flat_labels.data();
        const uint32_t* off=offsets.data();
        const uint32_t* sz=sizes.data();
        const int node_count=static_cast<int>(sizes.size());
        for(int t=0;t<target_count;++t){
            const int target=targets[t]; out_dists[t]=INT_MAX;
            if(target<=0||target>=node_count)continue;
            const uint32_t tl=sz[target]; if(tl==0)continue;
            const uint32_t lr=static_cast<uint32_t>(lca_ranks[t]); if(lr==0)continue;
            const Label* tp=all_labels+off[target];
            if(stats){++stats->target_queries; stats->source_label_length_sum+=src_len; stats->target_label_length_sum+=tl;}
            uint32_t i=0,j=0; int best=current_D_k; bool found=false;
            while(j<src_len && i<tl){
                const uint32_t r1=src_ptr[j].rank,r2=tp[i].rank;
                if(stats)++stats->merge_comparisons;
                if(r1<lr||r2<lr){if(stats)++stats->lca_early_stops;break;}
                if(stats){++stats->source_entries_scanned;++stats->target_entries_scanned;}
                if(r1==r2){
                    const int64_t sum=(int64_t)src_ptr[j].dist+(int64_t)tp[i].dist;
                    const int d=sum<INT_MAX?(int)sum:INT_MAX; if(stats)++stats->hub_matches;
                    if(d<best){best=d;found=true;} ++j;++i;
                }else if(r1>r2)++j;else++i;
            }
            out_dists[t]=found?best:INT_MAX;
        }
    }

    // V8.9 Scheme C production path: two-pointer merge over descending
    // source-label ranks and descending compact Hub ranks + Object-Coalesced
    // K-way merge. No per-Cell dense rank map; no per-label binary search.
    //
    // Key idea:
    //   V8.8 stores one heap entry per Hub stream. If the same Object is the
    //   current head of multiple streams, all those duplicate candidates enter
    //   the heap and are popped one by one.
    //
    //   V8.9 stores at most one heap entry per live Object.  All streams whose
    //   current head is Object o are linked in object_group_head[o].  The heap
    //   key is the minimum current-head distance among that object's streams.
    //   When o is popped, all streams pointing to o are advanced together and
    //   their next live heads are merged into their respective Object groups.
    //
    // Correctness invariant:
    //   (I1) Every active Hub stream is sorted by hub_to_object distance, hence
    //        by total d(s,h)+d(h,o).
    //   (I2) Every non-finalized Object o has exactly one heap node whose key is
    //        the minimum candidate among all currently visible stream heads for o.
    //   (I3) Every candidate hidden behind a stream head is no smaller than that
    //        head. Therefore the minimum heap key is the globally minimum
    //        candidate among all non-finalized Objects.
    //   (I4) When Object o is popped, its heap key is its exact minimum over all
    //        matched Hub occurrences. Otherwise a hidden smaller occurrence
    //        would have a preceding stream head with an even smaller candidate,
    //        contradicting the heap minimum.
    //   (I5) Once K distinct Objects are evaluated, heap.min >= D_k is safe:
    //        every unseen candidate is at least heap.min and cannot improve TOP-K.
    //
    // The heap therefore performs at most K Object-level pops. Duplicate
    // Hub/Object incidences are handled by stream coalescing/advancement rather
    // than by duplicate heap pops.
    template <class EvalFn>
    inline void query_1_to_N_hub_grouped(
        int src, const int* __restrict targets, const int* __restrict lca_ranks,
        int target_count, int* __restrict out_dists, int& current_D_k,
        const CellOneToManyIndex& index, EvalFn&& evaluate_fn,
        BatchQueryStats* stats = nullptr) const {
        if (target_count <= 0 || !index.valid || src <= 0 ||
            src >= static_cast<int>(sizes.size())) {
            for (int i = 0; i < target_count; ++i) out_dists[i] = INT_MAX;
            return;
        }

        if (stats) {
            ++stats->batch_queries;
            stats->target_queries += target_count;
        }
        for (int i = 0; i < target_count; ++i) out_dists[i] = INT_MAX;
        ensure_batch_workspace(target_count);

        // Static index path: HubObjectEntry already certifies target-label
        // membership. Keep the ABI-compatible argument but do not inspect it.
        (void)lca_ranks;

        const uint32_t src_len = sizes[src];
        const Label* sp = flat_labels.data() + offsets[src];
        if (stats) {
            stats->source_label_length_sum += src_len;
            ++stats->physical_source_passes;
        }

        const auto p1 = stats
            ? std::chrono::steady_clock::now()
            : std::chrono::steady_clock::time_point{};
        matched_hubs.clear();

        // Phase 1 (Scheme C): merge the two descending rank sequences.
        // Source labels are sorted by descending rank in update_local_node();
        // hub_ranks is emitted in the same descending order during index build.
        // A matched Hub's group id is its index in hub_ranks, so no reverse map
        // and no per-source-label binary search are required.
        uint32_t i = 0;
        uint32_t j = 0;
        const uint32_t hub_count = static_cast<uint32_t>(index.hub_ranks.size());
        while (i < src_len && j < hub_count) {
            const uint32_t sr = sp[i].rank;
            const uint32_t hr = index.hub_ranks[j];
            if (stats) {
                ++stats->merge_comparisons;
                ++stats->source_entries_scanned;
            }

            if (sr > hr) {
                // Source rank is earlier in descending order; advance source.
                ++i;
                continue;
            }
            if (sr < hr) {
                // Hub rank is earlier; this Hub cannot match current/later
                // source ranks until the Hub cursor advances.
                ++j;
                continue;
            }

            // sr == hr: source label matched Hub group j.
            const uint32_t group = j;
            const uint32_t begin = index.hub_offsets[group];
            const uint32_t end = index.hub_offsets[group + 1];
            if (begin < end) {
                matched_hubs.push_back({group, sp[i].dist, begin, -1, -1});
                if (stats) {
                    ++stats->matched_hub_candidates;
                    ++stats->matched_hub_groups;
                }
            }
            ++i;
            ++j;
        }

        if (stats) {
            stats->pass1_time_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - p1).count();
        }

        // Query-local object grouping state.
        std::fill(object_group_head.begin(),
                  object_group_head.begin() + target_count, -1);
        std::fill(object_heap_pos.begin(),
                  object_heap_pos.begin() + target_count, -1);
        object_heap_slots.clear();
        object_heap_dist.clear();

        // Indexed min-heap over OBJECT GROUPS.  The node position is stored in
        // object_heap_pos[slot], so a newly discovered better occurrence can
        // update an existing object key without inserting a duplicate heap node.
        auto heap_better = [&](int a_idx, int b_idx) {
            if (object_heap_dist[a_idx] != object_heap_dist[b_idx])
                return object_heap_dist[a_idx] < object_heap_dist[b_idx];
            return object_heap_slots[a_idx] < object_heap_slots[b_idx];
        };

        auto heap_swap_nodes = [&](int a, int b) {
            if (a == b) return;
            std::swap(object_heap_slots[a], object_heap_slots[b]);
            std::swap(object_heap_dist[a], object_heap_dist[b]);
            object_heap_pos[object_heap_slots[a]] = a;
            object_heap_pos[object_heap_slots[b]] = b;
        };

        auto heap_sift_up = [&](int pos) {
            while (pos > 0) {
                const int parent = (pos - 1) >> 1;
                if (heap_better(pos, parent)) {
                    heap_swap_nodes(pos, parent);
                    pos = parent;
                } else {
                    break;
                }
            }
        };

        auto heap_sift_down = [&](int pos) {
            const int n = static_cast<int>(object_heap_slots.size());
            while (true) {
                const int left = pos * 2 + 1;
                if (left >= n) break;
                const int right = left + 1;
                int best = left;
                if (right < n && heap_better(right, left)) best = right;
                if (heap_better(best, pos)) {
                    heap_swap_nodes(best, pos);
                    pos = best;
                } else {
                    break;
                }
            }
        };

        auto heap_insert_group = [&](uint32_t slot, int64_t distance) {
            const int pos = static_cast<int>(object_heap_slots.size());
            object_heap_slots.push_back(static_cast<int32_t>(slot));
            object_heap_dist.push_back(distance);
            object_heap_pos[slot] = pos;
            heap_sift_up(pos);
            if (stats) {
                ++stats->kway_heap_pushes;
                ++stats->kway_object_group_creations;
                stats->kway_heap_peak_size = std::max<long long>(
                    stats->kway_heap_peak_size,
                    static_cast<long long>(object_heap_slots.size()));
            }
        };

        auto heap_decrease_key = [&](uint32_t slot, int64_t distance) {
            const int pos = object_heap_pos[slot];
            if (pos < 0) return;
            if (distance >= object_heap_dist[pos]) return;
            object_heap_dist[pos] = distance;
            heap_sift_up(pos);
            if (stats) ++stats->kway_object_group_key_decreases;
        };

        auto heap_pop_group = [&]() -> int {
            const int root_slot = object_heap_slots[0];
            const int last = static_cast<int>(object_heap_slots.size()) - 1;
            object_heap_pos[root_slot] = -1;
            if (last == 0) {
                object_heap_slots.pop_back();
                object_heap_dist.pop_back();
                return root_slot;
            }
            object_heap_slots[0] = object_heap_slots[last];
            object_heap_dist[0] = object_heap_dist[last];
            object_heap_pos[object_heap_slots[0]] = 0;
            object_heap_slots.pop_back();
            object_heap_dist.pop_back();
            heap_sift_down(0);
            return root_slot;
        };

        auto link_stream = [&](uint32_t hi, uint32_t slot,
                               int64_t distance) {
            MatchedHub& mh = matched_hubs[hi];
            const int32_t old_head = object_group_head[slot];
            mh.object_prev = -1;
            mh.object_next = old_head;
            if (old_head >= 0)
                matched_hubs[static_cast<size_t>(old_head)].object_prev = static_cast<int32_t>(hi);
            object_group_head[slot] = static_cast<int32_t>(hi);

            if (old_head < 0) {
                // First visible stream for this Object.
                heap_insert_group(slot, distance);
            } else {
                // The Object is already represented once in the heap.  The new
                // stream can only keep the same key or lower it.
                if (stats) {
                    ++stats->kway_coalesced_streams;
                }
                heap_decrease_key(slot, distance);
            }
        };

        // Advance a stream from its current head.  All source/self and already
        // finalized objects are skipped before exposing the next live head.
        auto expose_next_live = [&](uint32_t hi) {
            MatchedHub& mh = matched_hubs[hi];
            const uint32_t end = index.hub_offsets[mh.group + 1];
            uint32_t p = mh.lower_bound;

            // The current head has just been consumed as part of an Object group.
            if (p < end) {
                ++p;
                if (stats) ++stats->kway_stream_entry_advances;
            }

            while (p < end) {
                const auto& e = index.entries[p];
                const uint32_t slot = e.target_slot;
                if (e.object == src) {
                    ++p;
                    if (stats) {
                        ++stats->kway_stream_entry_advances;
                        ++stats->kway_source_entry_skips;
                    }
                    continue;
                }
                if (slot >= static_cast<uint32_t>(target_count)) {
                    ++p;
                    if (stats) ++stats->kway_stream_entry_advances;
                    continue;
                }
                if (out_dists[slot] != INT_MAX) {
                    ++p;
                    if (stats) {
                        ++stats->kway_stream_entry_advances;
                        ++stats->kway_finalized_object_skips;
                    }
                    continue;
                }
                break;
            }

            mh.lower_bound = p;
            if (p >= end) return;

            const auto& e = index.entries[p];
            const int64_t raw = static_cast<int64_t>(mh.source_to_hub) +
                                static_cast<int64_t>(e.hub_to_object);
            const int64_t total = std::min<int64_t>(raw, INT_MAX);
            link_stream(hi, e.target_slot, total);
        };

        // Pass 2 includes initial stream exposure/group construction and the
        // subsequent object-level merge.
        const auto p2 = stats
            ? std::chrono::steady_clock::now()
            : std::chrono::steady_clock::time_point{};

        // Initial exposure: one head per matched Hub stream.
        for (uint32_t hi = 0; hi < matched_hubs.size(); ++hi) {
            MatchedHub& mh = matched_hubs[hi];
            const uint32_t begin = index.hub_offsets[mh.group];
            const uint32_t end = index.hub_offsets[mh.group + 1];
            uint32_t p = begin;

            while (p < end) {
                const auto& e = index.entries[p];
                const uint32_t slot = e.target_slot;
                if (e.object == src) {
                    ++p;
                    if (stats) {
                        ++stats->kway_stream_entry_advances;
                        ++stats->kway_source_entry_skips;
                    }
                    continue;
                }
                if (slot >= static_cast<uint32_t>(target_count)) {
                    ++p;
                    if (stats) ++stats->kway_stream_entry_advances;
                    continue;
                }
                if (out_dists[slot] != INT_MAX) {
                    ++p;
                    if (stats) {
                        ++stats->kway_stream_entry_advances;
                        ++stats->kway_finalized_object_skips;
                    }
                    continue;
                }
                break;
            }

            mh.lower_bound = p;
            mh.object_prev = -1;
            mh.object_next = -1;
            if (p >= end) continue;

            const auto& e = index.entries[p];
            const int64_t raw = static_cast<int64_t>(mh.source_to_hub) +
                                static_cast<int64_t>(e.hub_to_object);
            const int64_t total = std::min<int64_t>(raw, INT_MAX);
            link_stream(hi, e.target_slot, total);
        }

        while (!object_heap_slots.empty()) {
            const uint32_t slot = static_cast<uint32_t>(object_heap_slots[0]);
            const int64_t best_distance = object_heap_dist[0];
            // The heap root is the minimum among all live Object groups.  Since
            // every hidden candidate is behind an equal-or-larger stream head,
            // D_k stopping is exact.
            if (current_D_k < INT_MAX && best_distance >= current_D_k) {
                if (stats) ++stats->hub_dk_breaks;
                break;
            }

            const int popped_slot = heap_pop_group();
            (void)popped_slot;
            if (stats) {
                ++stats->kway_heap_pops;
                ++stats->kway_distinct_objects_popped;
            }

            if (slot >= static_cast<uint32_t>(target_count) ||
                object_group_head[slot] < 0) {
                // Should be unreachable for a valid index, but keep the query
                // total-defense-safe rather than dereferencing stale state.
                continue;
            }

            const int object = targets[slot];
            if (object != src && out_dists[slot] == INT_MAX) {
                const int total = best_distance < INT_MAX
                                      ? static_cast<int>(best_distance)
                                      : INT_MAX;
                out_dists[slot] = total;
                if (stats) {
                    ++stats->distance_updates;
                    ++stats->hub_object_scans;
                }
                evaluate_fn(object, total);
            }

            // Object is now finalized.  Every stream whose current head is this
            // Object can be advanced independently; all of them are processed as
            // one batch, preventing duplicate heap pops.
            int32_t hi = object_group_head[slot];
            object_group_head[slot] = -1;
            while (hi >= 0) {
                const int32_t next = matched_hubs[static_cast<size_t>(hi)].object_next;
                matched_hubs[static_cast<size_t>(hi)].object_prev = -1;
                matched_hubs[static_cast<size_t>(hi)].object_next = -1;
                if (next >= 0)
                    matched_hubs[static_cast<size_t>(next)].object_prev = -1;

                // The stream currently points to the just-finalized object.
                expose_next_live(static_cast<uint32_t>(hi));
                hi = next;
            }
        }

        if (stats) {
            stats->pass2_time_ns += std::chrono::duration_cast<std::chrono::nanoseconds>(
                std::chrono::steady_clock::now() - p2).count();
        }
    }

    // Legacy non-callback interface retained for standalone correctness tests.
    inline void query_1_to_N_inplace(
        int src, const int* __restrict targets, const int* __restrict lca_ranks,
        int target_count, int* __restrict out_dists, int current_D_k,
        const CellOneToManyIndex* index, BatchQueryStats* stats = nullptr) const {
        if (index && index->valid) {
            int mutable_dk = current_D_k;
            auto noop = [](int, int) {};
            query_1_to_N_hub_grouped(src, targets, lca_ranks, target_count, out_dists,
                                     mutable_dk, *index, noop, stats);
        } else {
            QueryStats qs;
            query_1_to_N_pairwise(src, targets, lca_ranks, target_count, out_dists, current_D_k, &qs);
            if (stats) {
                stats->batch_queries = 1;
                stats->target_queries = qs.target_queries;
                stats->source_label_length_sum = qs.source_label_length_sum;
                stats->target_label_length_sum = qs.target_label_length_sum;
                stats->source_entries_scanned = qs.source_entries_scanned;
                stats->target_entries_scanned = qs.target_entries_scanned;
                stats->matched_hub_groups = qs.hub_matches;
            }
        }
    }

};

#endif // CH_TNR_FAST2HOP_H