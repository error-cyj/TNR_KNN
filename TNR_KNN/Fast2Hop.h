#ifndef CH_TNR_FAST2HOP_H
#define CH_TNR_FAST2HOP_H

#include <vector>
#include <limits>
#include <cstdint>
#include <iostream>
#include <algorithm>

using namespace std;

class Fast2Hop {
public:

    struct Label {
        int hub;
        int dist;
    };
    vector<Label> flat_labels;
    vector<uint32_t> offsets;
    Fast2Hop() = default;
    void build_index(const vector<vector<pair<int, int>>>& local_ancestors, const vector<int>& order_find) {
        int num_nodes = local_ancestors.size() - 1;
        offsets.resize(num_nodes + 2, 0);
        uint32_t current_offset = 0;

        for (int i = 1; i <= num_nodes; ++i) {
            offsets[i] = current_offset;
            auto hubs = local_ancestors[i];
            sort(hubs.begin(), hubs.end(),
                 [&](const pair<int, int>& a, const pair<int, int>& b) {
                     return order_find[a.second] > order_find[b.second];
                 });

            for (const auto& hub : hubs) {
                flat_labels.push_back({hub.second, hub.first});
            }
            current_offset += hubs.size();
        }
        offsets[num_nodes + 1] = current_offset;
        cout << " [Local Phase] LCA-Optimized Fast2Hop Index Built!" << endl;
    }

    inline void query_1_to_N_inplace(int src, const int* targets, int target_count, int* out_dists, int current_D_k) {
        uint32_t src_start = offsets[src];
        uint32_t src_end = offsets[src + 1];
        const Label* labels = flat_labels.data();

        for (int t = 0; t < target_count; ++t) {
            int target = targets[t];
            int min_dist = current_D_k;
            bool found = false;

            uint32_t i = src_start;
            uint32_t j = offsets[target];
            uint32_t j_end = offsets[target + 1];

            while (i < src_end && j < j_end) {
                int hub_i = labels[i].hub;
                int hub_j = labels[j].hub;

                if (hub_i == hub_j) {
                    int d = labels[i].dist + labels[j].dist;
                    if (d < min_dist) {
                        min_dist = d;
                        found = true;
                    }
                    i++; j++;
                } else {
                    break;
                }
            }
            out_dists[t] = found ? min_dist : 2147483647;
        }

    }
};

#endif 