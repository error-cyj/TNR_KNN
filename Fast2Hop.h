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
    struct Label { int hub; int dist; };
    vector<vector<Label>> dynamic_labels;
    const vector<int>* rank_map = nullptr;

    Fast2Hop() = default;


    void update_local_node(int u, const vector<pair<int, int>>& new_ancestors, const vector<int>& order_find) {
        dynamic_labels[u].clear();
        for (auto& anc : new_ancestors) dynamic_labels[u].push_back({anc.second, anc.first});
        std::sort(dynamic_labels[u].begin(), dynamic_labels[u].end(),
                  [&](const Label& a, const Label& b) { return order_find[a.hub] > order_find[b.hub]; });
    }


    void build_index(vector<vector<pair<int, int>>>& local_ancestors, const vector<int>& order_find) {
        this->rank_map = &order_find;
        int num_nodes = local_ancestors.size() - 1;
        dynamic_labels.assign(num_nodes + 1, vector<Label>());

        for (int i = 1; i <= num_nodes; ++i) {
            update_local_node(i, local_ancestors[i], order_find);
        }
        cout << " [Local Phase] Academic Dynamic Fast2Hop Built!" << endl;
    }

    inline void query_1_to_N_inplace(int src, const int* targets, int target_count, int* out_dists, int current_D_k) {
        const vector<Label>& src_labels = dynamic_labels[src];
        const int* ranks = rank_map->data();
        size_t src_end = src_labels.size();

        for (int t = 0; t < target_count; ++t) {
            int target = targets[t];
            int min_dist = current_D_k;
            bool found = false;

            const vector<Label>& tgt_labels = dynamic_labels[target];
            size_t j_end = tgt_labels.size();
            size_t i = 0, j = 0;

            while (i < src_end && j < j_end) {
                int hub_i = src_labels[i].hub;
                int hub_j = tgt_labels[j].hub;

                if (hub_i == hub_j) {
                    int d = src_labels[i].dist + tgt_labels[j].dist;
                    if (d < min_dist) { min_dist = d; found = true; }
                    i++; j++;
                } else if (ranks[hub_i] > ranks[hub_j]) {
                    i++;
                } else {
                    j++;
                }

            }
            out_dists[t] = found ? min_dist : 2147483647;
        }
    }
};


#endif // CH_TNR_FAST2HOP_H