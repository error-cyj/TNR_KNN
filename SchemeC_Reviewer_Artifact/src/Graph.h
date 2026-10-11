//
// Created by chenYJ on 2024/11/5.
//

#ifndef TNR3_0_GRAPH_H
#define TNR3_0_GRAPH_H
#include <stdio.h>
#include <iostream>
#include <vector>
#include <queue>
#include <limits>
#include <unordered_map>
#include <cmath>
#include <chrono>
#include <set>
#include<sys/time.h>
#include "bits/stdc++.h"
#include <algorithm>

using namespace std::chrono;
using namespace std;

class Graph {

public:

    double get_time(void){
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return tv.tv_sec + tv.tv_usec * 1e-6;
    }

    int nums_object = 0;

    struct Node{
        bool is_Object=0 ;
    };
    static inline  double INF = numeric_limits<double>::infinity();
    vector < unordered_map <int , int > > graph;
    vector<Node> nodes;
    Graph()= default;

    void addObject(string file){
        ifstream in(file);
        assert(in.is_open());
        int nums ;
        in>>nums;
        for (int i = 0; i < nums; ++i) {
            int x;
            if (!(in >> x)) {
                std::cerr << "[Object Load Error] Truncated object file: " << file << '\n';
                break;
            }
            // Vertex IDs are 1-based. Reject 0, negatives and IDs outside the
            // loaded graph before indexing nodes; de-duplicate object records
            // so nums_object always equals the actual active-object count.
            if (x <= 0 || x >= static_cast<int>(nodes.size())) {
                std::cerr << "[Object Load Warning] Ignore invalid object vertex "
                          << x << " in " << file << '\n';
                continue;
            }
            if (!nodes[x].is_Object) {
                nodes[x].is_Object = true;
                ++nums_object;
            }
        }
    }
};

#endif //TNR3_0_GRAPH_H
