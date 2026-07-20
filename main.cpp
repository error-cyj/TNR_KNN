#include "TNR_2hop.h"
#include <iostream>
#include <fstream>
#include <vector>
#include <chrono>

using namespace std;
using namespace std::chrono;

int main() {
    // ==========================================
    // 0. Log File Redirection
    // ==========================================
    ofstream log_file("Runtime_Log.txt");
    streambuf* original_cout = cout.rdbuf();
    cout.rdbuf(log_file.rdbuf());

    // ==========================================
    // 0.1 Defensive Fix for objects.txt Format
    // ==========================================
    ifstream check_obj("objects.txt");
    if (check_obj.is_open()) {
        vector<int> temp_objs;
        int val;
        while (check_obj >> val) {
            temp_objs.push_back(val);
        }
        check_obj.close();

        if (!temp_objs.empty() && temp_objs[0] != temp_objs.size() - 1) {
            ofstream fix_obj("objects.txt");
            fix_obj << temp_objs.size() << "\n";
            for (int obj : temp_objs) {
                fix_obj << obj << "\n";
            }
            fix_obj.close();
        }
    }

    // 1. Algorithm Parameters
    int K = 20;
    TNR_2hop tnr(K);

    // 2. Preprocessing
    tnr.preprocess_KNN("", "objects.txt");

    // ==========================================
    // 3. Dynamic Threshold Calculation
    // ==========================================
    double object_density = (double)tnr.nums_object / tnr.node_nums;
    int trigger_threshold = (int)(K / object_density);
    
    // Graph Info Output (Strictly English)
    cout << "Graph Nodes: " << tnr.node_nums << "\n";
    cout << "Object Density: " << object_density << "\n";

    // 4. Read Queries
    ifstream q_file("queries.txt");
    if (!q_file.is_open()) return 1;
    vector<pair<int, int>> queries;
    int q_node, q_k;
    while (q_file >> q_node >> q_k) {
        queries.push_back({q_node, q_k});
    }
    q_file.close();

    // 5. Read Updates
    ifstream u_file("updates.txt");
    if (!u_file.is_open()) return 1;
    vector<pair<int, int>> updates;
    int old_obj, new_obj;
    while (u_file >> old_obj >> new_obj) {
        updates.push_back({old_obj, new_obj});
    }
    u_file.close();

    // 6. Initial Query Test (Strict Timing)
    long long total_initial_query_time_us = 0;
    for (const auto& q : queries) {
        vector<pair<int, int>> ans(q.second);
        auto start = steady_clock::now();
        tnr.quickQueryKNN(q.first, q.second, ans);
        auto end = steady_clock::now();
        total_initial_query_time_us += duration_cast<microseconds>(end - start).count();
    }
    if (!queries.empty()) {
        cout << "Initial Avg Query Time: " << (double)total_initial_query_time_us / queries.size() << " us\n";
    }

    // 8. Post-Update Query Test (Strict Timing)
    long long total_post_query_time_us = 0;
    for (const auto& q : queries) {
        vector<pair<int, int>> ans(q.second);
        auto start = steady_clock::now();
        tnr.quickQueryKNN(q.first, q.second, ans);
        auto end = steady_clock::now();
        total_post_query_time_us += duration_cast<microseconds>(end - start).count();
    }
    if (!queries.empty()) {
        cout << "Post-Update Avg Query Time: " << (double)total_post_query_time_us / queries.size() << " us\n";
    }

    // 7. Dynamic Updates (Strict Timing)
    long long total_update_time_us = 0;
    int update_count = 1;

    for (const auto& u : updates) {
        auto start_upd = steady_clock::now();
        tnr.Delete_Complete(u.first);
        tnr.Add_Complete(u.second);
        if (update_count % trigger_threshold == 0) {
            tnr.Apply_Index_Updates();
        }
        auto end_upd = steady_clock::now();
        total_update_time_us += duration_cast<microseconds>(end_upd - start_upd).count();
        update_count++;
    }

    if (!updates.empty()) {
        cout << "Avg Update Time: " << (double)total_update_time_us / updates.size() << " us\n";
    }

    cout.rdbuf(original_cout);

    return 0;
}
