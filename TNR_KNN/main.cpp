#pragma comment(linker,"/STACK:1024000000,1024000000")
#include <iostream>
#include <vector>
#include <chrono>
#include <fstream>
#include <malloc.h>
#include <string>
#include <iomanip>
#include "TNR_2hop.h"
using namespace std;
using namespace std::chrono;
const string BASE_PATH = "/c++/CH-TNR-KNN/";
struct UpdateOp {
    char u;
    int v;
};
void run_single_test(string dataset, string density, string k_label) {
    cout << "------------------------------------------------------" << endl;
    cout << "Dataset: " << dataset << ", Density = " << density << ", Query K_Label = " << k_label << endl;
    cout << "Time: " << std::chrono::system_clock::to_time_t(std::chrono::system_clock::now()) << endl;
    string Graph_data = BASE_PATH + "data/" + dataset + "/Ein.txt";
    ifstream in(Graph_data);
    if (!in.is_open()) {
        cout << "Error: Cannot open graph file: " << Graph_data << endl;
        return;
    }
    int nums_nodes, nums_edges;
    in >> nums_nodes >> nums_edges;
    in.close();
    string query_data = BASE_PATH + "Queary data/" + dataset + "/_k" + k_label + "_Query.txt";
    ifstream qu(query_data);
    if (!qu.is_open()) {
        cout << "Error: Cannot open query file: " << query_data << endl;
        return;
    }
    int MAXk = 0;
    vector<pair<int, int>> query;
    int all;
    qu >> all;
    while (all--) {
        int a, k;
        qu >> a >> k;
        query.push_back({a, k});
        MAXk = max(MAXk, k);
    }
    qu.close();
    TNR_2hop tnr(MAXk , 2);
    cout << "Max K value: " << MAXk << endl;
    string obj_file = BASE_PATH + "Object data/" + dataset + "/data_" + density + ".txt";
    string index_path = BASE_PATH + "data/" + dataset + "/";
    tnr.preprocess_KNN(index_path, obj_file);
    cout << "Load Complete. Nodes: " << tnr.nodes.size() << endl;
    cout << "Starting Query Phase..." << endl;
    int query_warmup = min(1000, (int)query.size() / 10);
    vector<pair<int, int>> ans(MAXk);
    cout << "Warming up with " << query_warmup << " queries..." << endl;
    for (int i = 0; i < query_warmup; ++i) {
        tnr.quickQueryKNN(query[i].first, query[i].second, ans) ;
    }
    double query_total_time = 0, query_cnt = 0;

    for (size_t i = query_warmup; i < query.size(); ++i) {
        auto start = steady_clock::now();
        tnr.quickQueryKNN(query[i].first, query[i].second ,ans);
        auto last = duration_cast<microseconds>(steady_clock::now() - start);
        query_total_time += last.count();
        query_cnt++;
    }
    double qutime = query_cnt > 0 ? query_total_time / query_cnt : 0;
    if(query_cnt > 0){
        cout << "Avg Query Time: " << qutime << " us" << endl;
    }
    cout << "Starting Insert Phase..." << endl;
    string ins_file = BASE_PATH + "Update data/" + dataset + "/data_" + density + "_Ins.txt";
    ifstream ins_in(ins_file);
    vector<UpdateOp> ins_ops;
    if (ins_in.is_open()) {
        char u; int v; int nums;
        ins_in >> nums;
        for (int i = 0; i < nums; ++i) {
            ins_in >> u >> v;
            ins_ops.push_back({u, v});
        }
        ins_in.close();
    }
    int ins_warmup = min(1000, (int)ins_ops.size() / 10);
    cout << "Warming up with " << ins_warmup << " inserts..." << endl;
    for (int i = 0; i < ins_warmup; ++i) {
        tnr.Add_Complete(ins_ops[i].v);
    }
    double ins_total_time = 0, ins_cnt = 0;
    for (size_t i = ins_warmup; i < ins_ops.size(); ++i) {
        auto start_ins = steady_clock::now();
        tnr.Add_Complete(ins_ops[i].v);
        duration<double, std::micro> last_ins = steady_clock::now() - start_ins;
        ins_total_time += last_ins.count();
        ins_cnt++;
    }
    double instime = ins_cnt > 0 ? ins_total_time / ins_cnt : 0;
    if(ins_cnt > 0) cout << "Avg Single Insert Time: " << std::fixed << std::setprecision(6) << instime << " us" << endl;

    cout << "Starting Delete Phase..." << endl;
    string del_file = BASE_PATH + "Update data/" + dataset + "/data_" + density + "_Del.txt";
    ifstream del_in(del_file);
    vector<UpdateOp> del_ops;
    if (del_in.is_open()) {
        char u; int v; int nums;
        del_in >> nums;
        for (int i = 0; i < nums; ++i) {
            del_in >> u >> v;
            del_ops.push_back({u, v});
        }
        del_in.close();
    }

    int del_warmup = min(1000, (int)del_ops.size() / 10);
    cout << "Warming up with " << del_warmup << " deletes..." << endl;
    for (int i = 0; i < del_warmup; ++i) {
        tnr.Delete_Complete(del_ops[i].v);
    }

    long double del_total_time = 0, del_cnt = 0;
    for (size_t i = del_warmup; i < del_ops.size(); ++i) {
        auto start_del = steady_clock::now();
        tnr.Delete_Complete(del_ops[i].v);
        duration<double, std::micro> last_del = steady_clock::now() - start_del;
        del_total_time += last_del.count();
        del_cnt++;
    }
    long double deltime = del_cnt > 0 ? del_total_time / del_cnt : 0;
    if(del_cnt > 0) cout << "Avg Single Delete Time: " << std::fixed << std::setprecision(6) << deltime << " us" << endl;

    double query_throughput = qutime > 0 ? 1000000.0 / qutime : 0;
    double ins_throughput = instime > 0 ? 1000000.0 / instime : 0;
    double del_throughput = deltime > 0 ? 1000000.0 / deltime : 0;
    double total_ops = query_cnt + ins_cnt + del_cnt;
    double sumtime_total = query_total_time + ins_total_time + del_total_time;
    double total_throughput = sumtime_total > 0 ? (total_ops / (sumtime_total / 1000000.0)) : 0;
    cout << fixed << setprecision(2);
    cout << "\n--- Performance & Throughput Results ---" << endl;
    cout << "Query Throughput : " << query_throughput << " OPS" << endl;
    cout << "Insert Throughput: " << ins_throughput << " OPS" << endl;
    cout << "Delete Throughput: " << del_throughput << " OPS" << endl;
    cout << "Overall Throughput: " << total_throughput << " OPS" << endl;
    cout << "\n--- Time Distribution ---" << endl;
    if (sumtime_total > 0) {
        cout << "Query time percent : " << (query_total_time / sumtime_total) * 100 << "%" << endl;
        cout << "Insert time percent: " << (ins_total_time / sumtime_total) * 100 << "%" << endl;
        cout << "Delete time percent: " << (del_total_time / sumtime_total) * 100 << "%" << endl;
    }
    cout << endl;
}

int main() {
    vector<string> datasets = { "NY", "BAY" ,"FLA" , "COL" ,"CAL" ,"LKS" , "E" ,"W" ,"CTR" ,"USA"};
    vector<string> densities = { "0.1" ,"0.01","0.001"  ,"0.0001"};
    vector<string> k_labels = {"10", "20", "30" ,"40"};
    string log_filename = "benchmark_log.txt";
    ofstream log_file(log_filename);
    if (!log_file.is_open()) {
        cerr << "Failed to open log file!" << endl;
        return 1;
    }
    streambuf* original_cout = cout.rdbuf();
    cout.rdbuf(log_file.rdbuf());
    cout << "Benchmark Started..." << endl;

    string fixed_k = "20";
    cout << "\n========== Phase 1: Fixed K = " << fixed_k << ", Varying Densities ==========\n" << endl;
    for (const string& ds : datasets) {
        for (const string& d : densities) {
            log_file.flush();
            run_single_test(ds, d, fixed_k);

#ifdef __linux__
            malloc_trim(0);
#endif
        }
    }
    string fixed_density = "0.001";
    cout << "\n========== Phase 2: Fixed Density = " << fixed_density << ", Varying K Values ==========\n" << endl;
    for (const string& ds : datasets) {
        for (const string& k : k_labels) {
            log_file.flush();
            run_single_test(ds, fixed_density, k);
#ifdef __linux__
            malloc_trim(0);
#endif
        }
    }
    cout << "All benchmarks finished." << endl;
    cout.rdbuf(original_cout);
    cout << "Tests finished. Results saved to " << log_filename << endl;
    return 0;
}