#pragma comment(linker,"/STACK:1024000000,1024000000")

#include <algorithm>
#include <cerrno>
#include <csignal>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <unordered_set>
#include <sys/resource.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <unistd.h>

#include "src/TNR_2hop.h"

using namespace std;
using namespace std::chrono;
namespace fs = std::filesystem;

// ============================================================================
// V8.9 FINAL EXPERIMENT DRIVER
//
// Purpose:
//   Produce the final machine-readable results used for plotting/evaluation.
//   The V8.9 query algorithm itself is unchanged; this file is only the
//   experiment driver and CSV exporter.
//
// Experiment matrix (paper-oriented default):
//   Phase 1: fixed K=20, varying density {0.1, 0.01, 0.001, 0.0001}
//   Phase 2: fixed density=0.001, varying K {10, 20, 30, 40}
//   Every parameter configuration runs one combined worker:
//       preprocess -> query benchmark -> inserts -> deletes -> exactly one Apply_Index_Updates()
//   There is NO separate Phase 3. Update metrics are emitted for every density/K case.
//
// Default datasets are the same set used by the previous main program:
//   NY, BAY, FLA, COL, CAL, CTR, E, LKS, W, USA
//
// Outputs:
//   <output_dir>/benchmark_v89_final.log
//   <output_dir>/density_results.csv
//   <output_dir>/k_results.csv
//   <output_dir>/update_results.csv (one row per density/K experiment)
//   <output_dir>/all_results.csv
// ============================================================================

struct Query {
    int src = 0;
    int k = 0;
};

struct UpdateOp {
    char type = 0;
    int vertex = 0;
};

struct Options {
    string base_path = "/c++/CH-TNR-KNN/";
    string output_dir = "v89_final_results";
    string phase = "all"; // all | density | k; updates are integrated in every case
    vector<string> datasets = {
        "NY", "BAY", "FLA", "COL", "CAL",
        "CTR", "E", "LKS", "W", "USA"
    };
    vector<string> densities = {
        "0.1", "0.01", "0.001", "0.0001"
    };
    vector<string> k_labels = {"10", "20", "30", "40"};
    string fixed_k = "20";
    string fixed_density = "0.001";
    int warmup = 100;
    int update_warmup_cap = 1000;
    bool resume = false;
};

static bool parse_nonnegative_int(const string& s, int& value) {
    try {
        size_t pos = 0;
        long long x = stoll(s, &pos);
        if (pos != s.size() || x < 0 || x > 1000000000LL) return false;
        value = static_cast<int>(x);
        return true;
    } catch (...) {
        return false;
    }
}

static vector<string> split_csv(const string& text) {
    vector<string> out;
    string token;
    stringstream ss(text);
    while (getline(ss, token, ',')) {
        if (!token.empty()) out.push_back(token);
    }
    return out;
}

static string normalize_base(string path) {
    if (!path.empty() && path.back() != '/' && path.back() != '\\') {
        path.push_back('/');
    }
    return path;
}

static bool parse_options(int argc, char** argv, Options& opt) {
    auto require_value = [&](int& i, const char* flag, string& out) -> bool {
        if (i + 1 >= argc) {
            cerr << "Missing value after " << flag << '\n';
            return false;
        }
        out = argv[++i];
        return true;
    };

    for (int i = 1; i < argc; ++i) {
        const string arg = argv[i];
        string value;

        if (arg == "--help" || arg == "-h") {
            cout <<
                "Usage: " << argv[0] << " [options]\n\n"
                "  --base PATH          CH-TNR-KNN root\n"
                "  --output-dir DIR     output directory\n"
                "  --phase NAME         all|density|k (updates are integrated in each case)\n"
                "  --datasets LIST      comma-separated dataset names\n"
                "  --densities LIST     comma-separated density values\n"
                "  --k-labels LIST      comma-separated K labels\n"
                "  --fixed-k K          K used by each density experiment\n"
                "  --fixed-density D    density used by each K experiment\n"
                "  --warmup N           query warmup count (default 100)\n"
                "  --update-warmup-cap N update warmup cap (default 1000)\n"
                "  --resume              append outputs and skip completed cases\n"
                "\nDefault phase matrix:\n"
                "  density: K=20, rho={0.1,0.01,0.001,0.0001}\n"
                "  k: rho=0.001, K={10,20,30,40}\n"
                "  EVERY case: query -> insert batch -> delete batch -> exactly one Apply_Index_Updates\n";
            return false;
        } else if (arg == "--base") {
            if (!require_value(i, "--base", opt.base_path)) return false;
        } else if (arg == "--output-dir") {
            if (!require_value(i, "--output-dir", opt.output_dir)) return false;
        } else if (arg == "--phase") {
            if (!require_value(i, "--phase", opt.phase)) return false;
        } else if (arg == "--datasets") {
            if (!require_value(i, "--datasets", value)) return false;
            opt.datasets = split_csv(value);
        } else if (arg == "--densities") {
            if (!require_value(i, "--densities", value)) return false;
            opt.densities = split_csv(value);
        } else if (arg == "--k-labels") {
            if (!require_value(i, "--k-labels", value)) return false;
            opt.k_labels = split_csv(value);
        } else if (arg == "--fixed-k") {
            if (!require_value(i, "--fixed-k", opt.fixed_k)) return false;
        } else if (arg == "--fixed-density") {
            if (!require_value(i, "--fixed-density", opt.fixed_density)) return false;
        } else if (arg == "--warmup") {
            if (!require_value(i, "--warmup", value)) return false;
            if (!parse_nonnegative_int(value, opt.warmup)) {
                cerr << "Invalid --warmup: " << value << '\n';
                return false;
            }
        } else if (arg == "--resume") {
            opt.resume = true;
        } else if (arg == "--update-warmup-cap") {
            if (!require_value(i, "--update-warmup-cap", value)) return false;
            if (!parse_nonnegative_int(value, opt.update_warmup_cap)) {
                cerr << "Invalid --update-warmup-cap: " << value << '\n';
                return false;
            }
        } else {
            cerr << "Unknown option: " << arg << '\n';
            return false;
        }
    }

    opt.base_path = normalize_base(opt.base_path);
    if (opt.phase != "all" && opt.phase != "density" &&
        opt.phase != "k") {
        cerr << "Invalid --phase: " << opt.phase << '\n';
        return false;
    }
    if (opt.datasets.empty()) {
        cerr << "Dataset list is empty.\n";
        return false;
    }
    if (opt.densities.empty()) {
        cerr << "Density list is empty.\n";
        return false;
    }
    if (opt.k_labels.empty()) {
        cerr << "K list is empty.\n";
        return false;
    }
    return true;
}

static bool file_exists(const string& path) {
    ifstream in(path);
    return in.good();
}

static vector<Query> load_queries(const string& path, int& max_k) {
    ifstream in(path);
    vector<Query> queries;
    max_k = 0;

    if (!in.is_open()) {
        cerr << "ERROR: Cannot open query file: " << path << '\n';
        return queries;
    }

    int count = 0;
    if (!(in >> count) || count < 0) {
        cerr << "ERROR: Invalid query count in: " << path << '\n';
        return queries;
    }

    queries.reserve(static_cast<size_t>(count));
    for (int i = 0; i < count; ++i) {
        Query q;
        if (!(in >> q.src >> q.k)) {
            cerr << "ERROR: Invalid query record " << i << " in: " << path << '\n';
            queries.clear();
            return queries;
        }
        queries.push_back(q);
        max_k = max(max_k, q.k);
    }
    return queries;
}

static vector<UpdateOp> load_updates(const string& path, char expected_type, bool& valid) {
    ifstream in(path);
    vector<UpdateOp> ops;
    valid = true;
    if (!in.is_open()) {
        valid = false;
        return ops;
    }

    int count = 0;
    if (!(in >> count) || count < 0) {
        valid = false;
        return ops;
    }
    ops.reserve(static_cast<size_t>(count));

    for (int i = 0; i < count; ++i) {
        UpdateOp op;
        if (!(in >> op.type >> op.vertex)) {
            valid = false;
            ops.clear();
            return ops;
        }
        const char canonical = static_cast<char>(std::tolower(static_cast<unsigned char>(op.type)));
        if (canonical != expected_type) {
            valid = false;
            ops.clear();
            return ops;
        }
        ops.push_back(op);
    }
    return ops;
}

struct StructureStats {
    long long total_tns = 0;
    long long primary_cells = 0;
    long long nonempty_cells = 0;
    long long access_tn_total = 0;
    double avg_access_tns = 0.0;
};

static StructureStats collect_structure_stats(const TNR_2hop& tnr) {
    StructureStats s;
    for (int u : tnr.core_tns) {
        if (u >= 0 && u < static_cast<int>(tnr.global_is_tn.size()) &&
            tnr.global_is_tn[u]) {
            ++s.total_tns;
        }

        if (u >= 0 && u < static_cast<int>(tnr.global_is_primary_tn.size()) &&
            tnr.global_is_primary_tn[u] &&
            u < static_cast<int>(tnr.tn_cells.size()) &&
            tnr.tn_cells[u].is_valid) {
            ++s.primary_cells;
            if (!tnr.tn_cells[u].local_nodes.empty()) {
                ++s.nonempty_cells;
            }
            s.access_tn_total += static_cast<long long>(tnr.tn_cells[u].access_tns.size());
        }
    }

    if (s.primary_cells > 0) {
        s.avg_access_tns = static_cast<double>(s.access_tn_total) /
                           static_cast<double>(s.primary_cells);
    }
    return s;
}

static vector<pair<int, int>> to_query_pairs(const vector<Query>& queries) {
    vector<pair<int, int>> out;
    out.reserve(queries.size());
    for (const Query& q : queries) out.push_back({q.src, q.k});
    return out;
}

static string csv_quote(const string& s) {
    string out = "\"";
    for (char c : s) {
        if (c == '"') out += "\"\"";
        else out += c;
    }
    out += '"';
    return out;
}

class CSVWriter {
public:
    ofstream out;

    explicit CSVWriter(const string& path, bool append = false)
        : out(path, append ? (ios::out | ios::app) : (ios::out | ios::trunc)) {}

    bool good() const { return out.is_open(); }

    template <typename T>
    void field(const T& value) {
        if (!first_) out << ',';
        out << value;
        first_ = false;
    }

    void field(const string& value) {
        if (!first_) out << ',';
        out << csv_quote(value);
        first_ = false;
    }

    void field(const char* value) { field(string(value)); }

    void end_row() {
        out << '\n';
        first_ = true;
    }

private:
    bool first_ = true;
};

static void write_query_header(CSVWriter& csv) {
    const char* headers[] = {
        "experiment", "dataset", "density", "k", "query_count", "warmup",
        "measured_queries", "max_k", "preprocess_sec",
        "total_tns", "primary_cells", "nonempty_cells", "avg_access_tns",
        "access_tn_total", "tn_object_index_entries",
        "f2h_index_cells", "f2h_index_entries", "f2h_index_build_ms",
        "clean_query_us"
    };
    for (const char* h : headers) csv.field(h);
    csv.end_row();
}

static void write_update_header(CSVWriter& csv) {
    const char* headers[] = {
        "dataset", "density", "k",
        "query_count", "query_warmup", "measured_queries",
        "baseline_query_us",
        "insert_count", "insert_warmup", "measured_inserts", "avg_insert_us",
        "delete_count", "delete_warmup", "measured_deletes", "avg_delete_us",
        "apply_calls", "apply_time_us",
        "avg_update_time_us",
        "index_maintenance_pending_at_end"
    };
    for (const char* h : headers) csv.field(h);
    csv.end_row();
}

static void write_query_header_if_needed(CSVWriter& csv, bool write_header) {
    if (write_header) write_query_header(csv);
}

static void write_update_header_if_needed(CSVWriter& csv, bool write_header) {
    if (write_header) write_update_header(csv);
}

static void write_query_result_exact(CSVWriter& csv,
                                     const string& experiment,
                                     const string& dataset,
                                     const string& density,
                                     const string& k_label,
                                     long long query_count,
                                     int warmup,
                                     int max_k,
                                     const TNR_2hop& tnr,
                                     const StructureStats& st,
                                     double preprocess_sec) {
    csv.field(experiment);
    csv.field(dataset);
    csv.field(density);
    csv.field(k_label);
    csv.field(query_count);
    csv.field(warmup);
    csv.field(tnr.last_clean_query_count);
    csv.field(max_k);
    csv.field(preprocess_sec);
    csv.field(st.total_tns);
    csv.field(st.primary_cells);
    csv.field(st.nonempty_cells);
    csv.field(st.avg_access_tns);
    csv.field(st.access_tn_total);
    csv.field(tnr.tn_object_index_entries);
    csv.field(tnr.local_f2h_index_cells_total);
    csv.field(tnr.local_f2h_index_entries_total);
    csv.field(tnr.local_f2h_index_build_time_ms);
    csv.field(tnr.last_clean_avg_query_us);
    csv.end_row();
}

static bool valid_update_vertex(const TNR_2hop& tnr, int vertex) {
    return vertex > 0 && vertex < static_cast<int>(tnr.nodes.size());
}

// Run update workload on the SAME freshly-built index after the query benchmark.
// Exact requested sequence per parameter configuration:
//   1) Query benchmark on clean index
//   2) All insert operations (insert warmup is not included in timing average)
//   3) All delete operations (delete warmup is not included in timing average)
//   4) Exactly one explicit Apply_Index_Updates() call at the end
//
// Average update time follows the agreed measurement convention:
//   avg_insert_us + avg_delete_us + single_apply_call_us.
static void run_update_batch_after_query(const Options& opt,
                                         const string& dataset,
                                         const string& density,
                                         const string& k_label,
                                         long long query_count,
                                         int query_warmup,
                                         double baseline_query_us,
                                         TNR_2hop& tnr,
                                         const vector<UpdateOp>& ins_ops,
                                         const vector<UpdateOp>& del_ops,
                                         CSVWriter& update_csv,
                                         ofstream& log) {
    const int insert_warmup = min(opt.update_warmup_cap,
                                  static_cast<int>(ins_ops.size()) / 10);
    const size_t insert_begin = static_cast<size_t>(insert_warmup);

    log << "\nStarting Insert Phase (after baseline query)...\n";
    log << "Insert operations available = " << ins_ops.size() << '\n';
    log << "Insert warmup               = " << insert_warmup << '\n';

    for (int i = 0; i < insert_warmup &&
                     i < static_cast<int>(ins_ops.size()); ++i) {
        const auto& op = ins_ops[static_cast<size_t>(i)];
        if (!valid_update_vertex(tnr, op.vertex)) continue;
        if (op.type == 'i' || op.type == 'I') {
            tnr.Add_Complete(op.vertex);
        }
    }

    double insert_total_us = 0.0;
    long long measured_inserts = 0;
    for (size_t i = insert_begin; i < ins_ops.size(); ++i) {
        const auto& op = ins_ops[i];
        if (!valid_update_vertex(tnr, op.vertex)) continue;
        if (!(op.type == 'i' || op.type == 'I')) continue;

        const auto start = steady_clock::now();
        tnr.Add_Complete(op.vertex);
        insert_total_us += duration<double, micro>(steady_clock::now() - start).count();
        ++measured_inserts;
    }
    const double avg_insert_us = measured_inserts > 0
        ? insert_total_us / static_cast<double>(measured_inserts) : 0.0;

    log << "Measured Inserts = " << measured_inserts << '\n';
    log << "Avg Insert Time  = " << avg_insert_us << " us\n";

    const int delete_warmup = min(opt.update_warmup_cap,
                                  static_cast<int>(del_ops.size()) / 10);
    const size_t delete_begin = static_cast<size_t>(delete_warmup);

    log << "\nStarting Delete Phase (still no index maintenance)...\n";
    log << "Delete operations available = " << del_ops.size() << '\n';
    log << "Delete warmup               = " << delete_warmup << '\n';

    for (int i = 0; i < delete_warmup &&
                     i < static_cast<int>(del_ops.size()); ++i) {
        const auto& op = del_ops[static_cast<size_t>(i)];
        if (!valid_update_vertex(tnr, op.vertex)) continue;
        if (op.type == 'd' || op.type == 'D') {
            tnr.Delete_Complete(op.vertex);
        }
    }

    double delete_total_us = 0.0;
    long long measured_deletes = 0;
    for (size_t i = delete_begin; i < del_ops.size(); ++i) {
        const auto& op = del_ops[i];
        if (!valid_update_vertex(tnr, op.vertex)) continue;
        if (!(op.type == 'd' || op.type == 'D')) continue;

        const auto start = steady_clock::now();
        tnr.Delete_Complete(op.vertex);
        delete_total_us += duration<double, micro>(steady_clock::now() - start).count();
        ++measured_deletes;
    }
    const double avg_delete_us = measured_deletes > 0
        ? delete_total_us / static_cast<double>(measured_deletes) : 0.0;

    log << "Measured Deletes = " << measured_deletes << '\n';
    log << "Avg Delete Time  = " << avg_delete_us << " us\n";

    // Exactly one explicit maintenance call, after all inserts and deletes.
    log << "\nCalling Apply_Index_Updates() exactly once after the update batch...\n";
    const auto apply_start = steady_clock::now();
    tnr.Apply_Index_Updates();
    const double apply_us = duration<double, micro>(steady_clock::now() - apply_start).count();

    // Preserve the user's agreed reporting formula. Apply is invoked once, so
    // apply_us is the measured duration of that one final maintenance call.
    const double avg_update_us = avg_insert_us + avg_delete_us + apply_us;
    const bool pending = tnr.Has_Pending_Maintenance();

    log << "Apply_Index_Updates Calls = 1\n";
    log << "Apply_Index_Updates Time  = " << apply_us << " us\n";
    log << "Average Update Time       = avg_insert + avg_delete + apply_call\n";
    log << "                          = " << avg_insert_us << " + "
        << avg_delete_us << " + " << apply_us << " = " << avg_update_us << " us\n";
    log << "Index Maintenance Pending = " << (pending ? "YES" : "NO") << '\n';

    update_csv.field(dataset);
    update_csv.field(density);
    update_csv.field(k_label);
    update_csv.field(query_count);
    update_csv.field(query_warmup);
    update_csv.field(tnr.last_clean_query_count);
    update_csv.field(baseline_query_us);
    update_csv.field(static_cast<long long>(ins_ops.size()));
    update_csv.field(insert_warmup);
    update_csv.field(measured_inserts);
    update_csv.field(avg_insert_us);
    update_csv.field(static_cast<long long>(del_ops.size()));
    update_csv.field(delete_warmup);
    update_csv.field(measured_deletes);
    update_csv.field(avg_delete_us);
    update_csv.field(1LL); // exactly one Apply_Index_Updates() call
    update_csv.field(apply_us);
    update_csv.field(avg_update_us);
    update_csv.field(pending ? 1 : 0);
    update_csv.end_row();
    update_csv.out.flush();
    log.flush();
}

static bool run_query_experiment(const Options& opt,
                                 const string& experiment,
                                 const string& dataset,
                                 const string& density,
                                 const string& k_label,
                                 CSVWriter& query_csv,
                                 CSVWriter& update_csv,
                                 ofstream& log) {
    const string graph_path = opt.base_path + "data/" + dataset + "/Ein.txt";
    const string query_path = opt.base_path + "Queary data/" + dataset +
                              "/_k" + k_label + "_Query.txt";
    const string object_path = opt.base_path + "Object data/" + dataset +
                               "/data_" + density + ".txt";
    const string ins_path = opt.base_path + "Update data/" + dataset +
                            "/data_" + density + "_Ins.txt";
    const string del_path = opt.base_path + "Update data/" + dataset +
                            "/data_" + density + "_Del.txt";
    const string index_path = opt.base_path + "data/" + dataset + "/";

    if (!file_exists(graph_path) || !file_exists(query_path) || !file_exists(object_path)) {
        log << "SKIP: missing graph/query/object input for dataset=" << dataset
            << " density=" << density << " k=" << k_label << '\n';
        return false;
    }

    int max_k = 0;
    const vector<Query> queries = load_queries(query_path, max_k);
    if (queries.empty()) {
        log << "SKIP: empty/invalid query set: " << query_path << '\n';
        return false;
    }

    // Load both update batches before preprocessing. Missing update files do not
    // invalidate query timing, but are clearly recorded and produce zero counts.
    bool ins_file_valid = true;
    bool del_file_valid = true;
    const bool ins_file_exists = file_exists(ins_path);
    const bool del_file_exists = file_exists(del_path);
    const vector<UpdateOp> ins_ops = ins_file_exists
        ? load_updates(ins_path, 'i', ins_file_valid) : vector<UpdateOp>{};
    const vector<UpdateOp> del_ops = del_file_exists
        ? load_updates(del_path, 'd', del_file_valid) : vector<UpdateOp>{};
    if (!ins_file_exists) log << "WARNING: insert file missing: " << ins_path << '\n';
    if (!del_file_exists) log << "WARNING: delete file missing: " << del_path << '\n';
    if (!ins_file_valid || !del_file_valid) {
        log << "SKIP: malformed or wrong-type update file (insert file must contain i/I records; "
               "delete file must contain d/D records).\n";
        return false;
    }
    if (ins_ops.empty() && del_ops.empty()) {
        log << "WARNING: both insert and delete operation sets are empty; "
               "Apply_Index_Updates() will still be called exactly once.\n";
    }

    int graph_nodes = 0;
    int graph_edges = 0;
    {
        ifstream in(graph_path);
        if (!(in >> graph_nodes >> graph_edges)) {
            log << "SKIP: invalid graph metadata: " << graph_path << '\n';
            return false;
        }
    }

    log << "\n------------------------------------------------------------\n";
    log << "Experiment = " << experiment << " (query + update metrics)\n";
    log << "Dataset    = " << dataset << '\n';
    log << "Density    = " << density << '\n';
    log << "K label    = " << k_label << '\n';
    log << "Graph      = " << graph_path << '\n';
    log << "Objects    = " << object_path << '\n';
    log << "Queries    = " << query_path << '\n';
    log << "Insert ops = " << ins_path << '\n';
    log << "Delete ops = " << del_path << '\n';
    // Fail fast on invalid query vertices/K values. The query implementation
    // uses 1-based graph IDs and a fixed LIMIT_K-sized Top-K heap.
    for (size_t qi = 0; qi < queries.size(); ++qi) {
        const Query& q = queries[qi];
        if (q.src <= 0 || q.src > graph_nodes || q.k <= 0 || q.k > 50) {
            log << "SKIP: invalid query record at index " << qi
                << " (src=" << q.src << ", K=" << q.k
                << "); expected 1 <= src <= " << graph_nodes
                << " and 1 <= K <= 50.\n";
            return false;
        }
    }
    auto validate_update_vertices = [&](const vector<UpdateOp>& ops, const char* label) {
        for (size_t i = 0; i < ops.size(); ++i) {
            if (ops[i].vertex <= 0 || ops[i].vertex > graph_nodes) {
                log << "SKIP: invalid " << label << " vertex at record " << i
                    << " (vertex=" << ops[i].vertex << "); expected 1.."
                    << graph_nodes << ".\n";
                return false;
            }
        }
        return true;
    };
    if (!validate_update_vertices(ins_ops, "insert") ||
        !validate_update_vertices(del_ops, "delete")) return false;

    log << "Graph      = " << graph_nodes << " nodes, " << graph_edges << " edges\n";
    log << "Query count= " << queries.size() << ", max K=" << max_k << '\n';
    log << "Sequence   = Query -> Insert batch -> Delete batch -> exactly one Apply_Index_Updates()\n";
    log.flush();

    TNR_2hop tnr(max_k);
    log << "Worker stage = preprocess_begin\n";
    log.flush();
    const auto t0 = steady_clock::now();
    tnr.preprocess_KNN(index_path, object_path);
    const double preprocess_sec = duration<double>(steady_clock::now() - t0).count();

    log << "Preprocess = " << fixed << setprecision(6) << preprocess_sec << " s\n";
    log << "Nodes      = " << tnr.nodes.size() << '\n';
    log.flush();

    const StructureStats structure = collect_structure_stats(tnr);
    tnr.print_statistics();

    // (1) Clean-only query benchmark before any update. The warmup and measured
    // queries come only from the real query file; there is no internal/profiling pass.
    log << "Starting clean-only query benchmark before inserts/deletes...\n";
    log.flush();
    const vector<pair<int, int>> query_pairs = to_query_pairs(queries);
    const double baseline_query_us = tnr.run_clean_query_benchmark(query_pairs, opt.warmup);
    int requested_warmup = opt.warmup;
    if (requested_warmup < 0) {
        requested_warmup = min(1000, static_cast<int>(queries.size()) / 10);
    }
    const int effective_query_warmup = max(0, min(
        requested_warmup, max(0, static_cast<int>(queries.size()) - 1)));
    const long long measured_query_count = tnr.last_clean_query_count;

    write_query_result_exact(query_csv, experiment, dataset, density, k_label,
                             static_cast<long long>(queries.size()),
                             effective_query_warmup, max_k, tnr, structure, preprocess_sec);
    query_csv.out.flush();
    log << "Warmup queries              = " << effective_query_warmup << '\n';
    log << "Measured clean queries      = " << measured_query_count << '\n';
    log << "Clean Avg Query Time        = " << baseline_query_us << " us\n";
    log << "Query profiling             = DISABLED\n";
    log.flush();

    // (2)-(4) On this same fresh index: insert, delete, then one final Apply.
    run_update_batch_after_query(opt, dataset, density, k_label,
                                 static_cast<long long>(queries.size()),
                                 effective_query_warmup, baseline_query_us,
                                 tnr, ins_ops, del_ops, update_csv, log);

    log.flush();
    return true;
}


// ---------------------------------------------------------------------------
// Process-isolated experiment runner
//
// Each (phase, dataset, density, K) executes in a separate child process.
// This is intentionally still single-threaded. The parent is a small scheduler
// and never constructs a road-network index, so memory from one trial is fully
// returned to the OS when that child exits. wait4() reports peak RSS even for
// a child terminated by SIGKILL on Linux.
// ---------------------------------------------------------------------------
struct MemorySnapshot {
    double host_ram_gb = -1.0;
    double cgroup_limit_gb = -1.0;
    double cgroup_current_gb = -1.0;
    long long oom_kill_count = -1;
};

static double bytes_to_gb(long long value) {
    return value >= 0 ? static_cast<double>(value) / (1024.0 * 1024.0 * 1024.0) : -1.0;
}

static long long read_numeric_file(const string& path) {
    ifstream in(path);
    string value;
    if (!(in >> value) || value == "max") return -1;
    try {
        size_t pos = 0;
        long long n = stoll(value, &pos);
        if (pos != value.size()) return -1;
        return n;
    } catch (...) {
        return -1;
    }
}

static string detect_memory_cgroup_dir() {
    // Prefer the process's actual cgroup path rather than assuming it lives at
    // the cgroup mount root. This matters under Slurm/systemd/container limits.
    ifstream in("/proc/self/cgroup");
    string line;
    while (getline(in, line)) {
        const size_t p = line.find("0::");
        if (p == 0) {
            string rel = line.substr(3);
            while (!rel.empty() && rel.front() == '/') rel.erase(rel.begin());
            const string v2 = string("/sys/fs/cgroup/") + rel;
            if (fs::exists(v2 + "/memory.max")) return v2;
        }
        const size_t first = line.find(':');
        if (first == string::npos) continue;
        const size_t second = line.find(':', first + 1);
        if (second == string::npos) continue;
        const string controllers = line.substr(first + 1, second - first - 1);
        const string rel = line.substr(second + 1);
        if (controllers.find("memory") != string::npos) {
            string clean = rel;
            while (!clean.empty() && clean.front() == '/') clean.erase(clean.begin());
            const string v1 = string("/sys/fs/cgroup/memory/") + clean;
            if (fs::exists(v1 + "/memory.limit_in_bytes")) return v1;
        }
    }
    return "/sys/fs/cgroup";
}

static const string& memory_cgroup_dir() {
    static const string path = detect_memory_cgroup_dir();
    return path;
}

static long long read_host_ram_kb() {
    ifstream in("/proc/meminfo");
    string key, unit;
    long long value = -1;
    while (in >> key >> value >> unit) {
        if (key == "MemTotal:") return value;
    }
    return -1;
}

static long long read_cgroup_oom_kill_count() {
    const string dir = memory_cgroup_dir();
    {
        ifstream in(dir + "/memory.events");
        string key;
        long long value = -1;
        while (in >> key >> value) {
            if (key == "oom_kill") return value;
        }
    }
    // cgroup v1 fallback. Some kernels expose OOM counters in memory.oom_control.
    {
        ifstream in(dir + "/memory.oom_control");
        string key;
        long long value = -1;
        while (in >> key >> value) {
            if (key == "oom_kill") return value;
        }
    }
    return -1;
}

static MemorySnapshot sample_memory() {
    MemorySnapshot m;
    const long long host_kb = read_host_ram_kb();
    if (host_kb >= 0) m.host_ram_gb = static_cast<double>(host_kb) / (1024.0 * 1024.0);
    const string dir = memory_cgroup_dir();
    long long limit = read_numeric_file(dir + "/memory.max");
    long long current = read_numeric_file(dir + "/memory.current");
    if (limit < 0) limit = read_numeric_file(dir + "/memory.limit_in_bytes");
    if (current < 0) current = read_numeric_file(dir + "/memory.usage_in_bytes");
    m.cgroup_limit_gb = bytes_to_gb(limit);
    m.cgroup_current_gb = bytes_to_gb(current);
    m.oom_kill_count = read_cgroup_oom_kill_count();
    return m;
}

static string safe_name(string s) {
    for (char& c : s) {
        if (!(isalnum(static_cast<unsigned char>(c)) || c == '-' || c == '_' || c == '.')) c = '_';
    }
    return s;
}

static string make_case_id(const string& experiment, const string& dataset,
                           const string& density, const string& k_label) {
    return experiment + "|" + dataset + "|" + density + "|" + k_label;
}

static unordered_set<string> load_completed_cases(const string& status_path) {
    unordered_set<string> completed;
    ifstream in(status_path);
    string line;
    bool first = true;
    while (getline(in, line)) {
        if (first) { first = false; continue; }
        if (line.empty()) continue;
        stringstream ss(line);
        string case_id, field, state;
        if (!getline(ss, case_id, ',')) continue;
        // state is column 12 (0-based index 11); discard columns 2..11.
        for (int col = 1; col < 11; ++col) {
            if (!getline(ss, field, ',')) { field.clear(); break; }
        }
        if (!getline(ss, state, ',')) continue;
        if (state == "OK") completed.insert(case_id);
    }
    return completed;
}

static void append_csv_row(const string& worker_csv_path, CSVWriter& destination) {
    ifstream in(worker_csv_path);
    string line;
    if (!getline(in, line)) return; // worker header
    while (getline(in, line)) {
        if (line.empty()) continue;
        destination.out << line << '\n';
        destination.out.flush();
        return; // a worker should emit exactly one row
    }
}

static void append_worker_log(const string& worker_log_path, ofstream& parent_log) {
    ifstream in(worker_log_path);
    if (!in.is_open()) return;
    parent_log << in.rdbuf();
    if (parent_log.tellp() > 0) parent_log << '\n';
    parent_log.flush();
}

static void write_status_header_if_needed(ofstream& out, bool write_header) {
    if (!write_header) return;
    out << "case_id,experiment,dataset,density,k,wall_sec,peak_rss_mb,host_ram_gb,"
           "cgroup_limit_gb,cgroup_current_before_gb,oom_kill_delta,state,exit_code,signal\n";
    out.flush();
}

static void write_status_row(ofstream& out, const string& case_id,
                             const string& experiment, const string& dataset,
                             const string& density, const string& k_label,
                             double wall_sec, double peak_rss_mb,
                             const MemorySnapshot& before, long long oom_delta,
                             const string& state, int exit_code, int signal_number) {
    out << case_id << ',' << experiment << ',' << dataset << ',' << density << ',' << k_label << ','
        << fixed << setprecision(6) << wall_sec << ',' << peak_rss_mb << ','
        << before.host_ram_gb << ',' << before.cgroup_limit_gb << ','
        << before.cgroup_current_gb << ',' << oom_delta << ',' << state << ','
        << exit_code << ',' << signal_number << '\n';
    out.flush();
}

struct IsolatedRunResult {
    string state = "FAILED";
    int exit_code = -1;
    int signal_number = 0;
    double wall_sec = 0.0;
    double peak_rss_mb = -1.0;
};

static IsolatedRunResult run_isolated_case(
        const Options& opt, const string& experiment, const string& dataset,
        const string& density, const string& k_label,
        CSVWriter& result_csv, CSVWriter* all_query_csv, CSVWriter& update_csv,
        ofstream& parent_log, ofstream& status_out,
        unordered_set<string>& completed_cases, unsigned long long& case_serial) {
    // New schema/version marker: old query-only cases must not be treated as
    // complete because they have no corresponding update_results.csv row.
    const string case_id = make_case_id("query_update_" + experiment, dataset, density, k_label);
    if (opt.resume && completed_cases.count(case_id)) {
        parent_log << "RESUME SKIP: already completed " << case_id << '\n';
        parent_log.flush();
        IsolatedRunResult skipped;
        skipped.state = "RESUME_SKIPPED";
        return skipped;
    }

    const string stem = opt.output_dir + "/.worker_" + to_string(++case_serial) + "_" +
                        safe_name(experiment + "_" + dataset + "_" + density + "_" + k_label);
    const string worker_log_path = stem + ".log";
    const string worker_csv_path = stem + ".csv";
    const string worker_update_csv_path = stem + ".update.csv";

    cout.flush();
    parent_log.flush();
    result_csv.out.flush();
    if (all_query_csv) all_query_csv->out.flush();
    update_csv.out.flush();
    status_out.flush();

    const MemorySnapshot memory_before = sample_memory();
    const auto wall_start = steady_clock::now();
    const pid_t pid = fork();
    if (pid < 0) {
        parent_log << "ERROR: fork failed for " << case_id << ", errno=" << errno << '\n';
        parent_log.flush();
        write_status_row(status_out, case_id, experiment, dataset, density, k_label,
                         0.0, -1.0, memory_before, -1, "FORK_FAILED", -1, 0);
        IsolatedRunResult failed;
        failed.state = "FORK_FAILED";
        return failed;
    }

    if (pid == 0) {
        // Child: open private outputs after fork, never write directly to the
        // parent's result CSVs. _exit() at the end ensures all process memory
        // is returned even if allocators retain large arenas.
        ofstream worker_log(worker_log_path, ios::out | ios::trunc);
        CSVWriter worker_csv(worker_csv_path);
        CSVWriter worker_update_csv(worker_update_csv_path);
        if (!worker_log.is_open() || !worker_csv.good() || !worker_update_csv.good()) _exit(101);

        streambuf* original_child_cout = cout.rdbuf(worker_log.rdbuf());
        bool completed = false;
        int child_exit = 0;
        try {
            worker_log << "Worker case = " << case_id << '\n';
            worker_log << "Worker PID  = " << static_cast<long long>(getpid()) << '\n';
            worker_log << "Worker isolation = one configuration per process\n";
            worker_log.flush();
            write_query_header(worker_csv);
            write_update_header(worker_update_csv);
            completed = run_query_experiment(opt, experiment, dataset, density,
                                             k_label, worker_csv, worker_update_csv, worker_log);
        } catch (const bad_alloc&) {
            worker_log << "WORKER_EXCEPTION = std::bad_alloc\n";
            worker_log.flush();
            child_exit = 4;
        } catch (const std::exception& e) {
            worker_log << "WORKER_EXCEPTION = " << e.what() << '\n';
            worker_log.flush();
            child_exit = 5;
        } catch (...) {
            worker_log << "WORKER_EXCEPTION = unknown\n";
            worker_log.flush();
            child_exit = 6;
        }

        cout.flush();
        cout.rdbuf(original_child_cout);
        worker_csv.out.flush();
        worker_csv.out.close();
        worker_update_csv.out.flush();
        worker_update_csv.out.close();
        worker_log.flush();
        worker_log.close();
        if (child_exit != 0) _exit(child_exit);
        _exit(completed ? 0 : 3); // 3 means SKIPPED / incomplete input, not a crash
    }

    int wait_status = 0;
    struct rusage usage {};
    pid_t waited;
    do {
        waited = wait4(pid, &wait_status, 0, &usage);
    } while (waited < 0 && errno == EINTR);

    IsolatedRunResult result;
    result.wall_sec = duration<double>(steady_clock::now() - wall_start).count();
    if (waited < 0) {
        result.state = "WAIT_FAILED";
    } else {
#if defined(__APPLE__)
        result.peak_rss_mb = static_cast<double>(usage.ru_maxrss) / (1024.0 * 1024.0);
#else
        // Linux reports ru_maxrss in KiB.
        result.peak_rss_mb = static_cast<double>(usage.ru_maxrss) / 1024.0;
#endif
        if (WIFEXITED(wait_status)) {
            result.exit_code = WEXITSTATUS(wait_status);
            if (result.exit_code == 0) result.state = "OK";
            else if (result.exit_code == 3) result.state = "SKIPPED";
            else result.state = "FAILED_EXIT";
        } else if (WIFSIGNALED(wait_status)) {
            result.signal_number = WTERMSIG(wait_status);
            result.state = "KILLED_SIGNAL";
        }
    }

    const MemorySnapshot memory_after = sample_memory();
    const long long oom_delta = (memory_before.oom_kill_count >= 0 && memory_after.oom_kill_count >= 0)
                                    ? memory_after.oom_kill_count - memory_before.oom_kill_count
                                    : -1;

    parent_log << "\n--- CHILD RESULT: " << case_id << " ---\n";
    append_worker_log(worker_log_path, parent_log);
    parent_log << "Worker status       = " << result.state << '\n'
               << "Worker exit code    = " << result.exit_code << '\n'
               << "Worker signal       = " << result.signal_number << '\n'
               << "Worker wall time    = " << fixed << setprecision(6) << result.wall_sec << " s\n"
               << "Worker peak RSS     = " << result.peak_rss_mb << " MiB\n"
               << "Host RAM            = " << memory_before.host_ram_gb << " GiB\n"
               << "Cgroup memory limit = " << memory_before.cgroup_limit_gb << " GiB (-1 means unavailable/unlimited)\n"
               << "Cgroup current pre  = " << memory_before.cgroup_current_gb << " GiB\n"
               << "Cgroup OOM kill delta = " << oom_delta << '\n';
    if (result.state == "KILLED_SIGNAL" && result.signal_number == SIGKILL && oom_delta > 0) {
        parent_log << "DIAGNOSIS = SIGKILL and cgroup oom_kill counter increased; likely cgroup OOM kill.\n";
    } else if (result.state == "KILLED_SIGNAL" && result.signal_number == SIGKILL) {
        parent_log << "DIAGNOSIS = SIGKILL observed; cgroup OOM counter did not confirm OOM (check kernel log / scheduler limits).\n";
    }
    parent_log << "--- END CHILD RESULT ---\n\n";
    parent_log.flush();

    if (result.state == "OK") {
        append_csv_row(worker_csv_path, result_csv);
        if (all_query_csv) append_csv_row(worker_csv_path, *all_query_csv);
        append_csv_row(worker_update_csv_path, update_csv);
        completed_cases.insert(case_id);
    }

    write_status_row(status_out, case_id, experiment, dataset, density, k_label,
                     result.wall_sec, result.peak_rss_mb, memory_before, oom_delta,
                     result.state, result.exit_code, result.signal_number);

    std::error_code ec;
    fs::remove(worker_log_path, ec);
    ec.clear();
    fs::remove(worker_csv_path, ec);
    ec.clear();
    fs::remove(worker_update_csv_path, ec);
    return result;
}

static bool file_is_empty_or_missing(const string& path) {
    std::error_code ec;
    if (!fs::exists(path, ec)) return true;
    return fs::file_size(path, ec) == 0 || ec;
}

static string read_first_line(const string& path) {
    ifstream in(path);
    string line;
    if (!in.is_open() || !getline(in, line)) return {};
    return line;
}

static bool contains_text(const string& value, const string& needle) {
    return value.find(needle) != string::npos;
}

static bool rebuild_all_results_csv(const string& density_path,
                                    const string& k_path,
                                    const string& all_path) {
    const string tmp_path = all_path + ".tmp";
    ifstream density_in(density_path);
    ifstream k_in(k_path);
    if (!density_in.is_open() || !k_in.is_open()) return false;
    ofstream out(tmp_path, ios::out | ios::trunc);
    if (!out.is_open()) return false;

    string line;
    bool got_header = false;
    if (getline(density_in, line)) {
        out << line << '\n';
        got_header = true;
        while (getline(density_in, line)) if (!line.empty()) out << line << '\n';
    }
    if (!getline(k_in, line)) {
        out.close();
        std::error_code ec;
        fs::remove(tmp_path, ec);
        return false;
    }
    if (!got_header) out << line << '\n'; // density file was unexpectedly empty
    while (getline(k_in, line)) if (!line.empty()) out << line << '\n';
    out.flush();
    out.close();

    std::error_code ec;
    fs::rename(tmp_path, all_path, ec);
    if (ec) {
        // Cross-filesystem rename is not expected, but fall back to replacement.
        fs::remove(all_path, ec);
        ec.clear();
        fs::rename(tmp_path, all_path, ec);
    }
    return !ec;
}

int main(int argc, char** argv) {
    Options opt;
    if (!parse_options(argc, argv, opt)) {
        return (argc > 1 && (string(argv[1]) == "--help" || string(argv[1]) == "-h")) ? 0 : 1;
    }

    fs::create_directories(opt.output_dir);
    const string log_path = opt.output_dir + "/benchmark_v89_final.log";
    const string density_csv_path = opt.output_dir + "/density_results.csv";
    const string k_csv_path = opt.output_dir + "/k_results.csv";
    const string update_csv_path = opt.output_dir + "/update_results.csv";
    const string all_csv_path = opt.output_dir + "/all_results.csv";
    const string status_csv_path = opt.output_dir + "/run_status.csv";

    // The query CSV schema intentionally removed all profiling-only fields.
    // Refuse to append to the previous schema, otherwise resume would silently
    // mix rows with different column counts and may skip cases from old runs.
    if (opt.resume) {
        const vector<string> query_paths = {density_csv_path, k_csv_path, all_csv_path};
        for (const string& path : query_paths) {
            if (file_is_empty_or_missing(path)) continue;
            const string header = read_first_line(path);
            if (!contains_text(header, "\"clean_query_us\"") ||
                contains_text(header, "profiling_overhead_pct") ||
                contains_text(header, "local_2hop_us")) {
                cerr << "ERROR: incompatible query CSV schema in " << path << '\n'
                     << "Use a fresh --output-dir for the Scheme C clean-query benchmark.\n";
                return 2;
            }
        }
        if (!file_is_empty_or_missing(update_csv_path)) {
            const string header = read_first_line(update_csv_path);
            if (!contains_text(header, "index_maintenance_pending_at_end")) {
                cerr << "ERROR: incompatible update CSV schema in " << update_csv_path << '\n'
                     << "Use a fresh --output-dir for the Scheme C clean-query benchmark.\n";
                return 2;
            }
        }
    }

    unordered_set<string> completed_cases;
    if (opt.resume) completed_cases = load_completed_cases(status_csv_path);

    const auto open_mode = opt.resume ? (ios::out | ios::app) : (ios::out | ios::trunc);
    ofstream log(log_path, open_mode);
    CSVWriter density_csv(density_csv_path, opt.resume);
    CSVWriter k_csv(k_csv_path, opt.resume);
    CSVWriter update_csv(update_csv_path, opt.resume);
    CSVWriter all_csv(all_csv_path, opt.resume);
    ofstream status_out(status_csv_path, open_mode);

    if (!log.is_open() || !density_csv.good() || !k_csv.good() ||
        !update_csv.good() || !all_csv.good() || !status_out.is_open()) {
        cerr << "ERROR: failed to open one or more output files under "
             << opt.output_dir << '\n';
        return 2;
    }

    write_query_header_if_needed(density_csv, file_is_empty_or_missing(density_csv_path));
    write_query_header_if_needed(k_csv, file_is_empty_or_missing(k_csv_path));
    write_query_header_if_needed(all_csv, file_is_empty_or_missing(all_csv_path));
    write_update_header_if_needed(update_csv, file_is_empty_or_missing(update_csv_path));
    write_status_header_if_needed(status_out, file_is_empty_or_missing(status_csv_path));

    streambuf* original_cout = cout.rdbuf(log.rdbuf());
    cout << "============================================================\n"
         << "V8.9 SCHEME C FINAL BENCHMARK DRIVER\n"
         << "============================================================\n"
         << "Base Path                  = " << opt.base_path << '\n'
         << "Output Directory           = " << opt.output_dir << '\n'
         << "Phase                      = " << opt.phase << '\n'
         << "Cell Capacity Factor      = " << TNR_2hop::CELL_CAPACITY_FACTOR << "K\n"
         << "Local Fast2Hop Height      = " << TNR_2hop::LOCAL_F2H_HEIGHT << '\n'
         << "Query Algorithm            = Scheme C Two-Pointer Hub Merge + Object-Coalesced K-way Merge\n"
         << "Query Warmup               = " << opt.warmup << '\n'
         << "Per-case sequence           = Query -> Insert -> Delete -> one Apply_Index_Updates()\n"
         << "Execution Isolation        = one child process per experiment configuration\n"
         << "Resume Mode                = " << (opt.resume ? "ENABLED" : "DISABLED") << '\n'
         << "Host RAM                   = " << sample_memory().host_ram_gb << " GiB\n"
         << "Cgroup Memory Limit        = " << sample_memory().cgroup_limit_gb << " GiB (-1 means unavailable/unlimited)\n"
         << "Datasets                   = ";
    for (size_t i = 0; i < opt.datasets.size(); ++i) {
        if (i) cout << ',';
        cout << opt.datasets[i];
    }
    cout << "\nDensity Values             = ";
    for (size_t i = 0; i < opt.densities.size(); ++i) {
        if (i) cout << ',';
        cout << opt.densities[i];
    }
    cout << "\nK Values                   = ";
    for (size_t i = 0; i < opt.k_labels.size(); ++i) {
        if (i) cout << ',';
        cout << opt.k_labels[i];
    }
    cout << "\n";
    log.flush();

    unsigned long long case_serial = 0;
    long long failed_cases = 0;
    long long skipped_cases = 0;
    const auto total_start = steady_clock::now();

    auto execute_query = [&](const string& exp, const string& dataset,
                             const string& density, const string& k) {
        IsolatedRunResult r = run_isolated_case(opt, exp, dataset, density, k,
                                                exp == "density" ? density_csv : k_csv,
                                                nullptr, update_csv, log, status_out,
                                                completed_cases, case_serial);
        if (r.state == "SKIPPED") ++skipped_cases;
        else if (r.state != "OK" && r.state != "RESUME_SKIPPED") ++failed_cases;
    };

    if (opt.phase == "all" || opt.phase == "density") {
        cout << "\n========== Phase 1: Fixed K=" << opt.fixed_k
             << ", Varying Density; each case includes query + insert + delete + one Apply (isolated workers) ==========\n";
        for (const string& dataset : opt.datasets) {
            for (const string& density : opt.densities) {
                cout.flush();
                execute_query("density", dataset, density, opt.fixed_k);
            }
        }
    }

    if (opt.phase == "all" || opt.phase == "k") {
        cout << "\n========== Phase 2: Fixed Density=" << opt.fixed_density
             << ", Varying K; each case includes query + insert + delete + one Apply (isolated workers) ==========\n";
        for (const string& dataset : opt.datasets) {
            for (const string& k_label : opt.k_labels) {
                cout.flush();
                execute_query("k", dataset, opt.fixed_density, k_label);
            }
        }
    }

    const double total_sec = duration<double>(steady_clock::now() - total_start).count();
    cout << "\n============================================================\n"
         << "V8.9 MEMORY-SAFE BENCHMARK FINISHED\n"
         << "Total wall time = " << fixed << setprecision(6) << total_sec << " s\n"
         << "Failed cases   = " << failed_cases << '\n'
         << "Skipped cases  = " << skipped_cases << '\n'
         << "Log            = " << log_path << '\n'
         << "Density CSV    = " << density_csv_path << '\n'
         << "K CSV          = " << k_csv_path << '\n'
         << "Update CSV     = " << update_csv_path << '\n'
         << "Status CSV     = " << status_csv_path << '\n'
         << "============================================================\n";
    cout.flush();
    cout.rdbuf(original_cout);
    log.close();
    density_csv.out.close();
    k_csv.out.close();
    update_csv.out.close();
    all_csv.out.close();
    status_out.close();
    const bool aggregate_ok = rebuild_all_results_csv(density_csv_path, k_csv_path, all_csv_path);
    if (!aggregate_ok) {
        cerr << "WARNING: failed to rebuild aggregate CSV: " << all_csv_path << '\n';
        ++failed_cases;
    }
    cout << "Finished. Inspect run_status.csv for per-case peak RSS and termination causes.\n";
    return failed_cases == 0 ? 0 : 3;
}
