# Scheme C (V8.9) — Reviewer Reproducibility Artifact

This package contains the **frozen Scheme C implementation** of the dynamic road-network Top-K nearest-object search algorithm and the original single-threaded benchmark driver. The algorithm is **not simplified or rewritten** for this artifact: the `src/` code and driver were copied from `V8_9_FINAL_SCHEME_C_EXPERIMENT.zip` to preserve experimental equivalence.

## Artifact status and included material

- `src/`: the six core C++17 implementation files: `ET.*`, `Graph.h`, `TNR_2hop.*`, `Fast2Hop.h`.
- `main_v89_final_benchmark.cpp`: exact frozen experiment driver (source filename retained).
- `test_v89_correctness.cpp`, `test_tnr_end_to_end.cpp`: standalone correctness checks, including a generated 40-vertex synthetic graph.
- `CMakeLists.txt`: release build and test definitions.
- `run_ny.sh`: NY-only density and/or K sweep wrapper with input checks.
- `results/provided/density_results_original.csv`: an **earlier provided** density-sweep CSV for provenance, **NOT final paper results**: its CAL 0.001/0.0001 timings were subsequently reported to be affected by server interference.
- `results/provided/original_validation_results.txt`: upstream validation notes.

**Not currently supplied:** the actual NY graph (`data/NY/Ein.txt`), the NY object/query/update datasets, and the corrected final CAL results. The package cannot honestly be called a complete real-NY reproduction artifact until the exact files used in the experiments are added. In particular, DO NOT substitute an arbitrary DIMACS New York graph for the original `Ein.txt`: preprocessing and timings depend on the graph content and indexing.

## Platform and build

Linux or WSL2, GCC/Clang supporting C++17, CMake >=3.16, adequate RAM for road-network preprocessing. The benchmark runner uses POSIX `fork`, `wait4`, and `/proc`; native Windows is not supported. Query/index execution is **single-threaded** (CMake `-j` controls build jobs only).

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_V89_CORRECTNESS_TEST=ON
cmake --build build -j2
./build/test_v89_correctness
./build/test_tnr_end_to_end
```

Expected correctness markers include `V8_9_DETERMINISTIC_OK`, `V8_9_HIGH_GLOBAL_RANK_OK`, `V8_9_RANDOMIZED_500_OK`. The end-to-end test compares Top-K results to Dijkstra before/after object updates. Tied objects at the Kth distance boundary may differ.

`-march=native` / LTO are enabled in the original benchmark CMake release configuration. **Record the compiler and CPU model**, since speed numbers may differ across machines.

## Original NY dataset layout

Place the **exact original input files** under this directory:

```text
data/NY/Ein.txt
Object data/NY/data_0.1.txt
Object data/NY/data_0.01.txt
Object data/NY/data_0.001.txt
Object data/NY/data_0.0001.txt
Queary data/NY/_k20_Query.txt
Queary data/NY/_k10_Query.txt
Queary data/NY/_k30_Query.txt
Queary data/NY/_k40_Query.txt
Update data/NY/data_0.1_Ins.txt
Update data/NY/data_0.1_Del.txt
Update data/NY/data_0.01_Ins.txt
Update data/NY/data_0.01_Del.txt
Update data/NY/data_0.001_Ins.txt
Update data/NY/data_0.001_Del.txt
Update data/NY/data_0.0001_Ins.txt
Update data/NY/data_0.0001_Del.txt
```

**Preserve the spelling** `Queary data`, as it is expected by the original driver.

Graph input is a plain-text, **1-based, undirected weighted edge list**:
first line `N M`, followed by edge lines `u v w` (integer endpoints and weights). The source reads all subsequent triples and inserts undirected edges. Do not assume this matches a raw `.gr` file from a public repository.

Query files: first line is the number of queries, then lines `source_vertex k`. Source vertex must be in `1..N` and K in `1..50`.

Update files: first line is number of operations, then `i vertex` for insertion files or `d vertex` for deletion files, using the same 1-based IDs. Missing update files are allowed by the driver but **do not reproduce update timing**.

Object files must be those prepared for this implementation and density; do not regenerate or change the random object assignment while claiming numerical reproduction.

## Running the NY experiments

```bash
# NY density sweep, K=20 (4 densities; updates are measured in each case)
bash run_ny.sh density

# NY K sweep, rho=0.001 (K=10,20,30,40; updates in each case)
bash run_ny.sh k

# Both sweeps
bash run_ny.sh all
```

Or invoke the frozen executable directly:

```bash
./build/v89_final_benchmark \
  --base "$PWD" --output-dir "$PWD/results/ny_density" \
  --phase density --datasets NY --warmup 100
```

The runner executes one sequential, isolated subprocess per configuration to release per-case memory. It **does not** run multiple benchmark cases concurrently.

## Timing and reported metrics

Each original query file supplies the queries: 100 warm-up records by default; the remaining records are measured exactly once in the clean query pass. `clean_query_us` is the per-query mean in **microseconds**. `preprocess_sec` is in **seconds**. `f2h_index_build_ms` is in **milliseconds** and includes local Cell-runtime initialization rather than isolated F2H index creation. The column `nonempty_cells` counts Cells with nonempty `local_nodes` (not necessarily Cells containing active objects).

For each case the order is:
1. preprocess and build the indices;
2. query warm-up and clean timed query pass;
3. process all inserts, then all deletes;
4. invoke `Apply_Index_Updates()` once.

`avg_update_time_us = avg_insert_us + avg_delete_us + apply_time_us`. In the original lazy policy, `Apply_Index_Updates()` may defer costly index reconstruction, so **do not call this a fully amortized rebuilding cost**. See the source for exact semantics.

Results appear in `results/ny_density` or `results/ny_k`: `density_results.csv`, `k_results.csv`, `update_results.csv`, `all_results.csv`, `run_status.csv`, and `benchmark_v89_final.log`. Preserve the full directory and its run status for independent checks.

## Algorithms / key entry points

- `src/ET.cpp`: weighted-road-network loading, elimination/tree construction, Cell/TN decomposition.
- `src/TNR_2hop.cpp`: F2H local index construction, query integration, TN global routing, dynamic object updates.
- `src/Fast2Hop.h`: compact descending `hub_ranks`, **Scheme C two-pointer Hub matching** in `query_1_to_N_hub_grouped`, and Object-Coalesced K-way candidate extraction.
- `main_v89_final_benchmark.cpp`: CLI, input loading, clean timing, update sequencing, process isolation, CSV reporting.

Scheme C matches descending source-label and Cell Hub rank streams in `O(L+H)` comparisons instead of allocating a per-Cell graph-sized reverse-rank array; the K-way merge retains one heap entry per candidate logical object as in the frozen implementation.

## Reproducibility caveats

- The included CSV is historical reference only; it contains superseded, abnormal CAL measurements. **Do not use it as a final-performance claim.**
- Re-run every included experimental row on a quiet machine and report median over independent **runs**, not the fastest run. Save CPU/compiler versions, concurrency, input checksums, per-run raw logs, and peak RSS.
- A NY-only artifact reproduces NY runs; it cannot independently reproduce other datasets without the corresponding exact data.
- This bundle preserves the author's current code; it is not a verified implementation of any third-party baseline.

## Provenance / verification

`SOURCE_SHA256SUMS.txt` records the hashes of the frozen source files. Running tests checks algorithmic correctness on synthetic cases; successful synthetic tests do not prove the same runtime on the unavailable real NY input.
