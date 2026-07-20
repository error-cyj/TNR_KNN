# CH_TNR: Dynamic k-Nearest Neighbor (kNN) Query Framework

## Prerequisites
* **C++17** compatible compiler (GCC 9+, Clang, MSVC)
* **CMake 3.15** or higher
* *(Optional)* **Jemalloc** allocator for maximum memory performance

## Data Preparation
Before running the program, ensure the following **four** input files are placed in your **Working Directory** (usually the project root folder):

1. **`graph.txt`**: The original graph data (e.g., rename your `Ein.txt` to `graph.txt`). The first line must be `Number_of_Nodes Number_of_Edges`.
2. **`objects.txt`**: Initial object distribution. (The program will automatically fix the format if the node count is missing from the first line).
3. **`queries.txt`**: Query workload file. Each line should contain `Source_Node_ID K_value`.
4. **`updates.txt`**: Dynamic update workload file. Each line should contain `Deleted_Old_Object_ID Added_New_Object_ID`.

## Compilation & Build
Open your terminal and run the following commands to build the project using CMake:

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
cmake --build . -j 10
```

## Execution
After a successful build, you can run the executable:

```bash
# On Windows
./CH_TNR.exe

# On Linux / macOS
./CH_TNR
```

*Note:* Standard console output will be redirected and saved to **`Runtime_Log.txt`** to eliminate I/O overhead during performance benchmarking. You can check this log file for detailed execution traces and latency metrics (measured in microseconds `us`).
