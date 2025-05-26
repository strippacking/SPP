# Efficient SAT and MaxSAT Techniques for the Two-Dimensional Strip Packing Problem

This repository contains the implementation and experimental evaluation of SAT and MaxSAT-based approaches for solving the Two-Dimensional Strip Packing Problem (2SPP), as presented in the paper:

**"Efficient SAT and MaxSAT techniques for solving the Two-Dimensional Strip Packing Problem"**  
*Authors: Tuyen Kieu Van, Duong Quy Le, and Khanh Van To*  
*VNU University of Engineering and Technology, Hanoi, Vietnam*  
*Submitted to: Pesquisa Operacional*

## Abstract

The NP-hard Two-Dimensional Strip Packing Problem (2SPP) demands efficient exact solutions. This paper presents SAT-based Order Encoding models incorporating item rotation and adapted symmetry-breaking (SB). It compares three height-minimization strategies: Non-Incremental SAT Bisection, Incremental SAT Bisection, and Direct Weighted Partial MaxSAT Optimization. Benchmarks on established 2SPP instances show our Non-Incremental SAT and Direct MaxSAT strategies significantly outperform earlier incremental techniques.

## Repository Structure

```
SPP/
├── README.md                    # This file
├── ws-apjor.tex                # Paper source (LaTeX)
├── SPP_Result.xlsx             # Complete experimental results summary
├── runlim                      # Time limit utility
├── run.sh                      # Batch execution script
├── tt-open-wbo-inc-Glucose4_1_static  # MaxSAT solver binary
│
├── inputs/                     # Dataset instances (41 benchmark problems)
│   ├── ins-1.txt              # HT01(c1p1) instance
│   ├── ins-2.txt              # HT02(c1p2) instance
│   ├── ...
│   └── ins-41.txt             # HT12(c4p3) instance
│
├── SAT/MaxSAT Implementations:
│   ├── SPP.py                 # Non-Incremental SAT (No SB)
│   ├── SPP_SB_C1.py          # Non-Incremental SAT (SB Config 1)
│   ├── SPP_SB_C2.py          # Non-Incremental SAT (SB Config 2)
│   ├── SPP_INC_SB_C1.py      # Incremental SAT (SB Config 1)
│   ├── SPP_INC_SB_C2.py      # Incremental SAT (SB Config 2)
│   ├── SPP_MS.py             # Direct MaxSAT (No SB)
│   ├── SPP_MS_SB_C1.py       # Direct MaxSAT (SB Config 1)
│   ├── SPP_MS_SB_C2.py       # Direct MaxSAT (SB Config 2)
│   ├── SPP_R.py              # Rotational SAT (No SB)
│   ├── SPP_R_SB.py           # Rotational SAT (with SB)
│   ├── SPP_INC_R_SB.py       # Incremental Rotational SAT
│   ├── SPP_MS_R.py           # Direct MaxSAT Rotational (No SB)
│   └── SPP_MS_R_SB.py        # Direct MaxSAT Rotational (with SB)
│
├── CP/MIP Solver Implementations:
│   ├── CPLEX_CP_C2.py        # CPLEX CP Optimizer (SB Config 2)
│   ├── CPLEX_CP_R_SB.py      # CPLEX CP Optimizer (Rotational)
│   ├── CPLEX_MIP_C2.py       # CPLEX MIP Solver (SB Config 2)
│   ├── CPLEX_MIP_R_SB.py     # CPLEX MIP Solver (Rotational)
│   ├── GUROBI_MIP_C2.py      # Gurobi MIP Solver (SB Config 2)
│   ├── GUROBI_MIP_R_SB.py    # Gurobi MIP Solver (Rotational)
│   ├── OR-TOOLS_CP_C2.py     # OR-Tools CP-SAT (SB Config 2)
│   ├── OR-TOOLS_CP_R_SB.py   # OR-Tools CP-SAT (Rotational)
│   ├── OR-TOOLS_MIP_C2.py    # OR-Tools MIP (SB Config 2)
│   └── OR-TOOLS_MIP_R_SB.py  # OR-Tools MIP (Rotational)
│
└── Results/                   # Output folders for each configuration
    ├── SPP_SB_C2/            # Results and visualizations
    ├── CPLEX_CP_C2/
    ├── OR-TOOLS_CP_C2/
    └── ... (other solver results)
```

## Dataset Description

This repository uses a comprehensive benchmark dataset of **41 standard 2SPP instances** aggregated from well-established sources in the strip packing literature. The dataset provides a robust testbed with varying problem sizes and characteristics.

### Dataset Composition

The benchmark instances are drawn from the following established datasets:

| Dataset | Instances | Source | Description |
|---------|-----------|---------|-------------|
| **HT** | 12 instances | Hopper & Turton (2001) | Standard benchmark instances with varying sizes |
| **BENG** | 10 instances | Bengtsson (1982) | Classical packing instances |
| **CGCUT** | 3 instances | Christofides & Whitlock (1977) | Cutting and packing problems |
| **GCUT** | 4 instances | Beasley (1985) | Large-scale instances |
| **NGCUT** | 12 instances | Beasley (1985) | Extended cutting problems |

### Instance Characteristics

- **Problem sizes**: Range from 7 to 200 rectangles
- **Strip widths**: Vary from 10 to 250 units
- **Rectangle dimensions**: Diverse width and height combinations
- **Optimal heights**: Known optimal solutions from literature for comparison

### Instance Naming Convention

| File | Instance Name | Items (n) | Width (W) | Lower Bound | Optimal Height |
|------|---------------|-----------|-----------|-------------|----------------|
| ins-1.txt | HT01(c1p1) | 16 | 20 | 20 | 20 |
| ins-2.txt | HT02(c1p2) | 17 | 20 | 20 | 20 |
| ins-3.txt | HT03(c1p3) | 16 | 20 | 20 | 20 |
| ... | ... | ... | ... | ... | ... |
| ins-29.txt | BENG01 | 20 | 25 | 30 | 30 |
| ins-30.txt | BENG02 | 40 | 25 | 57 | 57 |
| ... | ... | ... | ... | ... | ... |
| ins-41.txt | HT12(c4p3) | 49 | 60 | 60 | 60 |

### Instance File Format

Each instance file (`ins-X.txt`) follows this format:
```
W          # Strip width
n          # Number of rectangles
w1 h1      # Width and height of rectangle 1
w2 h2      # Width and height of rectangle 2
...
wn hn      # Width and height of rectangle n
```

### Dataset Sources and References

The instances were obtained from:
- **University of Bologna OR Group library**: Primary source for benchmark instances
- **Literature references**: 
  - Hopper, E., & Turton, B. C. (2001). An empirical investigation of meta-heuristic and heuristic algorithms for a 2D packing problem. *European Journal of Operational Research*
  - Bengtsson, B. E. (1982). Packing rectangular pieces—A heuristic approach. *The Computer Journal*
  - Christofides, N., & Whitlock, C. (1977). An algorithm for two-dimensional cutting problems. *Operations Research*
  - Beasley, J. E. (1985). Algorithms for unconstrained two-dimensional guillotine cutting. *Journal of the Operational Research Society*

### Optimal Solutions

The repository includes known optimal heights for both:
- **Non-rotational variant**: Fixed orientation of rectangles
- **Rotational variant**: 90-degree rotation allowed

These optimal values serve as targets for solution quality evaluation and are primarily sourced from:
- Kenmochi et al. (2009)
- Arahori et al. (2012)
- Wei & Liu (2011)
- Alvarez-Valdés et al. (2009)

**Note:** Detailed comparison of achieved vs. optimal heights for all configurations is available in `SPP_Result.xlsx`.  
- Wei & Liu (2011)
- Alvarez-Valdés et al. (2009)

## Key Features

### SAT/MaxSAT Methodologies

1. **Order Encoding**: Compact representation of geometric constraints using Boolean variables
2. **Rotation Support**: Explicit modeling of 90-degree rectangle rotations
3. **Symmetry Breaking**: Adapted constraints for rotational variants
4. **Multiple Strategies**:
   - Non-Incremental SAT Bisection
   - Incremental SAT Bisection  
   - Direct Weighted Partial MaxSAT Optimization

### Solver Integration

- **SAT Solver**: Glucose 4.2 for decision problems
- **MaxSAT Solver**: TT-Open-WBO-Inc for optimization
- **Comparison Baselines**: CPLEX, Gurobi, Google OR-Tools (CP and MIP)

### Symmetry Breaking Configurations

- **SB C1**: Domain Reduction + Large Rectangles + Same Rectangles
- **SB C2**: Large Rectangles + One Pair of Rectangles

## Requirements

### Dependencies
```bash
pip install pysat
pip install matplotlib
pip install pandas
pip install openpyxl
```

### External Solvers
- **Glucose 4.2**: Included in pysat package
- **TT-Open-WBO-Inc**: Binary included (`tt-open-wbo-inc-Glucose4_1_static`)
- **Commercial Solvers** (optional): CPLEX, Gurobi
- **OR-Tools**: `pip install ortools`

## Usage

### Running Individual Instances

```bash
# Non-Incremental SAT with Symmetry Breaking C2
python3 SPP_SB_C2.py 1  # Run instance 1 (HT01)

# Rotational SAT with Symmetry Breaking
python3 SPP_R_SB.py 1   # Run instance 1 with rotation

# Direct MaxSAT
python3 SPP_MS.py 1     # Direct MaxSAT optimization
```

### Batch Execution

```bash
# Run all instances with controller mode
python3 SPP_SB_C2.py    # Runs all 41 instances sequentially

# Using runlim for timeout control
./runlim --time-limit=1800 python3 SPP_SB_C2.py 1
```

### Results

Results are automatically saved to:
- **Excel files**: `[ConfigName].xlsx` with performance metrics for each configuration
- **Visualizations**: PNG files in respective result folders showing optimal packings
- **JSON checkpoints**: For timeout recovery and intermediate results
- **Complete Summary**: `SPP_Result.xlsx` - Comprehensive results across all methods

### SPP_Result.xlsx - Complete Experimental Results

The `SPP_Result.xlsx` file contains a comprehensive summary of all experimental results from our paper. This Excel file includes:

**Content Overview:**
- **Performance comparison** across all 22 solver configurations tested
- **Runtime analysis** for each of the 41 benchmark instances
- **Success rates** (number of instances solved optimally within time limit)
- **Optimal heights** achieved by each method for both rotational and non-rotational variants
- **Statistical analysis** and comparative metrics between SAT/MaxSAT and CP/MIP approaches

**Key Worksheets:**
- **Summary Sheet**: Overall performance ranking and success rates
- **Runtime Comparison**: Detailed execution times for all configurations
- **Height Analysis**: Optimal solutions found vs. known optimal values
- **Method Comparison**: Direct comparison between SAT, MaxSAT, CP, and MIP approaches
- **Instance Analysis**: Per-instance breakdown showing which methods succeeded/failed

**Usage:**
This file serves as the primary reference for reproducing the results presented in our paper and provides researchers with detailed performance data for comparative studies.

## Experimental Results Summary

Based on our comprehensive evaluation of 41 benchmark instances:

### Best Performing Configurations

| Method | Configuration | Instances Solved | Key Advantage |
|--------|---------------|------------------|---------------|
| Direct MaxSAT | SPP_MS, SPP_MS_SB_C2 | 26/41 | Highest success rate (non-rotational) |
| Non-Inc SAT | SPP_R_SB | 26/41 | Best for rotational problems |
| Non-Inc SAT | SPP_SB_C2 | 24/41 | Excellent runtime efficiency |
| OR-Tools CP | OR-TOOLS_CP_C2 | 28/41 | **Best overall performance** |

### Key Findings

1. **Non-Incremental SAT** and **Direct MaxSAT** significantly outperform Incremental SAT methods
2. **Google OR-Tools CP-SAT** achieves the highest number of optimal solutions (28/41)
3. **SAT/MaxSAT methods** are highly competitive, outperforming other CP solvers and all tested MIP solvers
4. **Symmetry breaking** generally improves performance
5. **Item rotation** increases encoding complexity but can yield better optimal heights for some instances

## Citation

If you use this code or dataset in your research, please cite:

<!-- ```bibtex
@article{van2024efficient,
  title={Efficient SAT and MaxSAT techniques for solving the Two-Dimensional Strip Packing Problem},
  author={Van, Tuyen Kieu and Le, Duong Quy and To, Khanh Van},
  journal={Pesquisa Operacional},
  year={2024},
  publisher={Brazilian Operations Research Society}
}
``` -->

## License

This project is made available for research and educational purposes. Please refer to individual solver licenses for commercial use restrictions.

## Contact

- **Tuyen Kieu Van**: tuyenkv@vnu.edu.vn
- **Duong Quy Le**: 21020560@vnu.edu.vn  
- **Khanh Van To**: khanhtv@vnu.edu.vn

**VNU University of Engineering and Technology**  
Hanoi, Vietnam

## Acknowledgments

This work has been supported by VNU University of Engineering and Technology under project number CN24.10.

The benchmark instances are sourced from established datasets in the strip packing literature, and we acknowledge the original authors and the University of Bologna OR Group for maintaining these valuable resources.
