# WAPA policy: A Microarchitecture- and Workload-Agnostic Universal SMT Scheduler

This repository contains the C++ implementation of the **WAPA policy**. WAPA is a T2C allocation policy for SMT processors grounded in Optimal Transport theory [[1]](https://dl.acm.org/doi/10.1007/978-3-031-99854-6_16). 

## ⚠️ Important Implementation Note

**Please note:** This code is not a standalone executable. It was developed to run inside a **proprietary custom scheduling framework**. 

* The code relies on specific data structures (e.g., `tasklist_t`, `Task`) and logging macros (`LOGINF`) provided by that host environment.
* It is provided here for research, educational purposes, and to demonstrate the algorithmic implementation of the WAPA policy and its mathematical auxiliary methods.
