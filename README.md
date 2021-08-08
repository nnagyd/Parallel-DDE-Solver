# Parallel DDE Solver
 Parallel Delay Differential Equation solver on CPU using AVX Vector instructions, based on [VCL](https://github.com/vectorclass/version2) library.

# Requirements
- C++14 or C++17 compiler
- [VCL](https://github.com/vectorclass/version2) library version 2

# Performance

- Faster than most other codes

- Using AVX instruction set and built-in loop unrolling

## Performance on the delayed logistic equation

![alt text](https://github.com/nnagyd/Parallel-DDE-Solver/blob/master/performance/logistic.png)

## Performance on the delayed Lorenz equation

![alt text](https://github.com/nnagyd/Parallel-DDE-Solver/blob/master/performance/lorenz.png)
