# QPSolver

![Unit Tests](https://github.com/govindchari/QPSolver/actions/workflows/unit_tests.yml/badge.svg)


A C++ solver for quadratic programs that utilizes Mehrotra's predictor-corrector algorithm.

Based on cvxgen paper includes
  - -O3 compiler optimization
  - Permuted LDLT decomposition
  - Iterative refinement
  - Sparse matrices

Debug:
  - Add #include <string.h> to top of QPSolver/tests/Epigraph/solvers/ecos/ecos_bb/ecos_bb.c
