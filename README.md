LP     {#mainpage}
==

LP is a simple and lightweight solver-independent C/C++ wraper to
Mixed-Integer Programming (MIP) solvers. It currently supports COIN-OR
CBC, GLPK, CPLEX and GUROBI. It can be called from any C 89 (or higher)
compatible compiler. I created LP because as someone used to the C based
CPLEX callable libeary, I was frustrated with the complex object oriented
framework that one had to master to successfully use the COIN-OR CBC
integer optimizer (OSI, CGL, CBC...).

Besides the standard
[Creation/Modification](group__groupCreateMod.html) and
[Optimization](group__groupOpt.html) routines LP provices access to
some ufesull features to Optimize Mixed Integer Programs: 

- Handling of solution pool
- MIPStart (entering an initial feasible solution)
- Callback for cut generation (currently only in CBC)

## Building

Just add lp.cpp and lp.h to your project and specify your mip solver
adding one of the following compilation directives (ex. in `GCC` `-DCBC`): 
- `CBC` to build using the COIN-OR Branch-and-Cut solver
- `GLPK` to build using the GNU Linear Programming Kit (GLPK) solver
- `CPX` to build using the IBM CPLEX MIP (c) Solver
- `GRB` to build using the Gurobi (c) MIP Solver

