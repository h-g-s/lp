LP     {#mainpage}
==

LP is a simple and lightweight solver-independent C/C++ wraper to
Mixed-Integer Programming (MIP) solvers. It currently supports COIN-OR
CBC, GLPK, CPLEX and GUROBI. It can be called from any C 89 (or higher)
compatible compiler. I created lp because as someone used to the C based
CPLEX callable libeary, I was frustrated with the complex object oriented
framework that one had to master to successfully use the COIN-OR CBC
integer optimizer (OSI, CGL, CBC...).

## Building

Just add lp.cpp and lp.h to your project and specify your mip solver
adding one of the following compilation directives: 
    - -DCBC to build using the COIN-OR Branch-and-Cut solver
    - -DGLPK to build using the GNU Linear Programming Kit (GLPK) solver
    - -DCPX to build using the IBM CPLEX MIP Solver
    - -DGRB to build using the Gurobi MIP Solver

