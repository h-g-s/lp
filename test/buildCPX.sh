g++ -g3  -fsanitize=address -DCPX ../lp.cpp test.cpp \
    -I/opt/ibm/ILOG/CPLEX_Studio1271/cplex/include/ilcplex/ \
    -L/opt/ibm/ILOG/CPLEX_Studio1271/cplex/lib/x86-64_linux/static_pic/ -lcplex -lm -lpthread -o test
