g++ -g3  -fsanitize=address -DGRB ../lp.cpp test.cpp \
    -I/opt/gurobi752/linux64/include/ \
    -L/opt/gurobi752/linux64/lib/ -lgurobi75 -lm -o test
