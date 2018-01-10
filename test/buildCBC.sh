g++ -DCBC \
    test.cpp ../lp.cpp  \
    -Wall -g3 \
    `pkg-config --cflags cbc` -o test `pkg-config --libs cbc`
