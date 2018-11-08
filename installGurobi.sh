g++ -shared -DGRB -Og -g -DDEBUG -fPIC\
    -L/opt/gurobi801/linux64/lib/ \
    -I/opt/gurobi801/linux64/include/ \
    lp.cpp -lgurobi80 -o lp-grb-linux64.so
sudo cp lp-grb-linux64.so /usr/lib/
sudo ldconfig
