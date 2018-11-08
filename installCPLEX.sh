g++ -shared -DCPX -Og -g -DDEBUG -fPIC\
    -L/opt/ibm/ILOG/CPLEX_Studio128/cplex/lib/x86-64_linux/static_pic/ \
    -I/opt/ibm/ILOG/CPLEX_Studio128/cplex/include/ilcplex/ \
    lp.cpp -lcplex -o lp-cpx-linux64.so
sudo cp lp-cpx-linux64.so /usr/lib/ 
sudo ldconfig
