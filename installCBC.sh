g++ -shared -DCBC -Og -g -DDEBUG -fPIC\
    `pkg-config --libs cbc` \
    `pkg-config --cflags cbc` \
    lp.cpp -o lp-cbc-linux64.so
sudo cp lp-cbc-linux64.so /usr/lib/
sudo ldconfig
