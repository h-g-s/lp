#!/bin/bash
./buildCBC.sh ; ./test
./buildGLPK.sh ; ./test
./buildCPX.sh ; ./test
./buildGRB.sh ; ./test
