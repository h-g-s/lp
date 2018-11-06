from lp import *
import sys
from typing import List

if len(sys.argv)==1:
	print('enter number of queens.')
	exit()

n = int(sys.argv[1])

mip = LinearProgram()

obj : List[int] = [1.0]*(n*n)
lb : List[float] = [0.0]*(n*n)
ub : List[float] = [1.0]*(n*n)
integer : List[bool] = [True]*(n*n)
cnames : List[str] = []

x =	[[0 for i in range(n)] for j in range(n)]

p = 0
for i in range(n):
	for j in range(n):
		cnames.append('x({},{})'.format(i,j))
		x[i][j] = p
		p += 1

mip.add_cols(n*n, obj, lb, ub, integer, cnames)

# one per line
for i in range(n):
	idx = [ x[i][j] for j in range(n) ]
	coef = [1.0]*n
	mip.add_row( idx, coef, 'line({})'.format(i), Sense.EQUAL, 1.0 )

for i in range(n):
	idx = [ x[j][i] for j in range(n) ]
	coef = [1.0]*n
	mip.add_row( idx, coef, 'column({})'.format(i), Sense.EQUAL, 1.0 )

ridx : List[int] = [0,1]
rcoef : List[float] = [1,1]
mip.add_row(ridx, rcoef, 'c1', Sense.EQUAL, 1.0 )
mip.write('test.lp')
