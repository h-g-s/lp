from lp import *
import sys
from typing import List

if len(sys.argv)==1:
	print('enter number of queens.')
	exit()

n = int(sys.argv[1])

mip = LinearProgram()

obj : List[float] = [1.0]*(n*n)
lb  : List[float] = [0.0]*(n*n)
ub  : List[float] = [1.0]*(n*n)
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

# one per column
for i in range(n):
	idx = [ x[j][i] for j in range(n) ]
	coef = [1.0]*n
	mip.add_row( idx, coef, 'column({})'.format(i), Sense.EQUAL, 1.0 )

# \ diagonal
p=0
for k in range(2-n,n-2+1):
	idx = [x[i][j] for i in range(0,n) for j in range(0,n) if i-j==k]
	coef = [1.0] * len(idx)
	mip.add_row(idx, coef, 'selDiag1({})'.format(p), Sense.LESS_OR_EQUAL, 1.0)
	p += 1

# / diagonal
p=0
for k in range(3,n+n):
	idx = [x[i][j] for i in range(0,n) for j in range(0,n) if i+j==k]
	coef = [1.0] * len(idx)
	mip.add_row(idx, coef, 'selDiag2({})'.format(p), Sense.LESS_OR_EQUAL, 1.0)
	p += 1

mip.optimize()

sol = mip.x()

for i in range(0,mip.cols()):
	print( 'cname {}'.format(mip.col_name(i)) )

for i in range(0, n):
	for j in range(0, n):
		if abs(sol[x[i][j]])>=0.99:
			sys.stdout.write(' O')
		else:
			sys.stdout.write(' .')

	sys.stdout.write('\n')


mip.write('test.lp')

