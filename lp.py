import ctypes as ct
from enum import Enum
from typing import List

class OptimizationStatus(Enum):
	Optimal         = 0
	Infeasible      = 1
	Unbounded       = 2
	Feasible        = 3
	IntInfeasible   = 4
	NoSolutionFound = 5
	Error           = 6

class Sense(Enum):
	LESS_OR_EQUAL    = 'L'
	GREATER_OR_EQUAL = 'G'
	EQUAL            = 'E'

class LinearProgram:
	""" low level class to create, modify and 
	optimize mixed-integer linear programs
	"""
	def __init__(self):
		self._plp = lp_create()


	def read(self, fileName : str):
		""" reads a Mixed-Integer Linear Program from an .lp or .mps file

		Args:
			fileName (str): file name
		"""
		lp_read( self._plp, fileName.encode('utf-8'))

	
	def write(self, fileName : str):
		lp_write( self._plp, fileName.encode('utf-8'))


	def add_col(self,
			obj      : float,      # objective function coefficient
			lb       : float,      # variable lower bound
			ub       : float,      # variable upper bound
			integer  : bool,       # variable upper bound
			name     : str,        # variable name
			nz       : int,        # number of non-zero entries for this column
			rowIdx   : List[int],  # index of rows where this column appears
			rowCoef  : List[float] # row coefficients for rows where this column appears
			):
		""" adds a new column (variable)
		"""
		if integer:
			charInt = ct.c_char(1)
		else:
			charInt = ct.c_char(0)

		vrowInt = (ct.c_int * nz)()
		for i in range(nz):
			vrowInt[i] = rowInt[i]

		vrowCoef = (ct.c_double * nz)()
		for i in range(nz):
			vrowCoef[i] = rowCoef[i]

		lp_add_col(self._plp, 
			ct.c_double(obj),
			ct.c_double(lb),
			ct.c_double(ub),
			charInt, name.encode('utf-8'),
			ct.c_int(len(rowIdx)), vrowInt, vrowCoef )


	def add_cols(self, 
		n : int, 
		obj : List[float],
		lb : List[float],
		ub : List[float],
		integer : List[bool],
		name : List[str] ):

		vobj = (ct.c_double*n)()
		vlb = (ct.c_double*n)()
		vub = (ct.c_double*n)()
		vint = (ct.c_char*n)()

		for i in range(n):
			vobj[i] = obj[i]

		for i in range(n):
			vlb[i] = lb[i]

		for i in range(n):
			vub[i] = ub[i]

		for i in range(n):
			vint[i] = ct.c_char( integer[i] )

		pnames = (ct.POINTER(ct.c_char)*n)()

		for i in range(0, n):
			pnames[i] = ct.create_string_buffer(name[i].encode('utf-8'))
			

		lp_add_cols(
			self._plp,
			ct.c_int(n),
			vobj,
			vlb,
			vub,
			vint,
			pnames )


	def add_row( self, 
			idx : List[int],
			coef : List[float], 
			name : str, 
			sense : str, rhs : float ):

		nz = len(idx)
		assert( len(idx)==len(coef) )

		vrowInt = (ct.c_int * nz)()
		for i in range(nz):
			vrowInt[i] = idx[i]

		vrowCoef = (ct.c_double * nz)()
		for i in range(nz):
			vrowCoef[i] = coef[i]

		if sense == Sense.EQUAL:
			chars = ct.c_char(ord('E'))
		elif sense == Sense.GREATER_OR_EQUAL:
			chars = ct.c_char(ord('G'))
		elif sense == Sense.LESS_OR_EQUAL:
			chars = ct.c_char(ord('L'))
		else:
			raise 'sense not recognized: {}'.format(sense)

		lp_add_row(
			self._plp,
			len(idx), vrowInt, vrowCoef, 
			ct.c_char_p(name.encode('utf-8')), chars, ct.c_double(rhs) )


	def cols(self) -> int:
		""" returns the number of columns (variables) in the
		mixed-integer linear program
		"""
		return lp_cols(self._plp)


	def rows(self) -> int:
		""" returns the number of rows (constraints) in the
		mixed-integer linear program
		"""
		return lp_rows(self._plp)


	def nz(self) -> int:
		"""returns the number of non-zeros in the constraint matrix"""
		return lp_nz(self._plp)

	
	def rhs(self, row : int) -> float:
		"""returns the right hand side of a given row"""
		return lp_rhs(self._plp, ct.c_int(row))


	def sense(self, row : int) -> str:
		"""returns the sense of a given constraint"""
		return chr(lp_sense(self._plp, ct.c_int(row))[0])


	def row_name(self, row : int) -> str:
		"""queries a row name"""
		return lp_row_name(self._plp, ct.c_int(row)).value


	def col_name(self, row : int) -> str:
		"""queries a col name"""
		return lp_col_name(self._plp, ct.c_int(row))


	def col_lb(self, col : int) -> float:
		"""queries a column lower bound"""
		return lp_col_lb(self._plp, ct.c_double(row))


	def col_ub(self, col : int) -> float:
		"""queries a column upper bound"""
		return lp_col_ub(self._plp, ct.c_double(row))


	def col_index(self, name : str) -> int:
		"""returns the column (variable) index of a given column name"""
		return lp_col_index(self._plp, ct.create_string_buffer(name.encode('utf-8')) )


	def row_index(self, name : str) -> int:
		"""returns the row (constraint) index of a given row name"""
		return lp_row_index(self._plp, ct.create_string_buffer(name.encode('utf-8')) )


	def obj_coef(self) -> List[float]:
		"""returns the objective function coefficients"""
		n = lp_cols(self._plp)
		objcoef = lp_obj_coef(self._plp)
		res : List[float] = [objcoef[i] for i in range(0,n)]
		return res


	def optimize(self, forceContinuous=False) -> int:
		"""  Optimizes your Mixed Integer Program. 
		Returns the problem status, which can be:
			0 : LP_OPTIMAL` : optimal solution found;
			1 : LP_INFEASIBLE` : the problem is infeasible; 
			2 : LP_UNBOUNDED` : the problem is unbounded;
			3 : LP_FEASIBLE` : a feasible solution was found, but its optimality was not proved;
			4 : LP_INTINFEASIBLE` : the lp relaxation is feasible but no integer feasible solution exists;
			5 : LP_NO_SOL_FOUND` :  optimization concluded without finding any feasible solution;
			6 : LP_ERROR` : the solver reported an error.

		Args:
			forceContinuous (str): if integrality constraints will be relaxed"""
		status = OptimizationStatus.Error
		if forceContinuous==False:
			status = lp_optimize(self._plp)
		else:
			status = lp_optimize_as_continuous(self._plp)

		return status


	def obj_value(self) -> float:
		"""objective value of your optimization """
		return lp_obj_value(self._plp)


	def x(self) -> List[float] :
		"""returns the vector of solution values for variables"""
		n = self.cols()
		res = [0.0]*n
		pr = lp_x(self._plp)
		for i in range(0,n):
			res[i] = pr[i]

		return res


	def __del__(self):
		lp_free(self._plp)

lplib = ct.CDLL('./lp-cbc-linux64.so')

lp_create = lplib.lp_create
lp_create.restype = ct.c_void_p

lp_read = lplib.lp_read
lp_read.argtypes = [ct.c_void_p, ct.c_char_p]

lp_write = lplib.lp_write_lp
lp_write.argtypes = [ct.c_void_p, ct.c_char_p]

lp_optimize = lplib.lp_optimize
lp_optimize.argtypes = [ct.c_void_p]
lp_optimize.restype = ct.c_int

lp_optimize_as_continuous = lplib.lp_optimize_as_continuous
lp_optimize_as_continuous.argtypes = [ct.c_void_p]
lp_optimize_as_continuous.restype = ct.c_int

lp_cols = lplib.lp_cols
lp_cols.argtypes = [ct.c_void_p]
lp_cols.restype = ct.c_int

lp_nz = lplib.lp_nz
lp_nz.argtypes = [ct.c_void_p]
lp_nz.restype = ct.c_int

lp_rhs = lplib.lp_rhs
lp_rhs.argtypes = [ct.c_void_p, ct.c_int]
lp_rhs.restype = ct.c_double

lp_sense = lplib.lp_sense
lp_sense.argtypes = [ct.c_void_p, ct.c_int]
lp_sense.restype = ct.c_char

lp_row_name = lplib.lp_row_name
lp_row_name.argtypes = [ct.c_void_p, ct.c_int]
lp_row_name.restype = ct.c_char_p

lp_col_name = lplib.lp_col_name
lp_col_name.argtypes = [ct.c_void_p, ct.c_int]
lp_col_name.restype = ct.c_char_p

lp_col_lb = lplib.lp_col_lb
lp_col_lb.argtypes = [ct.c_void_p, ct.c_int]
lp_col_lb.restype = ct.c_double

lp_col_ub = lplib.lp_col_ub
lp_col_ub.argtypes = [ct.c_void_p, ct.c_int]
lp_col_ub.restype = ct.c_double

lp_col_index = lplib.lp_col_lb
lp_col_index.argtypes = [ct.c_void_p, ct.c_char_p]
lp_col_index.restype = ct.c_int

lp_row_index = lplib.lp_row_index
lp_row_index.argtypes = [ct.c_void_p, ct.c_char_p]
lp_row_index.restype = ct.c_int

lp_obj_coef = lplib.lp_obj_coef
lp_obj_coef.argtypes = [ct.c_void_p]
lp_obj_coef.restype = ct.POINTER(ct.c_double)

lp_rows = lplib.lp_rows
lp_rows.argtypes = [ct.c_void_p]
lp_rows.restype = ct.c_int

lp_obj_value = lplib.lp_obj_value
lp_obj_value.argtypes = [ct.c_void_p]
lp_obj_value.restype = ct.c_float

lp_add_col = lplib.lp_add_col
lp_add_col.argtypes = [ct.c_void_p, ct.c_float, 
		ct.c_float, ct.c_float, ct.c_char, 
		ct.c_char_p, ct.c_int, ct.POINTER(ct.c_int),
		ct.POINTER(ct.c_double)]

lp_add_cols = lplib.lp_add_cols
lp_add_cols.argtypes = [ 
	ct.c_void_p,                 # lp
	ct.c_int,                    # n cols
	ct.POINTER(ct.c_double), # obj
	ct.POINTER(ct.c_double), # lb
	ct.POINTER(ct.c_double), # ub
	ct.POINTER(ct.c_char),   # integer
	ct.POINTER(ct.POINTER(ct.c_char)) ] # names

lp_add_row = lplib.lp_add_row
lp_add_row.argtypes = [ 
	ct.c_void_p, ct.c_int, ct.POINTER(ct.c_int),
	ct.POINTER(ct.c_double), ct.c_char_p, ct.c_char, 
	ct.c_double ]


lp_x = lplib.lp_x
lp_x.argtypes = [ct.c_void_p]
lp_x.restype = ct.POINTER(ct.c_double)


lp_free = lplib.lp_free
lp_free.argtypes = [ct.c_void_p]

#lp = LinearProgram()
#lp.read("a.lp")
#lp.optimize()


