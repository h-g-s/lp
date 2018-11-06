import ctypes
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
			charInt = ctypes.c_char(1)
		else:
			charInt = ctypes.c_char(0)

		vrowInt = (ctypes.c_int * nz)()
		for i in range(nz):
			vrowInt[i] = rowInt[i]

		vrowCoef = (ctypes.c_double * nz)()
		for i in range(nz):
			vrowCoef[i] = rowCoef[i]

		lp_add_col(self._plp, 
			ctypes.c_double(obj),
			ctypes.c_double(lb),
			ctypes.c_double(ub),
			charInt, name.encode('utf-8'),
			ctypes.c_int(len(rowIdx)), vrowInt, vrowCoef )

	def add_cols(self, 
		n : int, 
		obj : List[float],
		lb : List[float],
		ub : List[float],
		integer : List[bool],
		name : List[str] ):

		vobj = (ctypes.c_double*n)()
		vlb = (ctypes.c_double*n)()
		vub = (ctypes.c_double*n)()
		vint = (ctypes.c_char*n)()

		for i in range(n):
			vobj[i] = obj[i]

		for i in range(n):
			vlb[i] = lb[i]

		for i in range(n):
			vub[i] = ub[i]

		for i in range(n):
			vint[i] = ctypes.c_char( integer[i] )

		pnames = (ctypes.POINTER(ctypes.c_char)*n)()

		#charSpace : int = 0
		#for i in range(n):
	#		charSpace += len(name[i]) + 1


		for i in range(0, n):
			pnames[i] = ctypes.create_string_buffer(name[i].encode('utf-8'))
			

		lp_add_cols(
			self._plp,
			ctypes.c_int(n),
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

		vrowInt = (ctypes.c_int * nz)()
		for i in range(nz):
			vrowInt[i] = idx[i]

		vrowCoef = (ctypes.c_double * nz)()
		for i in range(nz):
			vrowCoef[i] = coef[i]

		if sense == Sense.EQUAL:
			chars = ctypes.c_char(ord('E'))
		elif sense == Sense.GREATER_OR_EQUAL:
			chars = ctypes.c_char(ord('G'))
		elif sense == Sense.LESS_OR_EQUAL:
			chars = ctypes.c_char(ord('L'))
		else:
			raise 'sense not recognized: {}'.format(sense)

		lp_add_row(
			self._plp,
			len(idx), vrowInt, vrowCoef, 
			ctypes.c_char_p(name.encode('utf-8')), chars, ctypes.c_double(rhs) )





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

	def optimize(self, forceContinuous=False) -> int:
		status = OptimizationStatus.Error
		if forceContinuous==False:
			status = lp_optimize(self._plp)
		else:
			status = lp_optimize_as_continuous(self._plp)

		return status


	def obj_value(self) -> float:
		return lp_obj_value(self._plp)


	def x(self) -> List[float] :
		n = self.cols()
		res = [0.0]*n
		pr = lp_x(self._plp)
		for i in range(0,n):
			res[i] = pr[i]

		return res


	def __del__(self):
		lp_free(self._plp)

lplib = ctypes.CDLL('./lp-cbc-linux64.so')

lp_create = lplib.lp_create
lp_create.restype = ctypes.c_void_p

lp_read = lplib.lp_read
lp_read.argtypes = [ctypes.c_void_p, ctypes.c_char_p]

lp_write = lplib.lp_write_lp
lp_write.argtypes = [ctypes.c_void_p, ctypes.c_char_p]

lp_optimize = lplib.lp_optimize
lp_optimize.argtypes = [ctypes.c_void_p]
lp_optimize.restype = ctypes.c_int

lp_optimize_as_continuous = lplib.lp_optimize_as_continuous
lp_optimize_as_continuous.argtypes = [ctypes.c_void_p]
lp_optimize_as_continuous.restype = ctypes.c_int

lp_cols = lplib.lp_cols
lp_cols.argtypes = [ctypes.c_void_p]
lp_cols.restype = ctypes.c_int

lp_rows = lplib.lp_rows
lp_rows.argtypes = [ctypes.c_void_p]
lp_rows.restype = ctypes.c_int

lp_obj_value = lplib.lp_obj_value
lp_obj_value.argtypes = [ctypes.c_void_p]
lp_obj_value.restype = ctypes.c_float

lp_add_col = lplib.lp_add_col
lp_add_col.argtypes = [ctypes.c_void_p, ctypes.c_float, 
		ctypes.c_float, ctypes.c_float, ctypes.c_char, 
		ctypes.c_char_p, ctypes.c_int, ctypes.POINTER(ctypes.c_int),
		ctypes.POINTER(ctypes.c_double)]

lp_add_cols = lplib.lp_add_cols
lp_add_cols.argtypes = [ 
	ctypes.c_void_p,                 # lp
	ctypes.c_int,                    # n cols
	ctypes.POINTER(ctypes.c_double), # obj
	ctypes.POINTER(ctypes.c_double), # lb
	ctypes.POINTER(ctypes.c_double), # ub
	ctypes.POINTER(ctypes.c_char),   # integer
	ctypes.POINTER(ctypes.POINTER(ctypes.c_char)) ] # names

lp_add_row = lplib.lp_add_row
lp_add_row.argtypes = [ 
	ctypes.c_void_p, ctypes.c_int, ctypes.POINTER(ctypes.c_int),
	ctypes.POINTER(ctypes.c_double), ctypes.c_char_p, ctypes.c_char, 
	ctypes.c_double ]


lp_x = lplib.lp_x
lp_x.argtypes = [ctypes.c_void_p]
lp_x.restype = ctypes.POINTER(ctypes.c_double)


lp_free = lplib.lp_free
lp_free.argtypes = [ctypes.c_void_p]

#lp = LinearProgram()
#lp.read("a.lp")
#lp.optimize()


