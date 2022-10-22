#!/usr/bin/env python3

from functools import reduce

from math import floor, ceil

from random import random

from .ctx import mst, Rational

def Frac_div(n, d, mst=mst) :
	mst = mst.accuracy_protect
	if mst is None :
		return n / d
	try :
		return n / mst(d)
	except :
		return mst(n) / d

def _gcd(a, b) :
	while b :
		a, b = b, a%b
	return a

def gcd(*n) :
	return reduce(_gcd, n)

def lcm(*n) :
	return reduce(lambda a, b, _div=Frac_div : _div(a*b, _gcd(a, b)), n)
	
def ratioSimp(n, use_lcm=False) :
	'''化简比例
	: n, list: 待化简比例
	: use_lcm, bool: 是否使用最小公倍数法化简
	return list
	'''
	if use_lcm :
		lcm_of_dens = lcm(*[Rational(str(i)).denominator for i in n])
		return [int(i*lcm_of_dens) for i in n]
		
	gcd_of_n = reduce(_gcd, n)
	return [int(i//gcd_of_n) for i in n]

def mayint(n) :
	if type(n) is int :
		return True
	if hasattr(n, 'denominator') :
		return n.denominator == 1
	if isinstance(n, Mat) :
		return all(map(mayint, n))
	if type(n) is dict :
		return all(map(mayint, n.values()))
	try :
		return int(n) == n
	except :
		return False
		
def zeros(*shape) :
	return Mat([[0]*shape[-1] for _ in range(shape[0])])
	
def ones(*shape) :
	return Mat([[1]*shape[-1] for _ in range(shape[0])])

def diag(elements) :
	if type(elements) in (tuple, list) :
		_el = elements 
	else :
		_el = tuple(elements)
		
	l = len(_el)
	new = [[0]*l for _ in range(l)]
	
	for piv, val in enumerate(_el) :
		new[piv][piv] = val
		
	return Mat(new)

def eye(*shape) :
	cols = shape[-1]
	
	new = [[0]*cols for _ in range(shape[0])]
	
	for piv, row in enumerate(new) :
		row[piv] = 1
		
	return Mat(new)

def randmat(*shape, min=0, max=10, type=None,
		## HACK: hand-optimized; turn globals into locals
		random=random,
		floor=floor,
		range=range,
		int=int,
		float=float,
		) :
	cols, m, n = \
		shape[-1], max + 1 - min, max - min
	uniform, randint = \
		lambda : n*random() + min, \
		lambda : floor(m * random() + min)
	
	if type is None :
		return Mat([
			[
				(randint if floor(random()+0.5) else uniform)() for j in range(cols)
			]
			for i in range(shape[0])
		])

	randnum = randint if type is int \
		else  uniform if type is float \
		else  lambda type=type, f='{:.3f}'.format: type(f(uniform()))
		
	return Mat([
		[
			randnum() for j in range(cols)
		]
		for i in range(shape[0])
	])
		
def hilbert(*shape, _div=Frac_div) :
	return Mat([
		[
			_div(1, i + j)
			for j in range(shape[-1])
		]
		for i in range(1, shape[0]+1)
	])


class Mat(object) :
	'''A simple matrix class.
	
	Examples
	==========
	
	>>> from fmat import *
	
	#--- statement initialization matrix ---
	
	>>> eye(5)
	Mat(
	┌ 1  0  0  0  0 ┐
	│ 0  1  0  0  0 │
	│ 0  0  1  0  0 │
	│ 0  0  0  1  0 │
	└ 0  0  0  0  1 ┘)
	
	Also zeros、ones、diag、randmat
	 and hilbert
	
	>>> a = Mat([[1,2,2,2], [2,4,6,8], [3,6,8,10]])
	Mat(
	┌ 1  2  2   2 ┐
	│ 2  4  6   8 │
	└ 3  6  8  10 ┘)
	
	#--- slice ---
	
	Index start from 0.
	
	>>> a[1, 2] # get value by slice, get row 1 col 2
	6
	>>> a[1:, 2]
	Mat(
	┌ 6 ┐
	└ 8 ┘)
	
	>>> a[1, 2] = 1000 # set value by slice
	>>> a
	Mat(
	┌ 1  2	 2   2 ┐
	│ 2  4  1000   8 │
	└ 3  6	 8  10 ┘)
	
	>>> a[1:, 2] = [[60], [80]] # set more values by slice
	>>> a
	Mat(
	┌ 1  2   2   2 ┐
	│ 2  4  60   8 │
	└ 3  6  80  10 ┘)
	
	#--- get matrix shape ---
	
	>>> a.shape
	(3, 4)
	
	>>> a.rows
	3
	>>> a.cols
	4
	
	#--- insert at the end of the matrix, and only ---
	
	>>> a.colinsert(eye(3))
	>>> a
	Mat(
	┌ 1  2  2   2  1  0  0 ┐
	│ 2  4  6   8  0  1  0 │
	└ 3  6  8  10  0  0  1 ┘)
	
	>>> del a[:, 4:]
	>>> a
	Mat(
	┌ 1  2  2   2 ┐
	│ 2  4  6   8 │
	└ 3  6  8  10 ┘)
	
	>>> a.rowinsert(eye(4))
	>>> a
	Mat(
	┌ 1  2  2   2 ┐
	│ 2  4  6   8 │
	│ 3  6  8  10 │
	│ 1  0  0   0 │
	│ 0  1  0   0 │
	│ 0  0  1   0 │
	└ 0  0  0   1 ┘)

	>>> del a[3:, :]
	>>> a
	Mat(
	┌ 1  2  2   2 ┐
	│ 2  4  6   8 │
	└ 3  6  8  10 ┘)
	
	#--- linear algebra ---
	
	 RREF
	  return tuple(rref of matrix, pivot_cols)
	>>> a.rref()
	(Mat(
	┌ 1  2  0  -2 ┐
	│ 0  0  1   2 │
	└ 0  0  0   0 ┘), (0, 2))

	 determinant
	>>> a.det()
	0
	
	 null space
	>>> a.nullspace()
	[Mat(
	┌ -2 ┐
	│  1 │
	│  0 │
	└  0 ┘), Mat(
	┌  2 ┐
	│  0 │
	│ -2 │
	└  1 ┘)]

	 column space
	>>> a.columnspace()
	[Mat(
	┌ 1 ┐
	│ 0 │
	└ 0 ┘), Mat(
	┌ 0 ┐
	│ 1 │
	└ 0 ┘)]

	 rank
	>>> a.rank()
	2
	
	 transpose
	>>> a.T
	┌ 1  2   3 ┐
	│ 2  4   6 │
	│ 2  6   8 │
	└ 2  8  10 ┘

	 conjugate transpose
	>>> a.H
	┌ 1  2   3 ┐
	│ 2  4   6 │
	│ 2  6   8 │
	└ 2  8  10 ┘
	
	and many other mathods...

	==========
	'''
	__slots__ = (
		'__data',
		'__shape',
		'__rows',
		'__cols',
	)
		
	def __init__(self, A=0, type=type, len=len) :
		'''Initialization matrix class.
		
		Examples
		
		>>> Mat().tolist()
		[]
		
		>>> Mat(3).tolist()
		[[], [], []]
		
		>>> Mat([1, 2, 3]).tolist()
		[[1 2 3]]
		
		>>> Mat([[1, 2], [3, 4]])
		Mat(
		┌ 1 2 ┐
		└ 1 2 ┘)
		
		>>> Mat([(1, 2), (3, 4)])
		Mat(
		┌ 1 2 ┐
		└ 1 2 ┘)
		
		>>> Mat(((1, 2), (3, 4)))
		Mat(
		┌ 1 2 ┐
		└ 1 2 ┘)
		
			Do not do like this, it will cause bugs
			
		>>> Mat(([1, 2], [3, 4]))
		Mat(
		┌ 1 2 ┐
		└ 1 2 ┘)
		
		>>> Mat(([1, 2], [3, 4]))._Mat__data # this is matrix's internal implementation
		([1, 2], [3, 4])
		'''
		if type(A) is int :
			self.__data, self.__shape, self.__rows, self.__cols = \
				[[] for _ in range(A)], (A, 0), A, 0
			return
		if type(A[0]) is list :
			self.__data, self.__shape = \
				A, (len(A), len(A[0]))
		elif type(A[0]) is tuple :
			self.__data, self.__shape = \
				[list(i) for i in A], (len(A), len(A[0]))
		else :
			self.__data, self.__shape = \
				[A], (1, len(A))
		self.__rows, self.__cols = self.__shape
			
	def __add__(self, other) :
		if isinstance(other, Mat) :
			return Mat([
				[k + u for k, u in zip(i, j)]
				for i, j in zip(self.__data, other.__data)
			])

		return Mat([
			[other + j for j in i]
			for i in self.__data
		])
		
	def __radd__(self, other) :
		return Mat([
			[other + j for j in i]
			for i in self.__data
		])
	
	def __sub__(self, other) :
		if isinstance(other, Mat) :
			return Mat([
				[k - u for k, u in zip(i, j)]
				for i, j in zip(self.__data, other.__data)
			])
		
		return Mat([
			[j - other for j in i]
			for i in self.__data
		])
			
	def __rsub__(self, other) :
		return Mat([
			[other - j for j in i]
			for i in self.__data
		])
	
	def __mul__(self, other, sum=sum, zip=zip, range=range, isinstance=isinstance) :
		if isinstance(other, Mat) :
			# TODO : Use Coppersmith-Winograd Algorithm?
			return Mat([
				[
					sum(k * u[j] for k, u in zip(i, other.__data))
					for j in range(other.__cols)
				]
				for i in self.__data
			])
	
		return Mat([
			[other * j for j in i]
			for i in self.__data
		])
		
	def __rmul__(self, other) :
		return Mat([
			[other * j for j in i]
			for i in self.__data
		])
	
	__matmul__ = __mul__
	__rmatmul__ = __rmul__
	
	def hada(self, other, zip=zip) :
		'''To compute the Hadamard product of two matrixes.'''
		return Mat([
			[k * u for k, u in zip(i, j)]
			for i, j in zip(self.__data, other.__data)
		])

	def kron(self, other) :
		'''To compute the Kronecker Product of two matrixes(A ⊗ B).
		Example
		>>> Mat([[1, 2], [3, 1]]).kron(Mat([[0, 3], [2, 1]]))
		Mat(
		┌ 0  3  0  6 ┐
		│ 2  1  4  2 │
		│ 0  9  0  3 │
		└ 6  3  2  1 ┘)
		'''
		n = self.__cols*other.__cols
		om, on = other.__shape
		new = []
		for i in range(self.__rows*other.__rows) :
			sdido, odimo = self.__data[i // om], other.__data[i % om]
			new.append([sdido[j // on] * odimo[j % on] for j in range(n)])
		return Mat(new)
	
	def __div__(self, other, zip=zip, map=map, _div=Frac_div) :
		if isinstance(other, Mat) :
			return Mat([
				list(map(_div, i, j))
				# [_div(k, u) for k, u in zip(i, j)]
				for i, j in zip(self.__data, other.__data)
			])

		return Mat([
			[_div(j, other) for j in i]
			for i in self.__data
		])
		
	__floordiv__ = __truediv__ = __div__
	
	def __rdiv__(self, other) :
		return Mat([
			[Frac_div(other, j) for j in i]
			for i in self.__data
		])
		
	__rfloordiv__ = __rtruediv__ = __rdiv__
	
	def __pow__(self, other) :
		if not other :
			return eye(self.__rows)
		
		j, k = 1, self.copy()
		if other < 0 :
			other = -other
			while other :
				if other & 1 : j *= k
				k, other = k*k, other>>1
			return j._inverse()
		
		while other :
			if other & 1 : j *= k
			k, other = k*k, other>>1
		return j
		
	def _inverse(self) :
		'''To compute the inverse of matrix.'''
		reduced, pivots = self.copy().colinsert(eye(self.__rows)).rref()
		del reduced[:, :self.__rows]
		return reduced
		
	def __neg__(self) :
		return Mat([
			[-j for j in i] for i in self.__data
		])
		
	def __getitem__(self, key) :
		if type(key) is not tuple :
			if self.__rows == 1 :
				if type(key) is slice :
					return Mat(self.__data[0][key])
				return self.__data[0][key]
				
			if self.__cols == 1 :
				if type(key) is slice :
					return Mat([
						[i[0]]
						for i in self.__data[key]
					])
				return self.__data[key][0]
				
			return self.__data[key]
		
		row, col = key
		if type(row) is type(col) is slice :
			return Mat([
				i[col] for i in self.__data[row]
			])
		if type(row) is slice :
			return Mat([
				[i[col]] for i in self.__data[row]
			])
		if type(col) is slice :
			return Mat(self.__data[row][col])
			
		return self.__data[row][col]
	
	def __setitem__(self, key, value) :
		if type(key) is not tuple :
			if self.__rows == 1 :
				key = (0, key)
			elif self.__cols == 1 :
				key = (key, 0)
			
		if isinstance(value, Mat) :
			value = value.__data
			
		if type(key) is tuple :
			row, col = key
			if type(row) is type(col) is slice :
				for i, j in zip(self.__data[row], value) :
					i[col] = j
			elif type(row) is slice :
				try :
					for i, j in zip(self.__data[row], value) :
						i[col] = j[0]
				except :
					# about this, to make a[:2, 2] = [1, 2]
					# could work as a[:2, 2] = [[1], [2]]
					for i, j in zip(self.__data[row], value) :
						i[col] = j
			else :
				self.__data[row][col] = value
		else :
			self.__data[key] = value
			
	def __delitem__(self, key) :
		if type(key) is not tuple :
			if self.__rows == 1 :
				key = (0, key)
			elif self.__cols == 1 :
				key = (key, 0)

		if type(key) is tuple :
			row, col = key
			if type(row) is type(col) is slice :
				for i in self.__data[row] :
					del i[col]
				self.__data = [i for i in self.__data if i]
			elif type(row) is slice :
				for i in self.__data[row] :
					del i[col]
			else :
				del self.__data[row][col]
		else :
			del self.__data[key]
		self.__shape = Mat(self.__data).__shape
		self.__rows, self.__cols = self.__shape
		
	def __eq__(self, other) :
		return self.__shape == other.__shape \
				and self.__data == other.__data
	
	def __len__(self) :
		if self.__rows == 1 :
			return self.__cols
		return self.__rows # do it like numpy

	@property
	def T(self) :
		return Mat(list(map(list, zip(*self.__data))))
		
	@property
	def H(self) :
		return Mat([[j.conjugate() for j in i] for i in zip(*self.__data)])

	@property
	def C(self) :
		return Mat([[j.conjugate() for j in i] for i in self.__data])
	
	def tolist(self) :
		return [i[:] for i in self.__data]
		
	def totype(self, type, map=map) :
		'''To format by `type` for each element of the matrix in place.
		Notice : Returns the modified matrix.
		
		Examples
		>>> from fmat import *
		
		>>> randmat(3)
		Mat(
		┌ 0.27956954001250645  2.3958028902679063                  0 ┐
		│                   2                  10  4.283462965997623 │
		└   3.421341518557158   7.048751872614475                  5 ┘)
		>>> _.totype(int)
		Mat(
		┌ 0   2  0 ┐
		│ 2  10  4 │
		└ 3   7  5 ┘)
		>>> _
		Mat(
		┌ 0   2  0 ┐
		│ 2  10  4 │
		└ 3   7  5 ┘)
		'''
		for i in self.__data : i[:] = map(type, i)
		return self
	
	def __iter__(self) :
		for i in self.__data :
			for j in i :
				yield j # do it like mpmath
	
	def __str__(self, str=str) :
		res = []
		maxlen = [0] * self.__cols
		for i in self.__data :
			res.append([])
			for j, k in enumerate(map(str, i)) :
				res[-1].append(k)
				maxlen[j] = max(len(k), maxlen[j])
		
		alignment = lambda row, ard='[]', sep='  ', nnjoin=' '.join, align=mst.alignment, maxlen=maxlen : \
			nnjoin((ard[0], sep.join(map(align, row, maxlen)), ard[-1]))
		# ┌ {} ┐	⎡ {} ⎤
		# │ {} │	⎢ {} ⎥
		# └ {} ┘	⎣ {} ⎦
		if mst.pretty_print :
			if self.__rows < 2 :
				return alignment(res.pop())
				
			res[0] = alignment(res[0], '┌┐')
			if self.__rows > 2 :
				res[1:-1] = [
					alignment(row, '││')
					for row in res[1:-1]
				]
			res[-1] = alignment(res[-1], '└┘')
			return '\n'.join(res)
		
		return ',\n '.join(map(lambda i : alignment(i, sep=', '), res)).join('[]')
		
	def __repr__(self) :
		return ''.join((
			self.__class__.__name__, '(', '\n'*int(self.__rows>1), str(self), ')'))
	
	def branch_cut(C, _A, _b, maximize=False, lp_solve=None) :
		'''! To Be Implemented !
		
		This method is used to solve ILP problem by branch and cut algorithm.
		It solve some :
			min C.T*X
			s.t. AX <= b.T
			x >= 0, x ∈ Integers
			
		Example
		<
			max 2x + y
			s.t. 
					 5y <= 15
				6x + 2y <= 24
				 x +  y <= 5
			x,y >= 0, x,y ∈ Integers
		>
		>>> C = Mat([2, 1])
		>>> A = Mat([[0, 5], [6, 2], [1, 1]])
		>>> b = Mat([15, 24, 5])
		
		>>> C.branch_cut(A, b, maximize=True)
		(Rational(8, 1), {1: Rational(4, 1), 2: Rational(0, 1)})
		'''
		if lp_solve is None : lp_solve = Mat.relax_simplex
		
		bestvalue = -0x7fffffff if maximize else 0x7fffffff
		nodes = [(_A.copy(), _b.copy())]
		optimal_solution = None
	
		while nodes :
			A, b = nodes.pop(0)
			while True : # python has no do while...unhappy...
				node = lp_solve(C, A, b, maximize)
				#print(node)
				cutting_planes_found = False
				
				if not node : break
				val, vars, stable = node
				if val <= bestvalue if maximize else val >= bestvalue :
					break
				if mayint(vars) :
					bestvalue, optimal_solution, _ = node
					break
	
				for i in range(1, C.__cols+1) :
					idn = vars.get(i, 0)
					if not mayint(idn) :
						idn = i - 1 # return to zero start
						break
						
				xi, bi = vars[idn + 1], b[idn]
				cutting_planes_A, cutting_planes_b = [0]*A.__cols, []
				j = A.__data[idn][idn]*xi
				cutting_planes_A[idn] = floor(j) - j
						
				cutting_planes_b.append(floor(bi) - bi)
				
				#print(cutting_planes_A, cutting_planes_b)
				if 0and cutting_planes_A :
					A.rowinsert(Mat(cutting_planes_A))
					b.colinsert(Mat(cutting_planes_b))
					cutting_planes_found = True
					
					continue
				
				nodes.extend(((
					Mat(A.copy().__data + [[int(i==idn) for i in range(A.__cols)]]),
					Mat(b.__data[0] + [floor(node[1][idn + 1])])
				), (
					Mat(A.copy().__data + [[-int(i==idn) for i in range(A.__cols)]]),
					Mat(b.__data[0] + [-ceil(node[1][idn + 1])])
				)))
				
				if not cutting_planes_found :
					break
					
		return bestvalue, optimal_solution

	def branch_bound(C, _A, _b, maximize=False, lp_solve=None) :
		'''This method is used to solve ILP problem by branch and bound algorithm.
		It solve some :
			min C.T*X
			s.t. AX <= b.T
			x >= 0, x ∈ Integers
			
		Example
		<
			max 2x + y
			s.t. 
					 5y <= 15
				6x + 2y <= 24
				 x +  y <= 5
			x,y >= 0, x,y ∈ Integers
		>
		>>> C = Mat([2, 1])
		>>> A = Mat([[0, 5], [6, 2], [1, 1]])
		>>> b = Mat([15, 24, 5])
		
		>>> C.branch_bound(A, b, maximize=True)
		(Rational(8, 1), {1: Rational(4, 1), 2: Rational(0, 1)})
		'''
		if lp_solve is None : lp_solve = Mat.simplex
		
		bestvalue = -0x7fffffff if maximize else 0x7fffffff
		nodes = [(_A.copy(), _b.copy())]
		optimal_solution = None
	
		while nodes :
			A, b = nodes.pop(0)
			node = lp_solve(C, A, b, maximize)
			
			if not node :
				continue
			if node[0] <= bestvalue if maximize else node[0] >= bestvalue :
				continue
			if mayint(node[1]) :
				bestvalue, optimal_solution = node
				continue
			
			for i in range(1, C.__cols+1) :
				idn = node[1].get(i, 0)
				if not mayint(idn) :
					idn = i - 1 # return to zero start
					break
			
			nodes.extend(((
				Mat(A.copy().__data + [[int(i==idn) for i in range(A.__cols)]]),
				Mat(b.__data[0] + [floor(node[1][idn + 1])])
			), (
				Mat(A.copy().__data + [[-int(i==idn) for i in range(A.__cols)]]),
				Mat(b.__data[0] + [-ceil(node[1][idn + 1])])
			)))
			
		return bestvalue, optimal_solution
	
	def simplex(C, A, b, maximize=False,
			Frac_div=Frac_div,
			enumerate=enumerate,
			min=min, zip=zip,
			range=range,
			list=list
		) :
		'''This method is used to solve LP problem by simplex algorithm.
		It solve some :
			min C.T*X
			s.t. AX <= b.T
			x >= 0
			
		Example
		<
			max 2x + y
			s.t. 
					 5y <= 15
				6x + 2y <= 24
				 x +  y <= 5
			x,y >= 0
		>
		>>> C = Mat([2, 1])
		>>> A = Mat([[0, 5], [6, 2], [1, 1]])
		>>> b = Mat([15, 24, 5])
		
		>>> C.simplex(A, b, maximize=True)
		(Rational(17, 2), {1: Rational(7, 2), 2: Rational(3, 2)})
		
		References :
			1.https://www.jianshu.com/p/dd0761a2fdfd
			2.https://www.hrwhisper.me/introduction-to-simplex-algorithm/
		'''
		M = Mat([0]).colinsert(-C if maximize else C)
		M.rowinsert(b.T.colinsert(A))
		
		def _pivot(mat, B, row, col) :
			dvb = mat.__data[row][col]
			mat.__data[row] = [Frac_div(i, dvb) for i in mat.__data[row]]
			for j in mat.__data[:row]+mat.__data[row+1:] :
				ji = j[col]
				for k,u in enumerate(mat.__data[row]) :
					j[k] = j[k] - u*ji
			B[row] = col
		
		def _simplex(mat, B, m, n) :
			while min(mat.__data[0][1:]) < 0 :
				for i,j in enumerate(mat.__data[0][1:], 1) :
					if j < 0 :
						col, j = (
							i, # use Bland's method to avoid degeneracy. use mat[0].argget(min) ok?
							0x7fffffff # temporary, the min value
						)
						for k,u in enumerate(mat.__data[1:], 1) :
							if u[col] > 0 :
								i = u[0] / u[col] # temporary min value
								if i < j :
									row, j = k, i # find the theta index
							else :
								continue
						break
				try :
					if mat.__data[row][col] <= 0 : return None  # the theta is ∞, the problem is unbounded
				except :
					# Cause unknown reasons, `row` won't be defined, but question is solvable
					break
				_pivot(mat, B, row, col)
			return (mat.__data[0][0] if maximize else -mat.__data[0][0],
				{i: j[0] for i,j in zip(B[1:m], mat.__data[1:m])}) # if i < n})
	
		m, n = M.__shape  # m - 1 is the number slack variables we should add
		B = list(range(n - 1, n + m - 1))  # add diagonal array
		M.colinsert(Mat([[0]*(m-1)] + eye(m-1).__data))  # combine them!
		if min(i[0] for i in M.__data[1:]) < 0 :  # is the initial basic solution feasible?
			minval = 0x7fffffff # temporary
			for i,j in enumerate(M.__data[1:], 1) :
				if j[0] < minval :
					row, minval = i, j[0] # find the index of min b
			temp, M.__data[0] = Mat([M.__data[0].copy()]), [0]*M.__cols  # set first row value to zero, and store the previous value
			M.colinsert(Mat([[1]] + [[-1] for _ in range(m - 1)]))
			_pivot(M, B, row, M.__cols - 1)
			if _simplex(M, B, m, n)[0] : return None  # if it is not equal to 0 that the problem has no answer
			
			if M.__cols - 1 in B :  # if the x0 in B, we should pivot it.
				_pivot(M, B, B.index(M.__cols - 1), [i for i,j in enumerate(M.__data[0][1:]) if j][0]+1)
				
			M = temp.rowinsert(M[1:, :-1])  # recover the first line
			for i, x in enumerate(B[1:], 1) :
				M0x = M.__data[0][x]
				M.__data[0] = [j - k*M0x for j,k in zip(M.__data[0], M.__data[i])]
		return _simplex(M, B, m, n)
		
	def relax_simplex(C, A, b, maximize=False) :
		'''This method is used to solve LP problem by simplex algorithm.
		And realize this just for the generate the cutting plane in
			the Gomory's cutting on integers linear programming.
		It will return best value, bases vars, the simplex table.
		And it's enter also can just not to be the relax form.
		Had better not call this method outside.
		
		Notice :
			return (best value, bases vars, simplex table)
		And! The simplex table is a kind of form like :
			┌ best_value   non_bases_vars ┐
			│							 │
			└ bases_vars the_coefficients ┘
		
		It solve some :
			min C.T*X
			s.t. AX <= b.T
			x >= 0
			
		Example
		<
			max 2x + y
			s.t. 
					 5y <= 15
				6x + 2y <= 24
				 x +  y <= 5
			x,y >= 0
		>
		>>> C = Mat([2, 1])
		>>> A = Mat([[0, 5], [6, 2], [1, 1]])
		>>> b = Mat([15, 24, 5])
		
		>>> C.relax_simplex(A, b, maximize=True)
		(Rational(17, 2), {1: Rational(7, 2), 2: Rational(3, 2)})
		
		References :
			1.https://www.jianshu.com/p/dd0761a2fdfd
			2.https://www.hrwhisper.me/introduction-to-simplex-algorithm/
		'''
		M = Mat([0]).colinsert(-C if maximize else C)
		M.rowinsert(b.T.colinsert(A))
		
		def _pivot(mat, B, row, col) :
			dvb = mat.__data[row][col]
			mat.__data[row] = [Frac_div(i, dvb) for i in mat.__data[row]]
			for j in mat.__data[:row]+mat.__data[row+1:] :
				ji = j[col]
				for k,u in enumerate(mat.__data[row]) :
					j[k] = j[k] - u*ji
			B[row] = col
		
		def _simplex(mat, B, m, n) :
			while min(mat.__data[0][1:]) < 0 :
				for i,j in enumerate(mat.__data[0][1:], 1) :
					if j < 0 :
						col, j = (
							i, # use Bland's method to avoid degeneracy. use mat[0].argget(min) ok?
							0x7fffffff # temporary, the min value
						)
						for k,u in enumerate(mat.__data[1:], 1) :
							if u[col] > 0 :
								i = u[0] / u[col] # temporary min value
								if i < j :
									row, j = k, i # find the theta index
							else :
								continue
						break
				
				try :
					if mat.__data[row][col] <= 0 : return None  # the theta is ∞, the problem is unbounded
				except :
					# Cause unknown reasons, `row` won't be defined, but question is solvable
					break
				_pivot(mat, B, row, col)
			return (mat.__data[0][0] if maximize else -mat.__data[0][0],
				{i: j[0] for i,j in zip(B[1:m], mat.__data[1:m])}, mat) # if i < n})
	
		m, n = M.__shape  # m - 1 is the number slack variables we should add
		B = list(range(n - 1, n + m - 1))  # add diagonal array
		M.colinsert(Mat([[0]*(m-1)] + eye(m-1).__data))  # combine them!
		if min(i[0] for i in M.__data[1:]) < 0 :  # is the initial basic solution feasible?
			minval = 0x7fffffff # temporary
			for i,j in enumerate(M.__data[1:], 1) :
				if j[0] < minval :
					row, minval = i, j[0] # find the index of min b
			temp, M.__data[0] = Mat([M.__data[0].copy()]), [0]*M.__cols  # set first row value to zero, and store the previous value
			M.colinsert(Mat([[1]] + [[-1] for _ in range(m - 1)]))
			_pivot(M, B, row, M.__cols - 1)
			if _simplex(M, B, m, n)[0] : return None  # if it is not equal to 0 that the problem has no answer
			
			if M.__cols - 1 in B :  # if the x0 in B, we should pivot it.
				_pivot(M, B, B.index(M.__cols - 1), [i for i,j in enumerate(M.__data[0][1:]) if j][0]+1)
				
			M = temp.rowinsert(M[1:, :-1])  # recover the first line
			for i, x in enumerate(B[1:], 1) :
				M0x = M.__data[0][x]
				M.__data[0] = [j - k*M0x for j,k in zip(M.__data[0], M.__data[i])]
		return _simplex(M, B, m, n)
		
	def nullspace(self, simp=False) :
		reduced, pivots = self.rref()
		
		free_vars = [i for i in range(self.__cols) if i not in pivots]
		basis = []
		
		for free_var in free_vars :
			# For each free variable, we will set it to 1 and others to 0
			# Then, we will use back substitutppion to solve the system
			vec 		  = [0] * self.__cols
			vec[free_var] = 1
	
			for row, piv_col in zip(reduced.__data, pivots) :
				vec[piv_col] -= row[free_var]
		
			basis.append(ratioSimp(vec) if simp else vec)

		return [Mat(b).T for b in basis]
		
	def columnspace(self, simp=False) :
		reduced, pivots = self.rref()
		basis = list(map(reduced.col, pivots))
		
		if simp :
			basis[:] = map(lambda i: Mat(i).T, map(ratioSimp, basis))
		
		return basis

	def _rref(self, abs=abs, enumerate=enumerate) :
		'''To get reduced row echelon form by Gaussian-Jordan Elimination'''
		new_mat, fi, pivots = \
			[i[:] for i in self.__data], 0, []
			
		for i in range(self.__cols) :
			try :
				# 尝试寻找第i列中第fi行到末行里绝对值最大的系数作为fi行主元
				colMax = maxRow = 0
				for j,k in enumerate(new_mat[fi:], fi) :
					ki = k[i]
					maxRow, colMax = (maxRow, colMax) if abs(colMax) > abs(ki) else (j, ki)
			except :
				break # 嗯...能跑就不要动...???
				
			if abs(colMax) <= mst.eps : continue # 主元为0是自由列,直接跳过
			pivots.append(i) # 主元列计数
			
			if maxRow != fi : # 交换行同时将主元化为一
				new_mat[fi], new_mat[maxRow] = \
					[Frac_div(j, colMax) for j in new_mat[maxRow]], new_mat[fi]
			else : # 否则直接将主元化为一
				new_mat[fi] = [Frac_div(j, colMax) for j in new_mat[fi]]
			
			for j in new_mat[:fi]+new_mat[fi+1:] :
				ji = j[i]
				for k,u in enumerate(new_mat[fi]) :
					j[k] = j[k] - u*ji
			
			fi += 1
			
		return Mat(new_mat), tuple(pivots) # do it like sympy
	
	def rref(self, abs=abs, map=map, list=list, enumerate=enumerate) :
		'''To get reduced row echelon form by Gaussian-Jordan Elimination'''
		new_mat, fi, pivots = \
			[i[:] for i in self.__data], 0, []
			
		try :
			for i in range(self.__cols) :
				# 尝试寻找第i列中第fi行到末行里绝对值最大的系数作为fi行主元
				colMax = maxRow = 0
				for j,k in enumerate(new_mat[fi:], fi) :
					ki = k[i]
					maxRow, colMax = (maxRow, colMax) if abs(colMax) > abs(ki) else (j, ki)
			
				if abs(colMax) <= mst.eps : continue # 主元为0是自由列,直接跳过
				pivots.append(i) # 主元列计数
				
				_div = lambda j, colMax=colMax, _d=Frac_div: \
					_d(j, colMax)
				if maxRow != fi : # 交换行同时将主元化为一
					new_mat[maxRow], new_mat[fi] = \
						new_mat[fi], list(map(_div, new_mat[maxRow]))
				else : # 否则直接将主元化为一
					new_mat[fi][:] = map(_div, new_mat[fi])
				
				for j in new_mat[:fi]+new_mat[fi+1:] :
					j[:] = map(lambda k,u,ji=j[i] : k-u*ji, j, new_mat[fi])
				
				fi += 1
		finally :
			return Mat(new_mat), tuple(pivots) # do it like sympy
	
	def rank(self, abs=abs, enumerate=enumerate) :
		'''To get rank by Gaussian Elimination'''
		new_mat, rk = \
			[i[:] for i in self.__data], 0
			
		for i in range(self.__cols) :
			try :
				# 尝试寻找第i列中第rk行到末行里绝对值最大的系数作为rk行主元
				colMax = maxRow = 0
				for j,k in enumerate(new_mat[rk:], rk) :
					ki = k[i]
					maxRow, colMax = (maxRow, colMax) if abs(colMax) > abs(ki) else (j, ki)
			except :
				break # 嗯...能跑就不要动...???
				
			if abs(colMax) <= mst.eps : continue # 主元为0是自由列,直接跳过
			
			if maxRow != rk : # 交换行同时将主元化为一
				new_mat[rk], new_mat[maxRow] = \
					[Frac_div(j, colMax) for j in new_mat[maxRow]], new_mat[rk]
			else : # 否则直接将主元化为一
				new_mat[rk] = [Frac_div(j, colMax) for j in new_mat[rk]]
			
			for j in new_mat[rk+1:] :
				ji = j[i]
				for k,u in enumerate(new_mat[rk]) :
					j[k] = j[k] - u*ji
			rk += 1 # rank count
			
		return rk

	def det(self, abs=abs) :
		'''To get determinant by Bareiss algorithm
		Example
		>>> from fmat import *
		
		>>> mst.accuracy_protect = None # don't use full accuracy to compute
		
		>>> Mat(
		... [[ 67,  68,  -3,  71],
		...  [ 75,  27, -21,  30],
		...  [104, -34,  46, 163],
		...  [165, 110, 144,  25]])
		Mat(
		┌  67   68   -3   71 ┐
		│  75   27  -21   30 │
		│ 104  -34   46  163 │
		└ 165  110  144   25 ┘)
		
		>>> _.det()
		139971819.0
		
		Reference :
			1.https://en.wikipedia.org/wiki/Bareiss%20algorithm
		'''
		if self.__rows != self.__cols :
			return 0
			
		new_mat, sym = \
			[i[:] for i in self.__data], 1
		
		for kp,k in enumerate(new_mat[:-1]) :
			piv = k[kp]
			if abs(piv) <= mst.eps :
				for i in new_mat[kp:] :
					if abs(i[kp]) > mst.eps :
						k[:], i[:], sym = \
							i[:], k[:], -sym
						break
				piv = k[kp]
				if abs(piv) <= mst.eps :
					return 0
			for i in new_mat[kp+1:] :
				for j in range(kp+1, self.__cols) :
					# The main of Bareiss algorithm
					i[j] = Frac_div(i[j]*piv - i[kp]*k[j], sym)
			sym = piv
		return new_mat[-1][-1]
		
	def argget(self, getting) :
		gets = getting(self)
		for i,j in enumerate(self.__data) :
			for k,u in enumerate(j) :
				if u == gets :
					return k + i*self.__cols
		
	def where(self, condition) :
		return Mat([
			[
				j for j, k in enumerate(i)
				if condition(k)
			]
			for i in self.__data
		])
		
	def row(self, key) :
		return Mat([self.__data[key]])
		
	def col(self, key) :
		return Mat([[i[key]] for i in self.__data])
			
	def rowswap(self, i, j) :
		self.__data[i], self.__data[j] = \
			self.__data[j], self.__data[i]
		return self
			
	def colswap(self, i, j) :
		for _i in self.__data :
			_i[i], _i[j] = _i[j], _i[i]
		return self
			
	def rowinsert(self, other) :
		'''To insert data in place at the end of each column.
		Notice : Returns the inserted inplacely matrix.
		
		Examples
		>>> from fmat import *
		
		>>> randmat(3, type=int)
		Mat(
		┌  9  1  10 ┐
		│ 10  1   9 │
		└  8  2   0 ┘)
		>>> _.rowinsert(randmat(1, 3, type=int))
		Mat(
		┌  9  1  10 ┐
		│ 10  1   9 │
		│  8  2   0 │
		└  2  6   7 ┘)
		>>> _
		Mat(
		┌  9  1  10 ┐
		│ 10  1   9 │
		│  8  2   0 │
		└  2  6   7 ┘)
		'''
		self.__data.extend(other.__data)
		self.__rows += other.__rows
		self.__shape = (self.__rows, self.__cols)
		return self
	
	def colinsert(self, other) :
		'''To insert data in place at the end of each row.
		Notice : Returns the inserted inplacely matrix.
		
		Examples
		>>> from fmat import *
		
		>>> randmat(3)
		Mat(
		┌                  5    6.44329976935601  0.6235375848663727 ┐
		│ 2.0644988208476134  1.6292535038425615                   5 │
		└ 0.9769444170797714                   5   8.283193736588784 ┘)
		>>> _.colinsert(randmat(3, 1))
		Mat(
		┌                  5    6.44329976935601  0.6235375848663727   8.558565813026624 ┐
		│ 2.0644988208476134  1.6292535038425615                   5  7.6257670984834345 │
		└ 0.9769444170797714                   5   8.283193736588784   5.511601362818372 ┘)
		>>> _
		Mat(
		┌                  5    6.44329976935601  0.6235375848663727   8.558565813026624 ┐
		│ 2.0644988208476134  1.6292535038425615                   5  7.6257670984834345 │
		└ 0.9769444170797714                   5   8.283193736588784   5.511601362818372 ┘)
		'''
		for i, j in zip(self.__data, other.__data) :
			i.extend(j)
		self.__cols += other.__cols
		self.__shape = (self.__rows, self.__cols)
		return self
	
	rows, cols, shape  = (
		property(lambda self :  self.__rows),
		property(lambda self :  self.__cols),
		property(lambda self : self.__shape))
		#property(lambda self : (self.__rows, self.__cols)))

	def copy(self) :
		return Mat([i[:] for i in self.__data])
		
	__deepcopy__ = __copy__ = copy
	