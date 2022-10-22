#!/usr/bin/env python3
#
# Name: test
# Author: BiaoLin Liao
#

__all__ = [
	'gt',
	'gprint',
	'main',
]

from .ctx import mst

@mst(
	#pretty_print = False,
	accuracy_protect = mst.accuracy_protect,
	eps = 1e-14,
)
def main() :
	from timeit import timeit
	from fractions import Fraction, Decimal
	from .numbers import Rational
	from .core import (
		gcd, lcm, ratioSimp, Mat,
		zeros, ones, diag, eye, randmat, hilbert)
	from pprint import pprint
	
	a = Mat(
		[[1, 2, 2,  2],
		 [2, 4, 6,  8],
		 [3, 6, 8, 10]])
	#a = randmat(20, min=-50, max=200, type=int)
	
	print('a.__sizeof__()', a.__sizeof__(), '', sep='\n')
	print('a:', a, '', sep='\n')
	print('rows: {}, cols: {}'.format(a.rows, a.cols), '', sep='\n')
	print('shape: {}'.format(a.shape), '', sep='\n')
	print('rank: {}'.format(a.rank()), '', sep='\n')
	print('plus one:', a+1, '', sep='\n')
	print('minus one:', a-1, '', sep='\n')
	print('2 mul and mul 2:', 2*a*2, '', sep='\n')
	print('get a.T:', a.T, '', sep='\n')
	print('get a.H:', a.H, '', sep='\n')
	print('mul a.T:', a*a.T, '', sep='\n')
	print('get a.kron(a):', a.kron(a), '', sep='\n')
	print('get a.rref():', a.rref(), '', sep='\n')
	print('get a.columnspace():', a.columnspace(), '', sep='\n')
	print('get a.nullspace():', a.nullspace(), '', sep='\n')
	print('get a.rref() or a.nullspace() don\'t change a:', a, '', sep='\n')
	print('get a⁻¹ (a**-1):', a**-1, '', sep='\n')
	print('try a * a⁻¹ (a**-1):', a * a**-1, '', sep='\n')
	print('get a.det():', a.det(), '', sep='\n')
	print('get randmat(3):', randmat(3), '', sep='\n')
	
	def getOptimizationProblem(C, A, b, maximize=False, ILP=False) :
		g = 'max z = ' if maximize else 'min z = '
		SubscriptDict = {i: j for i, j in zip('0123456789.', '₀₁₂₃₄₅₆₇₈₉⡀')}
		toSubscript, checkCoe = (
			lambda n : ''.join(SubscriptDict[i] for i in str(n)),
			lambda n : '' if n == 1 else '-' if n == -1 else n)
		
		f = '+ '.join(f'{checkCoe(j)}x{toSubscript(i)}' for i,j in enumerate(C, 1))
		c = 's.t.\t' + \
			'\n\t\t'.join(
			'+ '.join(f'{checkCoe(u)}x{toSubscript(k)}' for k,u in enumerate(i, 1) if u)\
			+ f' ≤ {j}' for i,j in zip(A.tolist(), b))
		r = 'xₙ ∈ Integers ∧ '*int(ILP) + 'xₙ ≥ 0'
		
		return g + f'{f}\n{c}\n\t\t{r}\n'.replace('+ -', '- ')
		
	C = Mat([[2, 1]])
	A =Mat([[0, 5],
			[6, 2],
			[1, 1]])
	b = Mat([15, 24, 5])
	s = C.simplex(A, b, maximize=True)
	print(getOptimizationProblem(C, A, b, maximize=True, ILP=False),
			'solve the LP problem:', s, '', sep='\n')
			
	C = Mat([[3, 1, 3]])
	A =Mat([[  2, 2, 1],
			[1.5, 3, 3],
			[  2, 1, 1]])
	b = Mat([30, 25, 20])
	s = C.relax_simplex(A, b, maximize=True)
	print(getOptimizationProblem(C, A, b, maximize=True, ILP=False),
			'solve the LP problem:', s, '', sep='\n')
	
	C = Mat([[5, 8]])
	A =Mat([[1, 1],
			[5, 9]])
	b = Mat([6, 45])
	s = gt(mst(accuracy_protect=None)(C.branch_bound), 1)(A, b, True)
	print(getOptimizationProblem(C, A, b, True, True),
			'solve the ILP problem:', s, '', sep='\n')

	C = Mat([[1, 2, 3]])
	A =Mat([[-1,  1, -1],
			[ 1,  1,  2],
			[ 0, -1,  1]])
	b = Mat([-4, 8, -2])
	s = C.branch_cut(A, b, False)
	print(getOptimizationProblem(C, A, b, False, True),
			'solve the ILP problem:', s, '', sep='\n')
			
	C = Mat([[-5, -4, -6]])
	A =Mat([[1, -1, 1],
			[3,  2, 4],
			[3,  2, 0]])
	b = Mat([20, 42, 30])
	s = C.simplex(A, b, maximize=False)
	print(getOptimizationProblem(C, A, b, False, False),
			'solve the LP problem:', s, '', sep='\n')
	
	C = Mat([[-1, 1]])
	A=-Mat([[-2, 7],
			[-2, 5],
			[ 0, 1],
			[ 1, 0],
			[ 0, 1]])
	b = -Mat([1]*5)
	s = C.simplex(A, b, maximize=False)
	print(getOptimizationProblem(C, A, b, False, False),
			'solve the LP problem:', s, '', sep='\n')
			
	return
	A = Mat(
		[[ -386043000,   659181000],
		 [ -108067536,   203776560],
		 [ -117749241,           0],
		 [ -617668800,  1054689600],
		 [   29253552,  -112360080],
		 [ -377105244,   513014610],
		 [  -27016884,    50944140],
		 [ -188552622,   256507305],
		 [-2099840160,  3711806196],
		 [    9751184,   -37453360],
		 [  -54033768,   101888280],
		 [ -386043000,   659181000],
		 [-1235337600,  2109379200],
		 [  -78499494,           0],
		 [-2099840160,  3711806196],
		 [ -193021500,   329590500],
		 [ -486800223,   831227241],
		 [  -64340500,   109863500],
		 [ -425292747,   659181000],
		 [ -814439310,  1125299550],
		 [ -864540288,  1630212480],
		 [  -78499494,           0],
		 [-1049920080,  1855903098],
		 [    2437796,    -9363340],
		 [-2099840160,  3711806196],
		 [  -78499494,           0],
		 [ -108067536,   203776560],
		 [  -12390108,    -5235900],
		 [-1508420976,  2052058440],
		 [ -188552622,   256507305],
		 [ -128681000,   219727000],
		 [ -386043000,   659181000],
		 [-1235337600,  2109379200],
		 [  -39249747,           0],
		 [-1003711800,  1713870600],
		 [  -54033768,   101888280],
		 [-2438255160,  5108211180],
		 [ -386043000,   659181000],
		 [          0, -1569989880],
		 [ -117749241,           0]])[22:, :]
	C = ones(1, A.cols)
	b =-ones(1, A.rows)
	
	print(getOptimizationProblem(C, A, b, False, True))

	print(gt(C.branch_cut)(A, b, maximize=False), '\n')
	
	print(gt(C.branch_bound)(A, b, maximize=False))

from os import write as osWt
from time import perf_counter
from sys import stdout

class gt(object) :
	__slots__ = (
		'__t_s',
		'__tag',
		'__lp',
	)
	
	def __init__(self, tag=None, loop=1) :
		self.__tag, self.__lp = \
			('GET TIME' if tag is None else tag), loop
		
	def __call__(self, *args, **kwds) :
		sloop, func = self.__lp>1, self.__tag
		try :
			self.__t_s = perf_counter()
			res = func(*args, **kwds)
		except :
			self.__t_s = perf_counter() - self.__t_s
			try : n = self.__tag.__name__
			except : n = self.__tag
			gt.putout('%s: %.9fs' % (n, self.__t_s))
			raise
		if sloop :
			for _ in range(self.__lp-1) :
				func(*args, **kwds)
		self.__t_s = perf_counter() - self.__t_s
		try : n = self.__tag.__name__
		except : n = self.__tag
		gt.putout('%s: %.9fs' % (n, self.__t_s))
		return res

	def __enter__(self) :
		gt.putout('↓%s↓' % self.__tag)
		self.__t_s = perf_counter()
		return self.__tag
	
	def __exit__(self, exc_type, exc_val, exc_tb) :
		self.__t_s = perf_counter() - self.__t_s
		if self.__lp != 1 :
			raise NotImplemented
		if exc_type is None :
			gt.putout('↑%s→ %.9fs' % (self.__tag, self.__t_s))
		else :
			gt.putout('↑%s→ ERROR' % self.__tag)
		return False
	
	@classmethod
	def putout(cls, *objs, sep=' ', end='\n', flush=False, func=str) :
		osWt(1, sep.join(map(func, objs)).encode()+end.encode())
		if flush : stdout.flush()
	
	@classmethod
	def seeup(cls, obj, end='\n', flush=False, func=str) :
		osWt(1, func(obj).encode()+end.encode())
		if flush : stdout.flush()
		return obj
	
	gt = lambda s : perf_counter() - self.__t_s
	
	def reset(self) :
		self.__t_s = perf_counter()
		
	__name__ = \
		property(lambda s :
			s.__class__.__name__
			if type(s.__tag) in (None, str)
			else s.__tag.__name__)

gprint = gt.putout
