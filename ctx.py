#/!usr/bin/env python3

__all__ = [
	'mst',
	'MatCtx'
]

from .numbers import Rational
from functools import wraps as funcWraps


class MatCtx(object) :
	'''A context class for Mat class that can set some settings about math or others.
	
	.pretty_print=True : bool
		- If True, __str__ and __repr__ method
			will use unicode characters to format matrix.
		- If False, __str__ and __repr__ method
			will return a beautified list.
			
	.alignment=str.rjust : a method of str object
		- It defined how matrix alignment.
		- Can be str.ljust, str.center and str.rjust.
	
	.accuracy_protect=fmat.numbers.Rational : class
		- It is a data object that used to 
			making sure the highly computing accuracy.
		- When the matrix class is doing some complex
			compute(such as divide), it will try to 
			transformation the computing data object
			to accuracy_protect object, then compute 
			under it to ensure the computing accuracy.
		- If set it to None, that will not use any to
			change the data object for saving accuracy.
		- Had better not change it, usually cause 
			unknown bugs.
			
	.eps=0 : number
		- If compute got any error, try to set it to
			a small number.
			
				Setting the contexts :
		Fmat uses a global working contexts.
		The working precision is controlled by a context
			object called mst. When the contexts has been
			set, all Mat operations are carried out at
			that contexts.
			
	Examples
	
	>>> Mat([[1, 2], [3, 4]])
	Mat(
	┌ 1  2 ┐
	└ 3  4 ┘)

	>>> mst.pretty_print = False
	
	>>> Mat([[1, 2], [3, 4]])
	Mat(
	[[1, 2],
	 [3, 4]])
	
	>>> mst
	<--
		pretty_print = False
		alignment = <method 'rjust' of 'str' objects>
		accuracy_protect = <class 'fmat.numbers.Rational'>
		eps = 0
	-->
	
	Temporarily changing the context settings :
		It is often useful to change the setting during
			only part of a calculation. A way to temporarily
			increase the context settings and then restore
			it is as follows :
				mst.eps +=3
				# do something here
				mst.eps -= 3
				
		Meanwhile, mst also provides two better ways to
			change local context settings as shown in
			the following examples.
			
	Examples
	
		1. By with ...
		
	>>> with mst :
	... 	mst.accuracy_protect = None
	... 	Mat([7]) / 2
	... 
	Mat(
	[3.5])
	>>> Mat([7]) / 2
	Mat(
	[7/2])
	
		2. By function decorators
		
	>>> @mst(accuracy_protect=None)
	... def f() :
	... 	return Mat([7]) / 2
	...
	>>> f()
	Mat(
	[3.5])
	
	>>> Mat([7]) / 2
	Mat(
	[7/2])
	'''
	__slots__ = ('pretty_print', 'alignment', 
		'accuracy_protect', 'eps', '_orig')
		
	def __init__(ctx) :
		ctx.pretty_print = True
		ctx.alignment = str.rjust
		ctx.accuracy_protect = Rational
		ctx.eps = 0
		
		ctx._orig = None
		
	def __call__(ctx, **kwargs) :
		'''It is often useful to change the precision during only part of a calculation.
		This define a way to temporarily increase the precision and then restore.
		It works as follows :
			>>> @mst(accuracy_protect=None)
			... def f() :
			... 	return Mat([7]) / 2
			...
			>>> f()
			Mat(
			[3.5])
			
			>>> Mat([7]) / 2
			Mat(
			[7/2])
		'''
		def w(f) :
			def g(*args, **kwds) :
				origin = {}
				for k,v in kwargs.items() :
					origin[k] = getattr(ctx, k)
					setattr(ctx, k, v)
				try :
					return f(*args, **kwds)
				finally :
					for k,v in origin.items() :
						setattr(ctx, k, v)
			return funcWraps(f)(g)
		return w
		
	def __enter__(ctx) :
		'''The one other way.
		Example
		
		>>> mst.accuracy_protect
		<class 'fmat.numbers.Rational'>
		>>> with mst :
		...  mst.accuracy_protect = float
		...  Mat([7])/2
		... 
		Mat(
		[3.5])
		>>> mst.accuracy_protect
		<class 'fmat.numbers.Rational'>
		'''
		#class __ctx(object) :
		#	__slots__ = ctx.__slots__
		#	def __init__(self) :
		#		for i in self.__slots__ :
		#			setattr(self, i, getattr(ctx, i))
		#ctx._orig = __ctx()
		
		ctx._orig = tuple(getattr(ctx, i) for i in ctx.__slots__)
		
	def __exit__(ctx, exc_type, exc_val, exc_tb) :
		#for i in ctx._orig.__slots__ :
		#	setattr(ctx, i, getattr(ctx._orig, i))
		
		for i,j in zip(ctx.__slots__, ctx._orig) :
			setattr(ctx, i, j)
		return False
		
	__repr__ = __str__ = lambda ctx : \
		'<--\n' + '\n'.join(
			'\t{} = {}'.format(i, getattr(ctx, i))
			for i in ctx.__slots__
		) + '\n-->'
	
mst = MatCtx()
