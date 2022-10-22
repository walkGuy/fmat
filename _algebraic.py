#!/usr/bin/env python3

__all__ = [
	'Rational',
	'Radical',
	'RealNumber',
]

import re
from math import gcd, floor

def to_ratio(integer, start=None, loop=None) :
	'''
	Examples
		Turns 6.97 to a integer ratio
	
	>>> to_ratio(6, 97)
	(697, 100)
				 .  .
		Turns 23.6877 to a integer ratio
		
	>>> to_ratio(23, None, 6877) # way 1
	(236854, 9999)
	
	>>> to_ratio(23, loop=6877) # way 2
	(236854, 9999)
	
					.  .
		Turns 15.28615886 to a integer ratio
		
	>>> to_ratio(15, 2861, 5886)
	(6793169, 444400)
	'''
	s = str(start)
	if start is None :
		s, start = '', 0
	d = 10**len(s)
	if loop :
		loop = str(loop)
		d *= 10**len(loop)-1
		n = int(s+loop)-start + integer*d
	else :
		n = start + integer*d
	try :
		s = gcd(n, d)
	except :
		s, b = n, d
		while b : s, b = b, s%b
	return (n//s, d//s)
	
_RATIONAL_FORMAT = re.compile(r'''
	\A\s*					  # optional whitespace at the start, then
	(?P<sign>[-+]?)			# an optional sign, then
	(?=\d|\.\d)				# lookahead for digit or .digit
	(?P<num>\d*)			   # numerator (possibly empty)
	(?:						# followed by
	   (?:/(?P<denom>\d+))?	# an optional denominator
	|						  # or
	   (?:\.(?P<decimal>\d*))? # an optional fractional part
	   (?:E(?P<exp>[-+]?\d+))? # and optional exponent
	)
	\s*\Z					  # and optional whitespace to finish
''', re.VERBOSE | re.IGNORECASE)


class Rational(object) :
	__slots__ = (
		'main', # A tuple[int] or a integer number.
		'_ratio' # True or False.
	)
	
	@classmethod
	def _new_f(cls, n: int) :
		self = super().__new__(cls)
		self.main, self._ratio = n, False
		return self
		
	@classmethod
	def _new_t(cls, n, d, simplify=False) :
		self = super().__new__(cls)
		
		if simplify :
			if type(n) is int is type(d) :
				c = gcd(n, d)
				n, d = n//c, d//c
			else :
				a, b = n, d
				while b : a, b = b, a%b
				n, d = int(n/a), int(d/a)
				
		self.main, self._ratio = \
			(n, d), True
		return self
		
	@classmethod
	def _new_n(cls, n) : # n could not be str
		self = super().__new__(cls)
		if type(n) is int :
			self.main, self._ratio = n, False
			return self
		
		try :
			n, d = n.as_integer_ratio()
		except :
			a, d = n, 1
			while d : a, d = d, a%d
			n, d = int(n/a), int(1/a)
		
		self.main, self._ratio = \
			(n*d, False) if abs(d) == 1 else ((n, d), True)
		return self
		
	def __new__(cls, n, d=None, *, simplify=True) :
		self = super().__new__(cls)
			
		if d is None :
			# the only n entered status
			if type(n) is int :
				# if is int that just set n to main
				self.main, self._ratio = n, False
				return self
			elif type(n) is str :
				# if n is str
				m = _RATIONAL_FORMAT.match(n)
				if m is None :
					raise ValueError('Invalid literal for Rational: %r' % n)
				n, denom = int(m.group('num') or '0'), m.group('denom')
				if denom :
					d = int(denom)
				else :
					decimal, exp = \
						m.group('decimal'), m.group('exp')
					
					if decimal :
						scale = 10**len(decimal)
						n, d = n*scale + int(decimal), scale
					else :
						d = 1
						
					if exp :
						exp = int(exp)
						if exp >= 0 :
							n *= 10**exp
						else:
							d *= 10**-exp
				if m.group('sign') == '-' :
					n = -n
				
				if simplify :
					c = gcd(n, d)
					n, d = n//c, d//c
			else :
				# if n is not integer or str
				try :
					# if n is float that can just get ratio
					n, d = n.as_integer_ratio()
				except :
					# it could not work that gcd by Euclid method
					a, d = n, 1
					while d : a, d = d, a%d
					n, d = int(n/a), int(1/a)
		else :
			# both n and d are entered
			if type(n) is int is type(d) and simplify :
				# both them are integer
				c = gcd(n, d)
				n, d = n//c, d//c
			else :
				# gcd by Euclid method if is not integer
				a, b = n, d
				while b : a, b = b, a%b
				n, d = int(n/a), int(d/a)
				
		if d < 0 :
			if d == -1 :
				self.main, self._ratio = -n, False
			else :
				self.main, self._ratio = (-n, -d), True
		elif d == 1 :
			self.main, self._ratio = n, False
		elif not d :
			raise ZeroDivisionError('Rational(%s, 0)' % numerator)
		else :
			self.main, self._ratio = (n, d), True
			
		return self
	
	def __add__(self, other) :
		if self._ratio :
			ds = self.main[1]
			if isinstance(other, Rational) :
				if other._ratio :
					do = other.main[1]
					return Rational(
						do * self.main[0]  + ds * other.main[0],
						ds * do)
				return Rational(
					self.main[0]  + ds * other.main,
					ds)
			return Rational(
				self.main[0]  + ds * other,
				ds)
		if isinstance(other, Rational) :
			if other._ratio :
				do = other.main[1]
				return Rational(
					do * self.main  + other.main[0],
					do)
			return Rational._new_n(self.main + other.main)
		return Rational._new_n(self.main + other)
	__radd__ = __add__
		
	def __sub__(self, other) :
		if self._ratio :
			ds = self.main[1]
			if isinstance(other, Rational) :
				if other._ratio :
					do = other.main[1]
					return Rational(
						do * self.main[0]  - ds * other.main[0],
						ds * do)
				return Rational(
					self.main[0]  - ds * other.main,
					ds)
			return Rational(
				self.main[0]  - ds * other,
				ds)
		if isinstance(other, Rational) :
			if other._ratio :
				do = other.main[1]
				return Rational(
					do * self.main  - other.main[0],
					do)
			return Rational._new_n(self.main - other.main)
		return Rational._new_n(self.main - other)
	__rsub__ = lambda self, other : -self - other
	
	def __mul__(self, other) :
		if self._ratio :
			if isinstance(other, Rational) :
				if other._ratio :
					return Rational(
						self.main[0]*other.main[0],
						self.main[1]*other.main[1])
				return Rational(
					self.main[0]*other.main,
					self.main[1])
			return Rational(
				self.main[0]*other,
				self.main[1])
		if isinstance(other, Rational) :
			if other._ratio :
				return Rational(
					self.main*other.main[0],
					other.main[1])
			return Rational._new_n(other.main*self.main)
		return Rational._new_n(self.main*other)
	__rmul__ = __mul__
		
	def __truediv__(self, other) :
		if self._ratio :
			if isinstance(other, Rational) :
				if other._ratio :
					return Rational(
						self.main[0]*other.main[1],
						self.main[1]*other.main[0])
				return Rational(
					self.main[0],
					self.main[1]*other)
			return Rational(
				self.main[0],
				self.main[1]*other)
		if isinstance(other, Rational) :
			if other._ratio :
				return Rational(
					self.main*other.main[1],
					other.main[0])
			return Rational(self.main, other.main)
		return Rational(self.main, other)
		
	__rtruediv__ = lambda self, other : \
		Rational._new_n(other) / self
	
	__floordiv__ = lambda self, other : \
		floor(self / other)
		
	__rfloordiv__ = lambda self, other : \
		floor(other / self)
	
	__mod__ = lambda self, other : \
		self - other * floor(self / other)
	
	__rmod__ = lambda self, other : \
		Rational._new_n(other) % self
	
	def __pow__(self, other) :
		if isinstance(other, Rational) :
			if other._ratio :
				return self ** (other.main[0] / other.main[1])
			return self ** other.main
			
		if self._ratio :
			if other < 0 :
				return Rational(
					self.main[1]**-other,
					self.main[0]**-other)
			return Rational(
				self.main[0]**other,
				self.main[1]**other)
				
		if other < 0 :
			return Rational(1, self.main**-other)
		return Rational._new_n(self.main**other)

	__rpow__ = lambda self, other : \
		Rational._new_n(other)**self
	
	def __pos__(self) :
		if self._ratio :
			return Rational._new_t(*self.main)
		return Rational._new_f(self.main)
		
	def __neg__(self) :
		if self._ratio :
			n, d = self.main
			return Rational._new_t(-n, d)
		return Rational._new_f(-self.main)
	
	def __trunc__(self) :
		if self._ratio :
			if self.main[0] < 0 :
				return self.__ceil__()
			return self.__floor__()
		return self.main
	
	def __floor__(self) :
		'''Will be math.floor(a) in Python 3.'''
		if self._ratio :
			return (lambda a,b : a // b)(*self.main)
		return self.main

	def __ceil__(self) :
		'''Will be math.ceil(a) in Python 3.'''
		# The negations cleverly convince floordiv to return the ceiling.
		if self._ratio :
			return -(lambda a,b : -a // b)(*self.main)
		return self.main
		
	def __round__(self, ndigits=None):
		'''Will be round(self, ndigits) in Python 3.

		Rounds half toward even.
		'''
		if ndigits is None:
			floor, remainder = divmod(self.numerator, self.denominator)
			if (remainder << 1) < self.denominator \
					or not (floor & 1) :
				return floor
			return floor + 1
		shift = 10**abs(ndigits)
		if ndigits > 0 :
			return Rational(round(self * shift), shift)
		else:
			return Rational._new_n(round(self / shift) * shift)

	def __abs__(self) :
		if self._ratio :
			if self.main[0] < 0 :
				return Rational._new_t(-self.main[0], self.main[1])
			return self
		return abs(self.main)
		
	__hash__ = lambda self : hash(self.main)
	
	def _compare(self, other, op) :
		if isinstance(other, Rational) :
			if self._ratio :
				if other._ratio :
					return op(self.main[0]*other.main[1],
						other.main[0]*self.main[1])
				return op(
					self.main[0],
					other.main*self.main[1])
			if other._ratio :
				return op(
					self.main*other.main[1],
					other.main[0])
			return op(self.main, other.main)
			
		if self._ratio :
			return op(self, Rational._new_n(other))
		return op(self.main, other)
		
	__eq__ = lambda self, other : \
		self._compare(other, lambda a,b : a == b)
	
	__lt__ = lambda self, other : \
		self._compare(other, lambda a,b : a < b)

	__gt__ = lambda self, other : \
		self._compare(other, lambda a,b : a > b)
		
	__le__ = lambda self, other : \
		self._compare(other, lambda a,b : a <= b)

	__ge__ = lambda self, other : \
		self._compare(other, lambda a,b : a >= b)

	__bool__ = lambda self : bool(self.numerator)
	
	__float__ = lambda self : \
		self.main[0]/self.main[1] \
		if self._ratio else \
		float(self.main)
		
	def as_integer_ratio(self) :
		if self._ratio :
			return self.main
		return (self.main, 1)
	
	def __str__(self) :
		if self._ratio : return '{}/{}'.format(*self.main)
		else : return str(self.main)
	
	def __repr__(self) :
		return self.__class__.__name__ + '(' + str(self) + ')'
	
	def __n_get(self) :
		if self._ratio :
			return self.main[0]
		return self.main
		
	def __n_set(self, value) :
		if self._ratio : self.main[0] = value
		else : self.main = value
		
	def __d_get(self) :
		if self._ratio :
			return self.main[1]
		return 1
		
	def __d_set(self, value) :
		if abs(value) == 1 :
			self.main, self_ratio = \
				value*self.numerator, False
		elif self._ratio :
			self.main[1] = value
		else :
			self.main, self._ratio = \
				(self.main, value), True
		
	numerator, denominator = \
		property(__n_get, __n_set), \
		property(__d_get, __d_set)
		
	# support for pickling, copy, and deepcopy

	def __reduce__(self) :
		return (self.__class__, (str(self),))

	def __copy__(self) :
		if type(self) is Rational:
			return self
		return self.__class__(self.numerator, self.denominator)

	def __deepcopy__(self, memo) :
		if type(self) is Rational:
			return self
		return self.__class__(self.numerator, self.denominator)


class Radical(Rational) :
	__slots__ = (
		'main',
		'_status'
	)
	...

class RealNumber(object) :
	__slots__ = (
		'main',
		'_status'
	)
	
	