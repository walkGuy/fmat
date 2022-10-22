#!/usr/bin/env python3

__version__ = '1.0.0'

__title__ = 'fmat'
__description__ = 'A Simple Matrix Realization'

__author__ = 'Baolin Liao'
__email__ = 'Not yet'

__license__ = 'LGPL v3'
__copyright__ = 'Copyright (c) 2022 Baolin Liao'

from .core import (
	mst, gcd, lcm, ratioSimp, Mat,
	zeros, ones, diag, eye, randmat, hilbert)

__all__ = [
	'mst',
	'gcd',
	'lcm',
	'ratioSimp',
	
	'Mat',
	'zeros',
	'ones',
	'diag',
	'eye',
	'randmat',
	'hilbert',
	
	'__version__',
	'__title__',
	'__description__',
	'__author__',
	'__email__',
	'__license__',
	'__copyright__'
]
