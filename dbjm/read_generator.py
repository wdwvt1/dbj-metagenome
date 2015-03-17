#!/usr/bin/env python
from numpy import array, arange
from numpy.random import rand, randint


def generate_read_length(lb, ub):
    return randint(lb, ub)

def base_freq_fxn(a,c,t,g):
    '''return function.'''
    def _x(x):
        if x < a:
            return 'A'
        elif a <= x < a + c:
            return 'C'
        elif a + c <= x < a + c + t:
            return 'T'
        elif a + c + t <= x:
            return 'G'
    return _x

def generate_read(arr, base_freq_fxn):
    '''make a read from an array of floats in [0, 1].'''
    # apply_along_axis(f, 0, q)
    # truth value error
    return ''.join([base_freq_fxn(i) for i in arr])

def _make_equiprobable_reads(rl_lb, rl_ub, num_reads):
    '''convenience function for generating reads with equal nt frequencies.'''
    _f = base_freq_fxn(.25, .25, .25, .25)
    reads = [generate_read(rand(generate_read_length(rl_lb, rl_ub)), _f) for _
             in range(num_reads)]
    return reads