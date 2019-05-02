# -*- coding: utf-8 -*-
"""
Created on Thu May  2 12:09:13 2019
@author: Suksmono
Prime factorization
"""
from qsymbolic import *

# create an object qsymaqc 
mypb=qsymb.qsymaqc()

# fill EQ with expression in STRING format
mypb.EQ='(15-(1-2*x0+4*x1)*(1+2*x2))**2'

mypb.delta=10 # qubo-type problem
mypb.var_type='q' # qubo-type problem

# extract qubo/ising coefficients
hi, Jij = mypb.extract_qs_coeffs_from_string() 

## to be processed by neal, wildqat, qbsolve, dwave, or other quantum machines
#print('\nhi=',mypb.ising_hi)
#print('\nJij=',mypb.ising_Jij)

###
'''
Solve the problem using wildqat
'''
print('\nSolve using wildqat ... !')
import wildqat as wq
b = wq.opt()
# fill in qubo coeffs
b.qubo =Jij.tolist() 
print(b.sa())