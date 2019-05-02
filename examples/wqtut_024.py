# -*- coding: utf-8 -*-
"""
Created on Thu May  2 12:09:13 2019
@author: Suksmono
Protein folding
"""
from qsymbolic import *

# create an object qsymaqc 
mypb=qsymb.qsymaqc()

# fill EQ with expression in STRING format
mypb.EQ='-1-4*q3+9*q1*q3+9*q2*q3-16*q1*q2*q3'

mypb.delta=10 # qubo-type problem
mypb.var_type='q' # qubo-type problem

# extract qubo/ising coefficients
hi, Jij = mypb.extract_qs_coeffs_from_string() 
## embed h into qubo matrix
J=Jij.tolist()
for m in range(0,len(hi)):
    J[m][m]=hi[m]

## to be processed by neal, wildqat, qbsolve, dwave, or other quantum machines
#print('\nhi=',mypb.ising_hi)
#print('\nJij=',mypb.ising_Jij)

###
'''
Solve the problem using wildqat
'''
### ---    
print('\nSolve using wildqat ... !')
import wildqat as wq
b = wq.opt()
# fill in qubo coeffs
b.qubo =J 
print(b.sa())