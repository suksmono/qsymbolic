# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 09:05:27 2019

@author: User
"""

# test using wildcat
from qsymbolic import *
import itertools as itr
from sympy import *
import neal
import numpy as np

# Hamiltonian in text/string form
#E0 = '2 + 3*x0 + 4*x0*x1 + 5*x0*x1*x2 + 7*x0*x1*x2*x3 + 5*x4- x1**2 + 2*x3**2'
#E0='2 + x0*x1*x2*x3'

mypb=qsymb.qsymaqc()
mypb.EQ='2 + s0*s1*s2*s3'
#mypb.var_type='q'
mypb.solve()
#
if(mypb.err_code<1):
    # minimum energies and configurations
    print('Minimum Energy: {:06.2f}'.format(mypb.min_energy),\
          '->should has been: {:06.2f}'.format(-mypb.ising_b))
    print('standardized variables', mypb.var_standard)
    print('original variables', mypb.var_original)
    print('*Minimum configurations:',mypb.min_configuration.tolist()[0][0:4])
    
'''
#### TRY WILDQAT
'''    
from wildqat import *
#
J=mypb.ising_Jij
for i in range(0,len(mypb.ising_hi)-1):
        #print(i,'-',j)
        J[i,i]=mypb.ising_hi[i]
####
a = opt()
a.qubo = J
result=a.run() #=> 
print('*Wildqat solution: ', result[0:4])