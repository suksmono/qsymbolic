"""
Created on Sat Apr 27 10:07:55 2019
@author: suksmono
"""

from qsymbolic import *
#import qsymbolic as qs
import itertools as itr
from sympy import *
import neal
import numpy as np

# Hamiltonian in text/string form
#E0 = '2 + 3*x0 + 4*x0*x1 + 5*x0*x1*x2 + 7*x0*x1*x2*x3 + 5*x4- x1**2 + 2*x3**2'
E0='2 + x0*x1*x2*x3'

# create an object qsymaqc 
mypb=qsymb.qsymaqc()

# set initial values
mypb.EQ=E0
#mypb.var_type='q'
mypb.verbose=True
mypb.N_SWEEPS=1000*1*10
# solve it
mypb.solve()
#
if(mypb.err_code<1):
    # minimum energies and configurations
    print('Minimum Energy: {:06.2f}'.format(mypb.min_energy),\
          '->should has been: {:06.2f}'.format(-mypb.ising_b))
    print('Minimum Configurations:',mypb.min_configuration)
    print('standardized variables', mypb.var_standard)
    print('original variables', mypb.var_original)



