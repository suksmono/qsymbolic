"""
------------------------------------------------------------
 symbolic compution of factorization problem
 Created by: Suksmono@{STEI-ITB, MDR Inc.}
 problem adopted from:
 https://qiita.com/YuichiroMinato/items/9c9fba2b5bc978ec9e73
 the idea: 
 (1) a user enter his/her problem
 (2) server do symbolic calculation and simulation: 
     (a)Pq inc k-body -> Ps with only 2-bdy in s-domain
     (b) extract coefficients {b, hi, Jij}
           -> fed to simulator
           -> deliver to user for D-Wave input
 (3) server send back simulation results and parameter to user
 (4) user run his problem (b, hi, Jij) in D-wave, compare the
     result from simulator
------------------------------------------------------------
 """
from sympy import *
import numpy as np
from qsymbolic import *
'''
we will calculate the Hamiltonian
H=(Nc- p*q)^2,where p=(1+2*q0+4*q1),q=(1+2*q2);Nc=15, 21 ...
------------------------------------------------------------
    1. USER WRITES THE PROBLEM
------------------------------------------------------------
'''
Nc=15 # a composite number to be factored


E0='(-(2*x2 + 1)*(2*x0 + 4*x1 + 1) + 15)**2'
mypb=qsymb.qsymaqc()
mypb.EQ=E0
#mypb.N_SWEEPS=1000000
mypb.var_type='q'
mypb.solve()
if(mypb.err_code<1):
    # minimum energies and configurations
    print('Minimum Energy: {:06.2f}'.format(mypb.min_energy),\
          '->should has been: {:06.2f}'.format(-mypb.ising_b))
    print('Minimum Configurations:',mypb.min_configuration)
    #
    q0=int(mypb.min_configuration.tolist()[0][0])
    q1=int(mypb.min_configuration.tolist()[0][1])
    q2=int(mypb.min_configuration.tolist()[0][2])
    #
    P=2*q2 + 1
    Q=2*q0 + 4*q1 + 1
    print(mypb.var_original[0],'=', q0,'; ',\
          mypb.var_original[1],'=', q1,'; ',\
          mypb.var_original[2],'=', q2 )  
    print('Prime factors of ', Nc, ' are', P, 'and', Q)

