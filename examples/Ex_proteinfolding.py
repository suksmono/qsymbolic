"""
------------------------------------------------------------
 symbolic computation of protein folding problem
 Created by: Suksmono@{STEI-ITB, MDR Inc.}
------------------------------------------------------------
 """
from sympy import *
import numpy as np
from qsymbolic import *
#from qsymbolic import *
#from prb_simulation import *

# define used symbols
NQ=4
q=symbols('q0:%d'%NQ)

'''
------------------------------------------------------------
    1. USER WRITES THE PROBLEM
------------------------------------------------------------
'''
# Hamiltonian of the problem
Hq = -1 -4*q[2] + 9*q[0]*q[2] + 9*q[1]*q[2] \
     - 16*q[0]*q[1]*q[2]
 
'''
------------------------------------------------------------
    2. DO SYMBOLIC COMPUTATION AND TRANSFORM THE PROBLEM 
       FROM Q-DOMAIN TO S-DOMAIN
------------------------------------------------------------
'''

# define substition of k-body -> 2-body interaction
qSub=[[1,2,3]] # q1*q2->q3

H2b, H2s = qsym.q2s_symbolic(Hq, qSub, NQ)
print('H2b= \n',H2b);
print('H2s= \n',H2s);

'''
------------------------------------------------------------
    3. EXTRACT ISING PARAMETERS FROM SYMBOLIC SOLUTION
------------------------------------------------------------
'''
b, hi, Jij = qsym.isingCoeffs(H2s,NQ)

'''
------------------------------------------------------------
    4. OBTAIN A SOLUTION BY SIMULATION (S.A.)
------------------------------------------------------------
'''

SUBITER=100
MAXITER=1000
# calculate the solution using simulated annealing prb_sa
threshold, vMax=qsym.minMaxHq(Hq, NQ)
vSpin=qsimu.prb_sa(b, hi, Jij, NQ, SUBITER, MAXITER, threshold)

'''
 ------------------------------------------------------------
     5. DISPLAY THE RESULTS
 ------------------------------------------------------------
 '''
# qubits
qSpin=qsym.vs2q(vSpin)
print('\n******************************')
print('            RESULTS')
print('******************************\n')
print('q0=', int(qSpin[0]),',q1=', int(qSpin[1]), \
      ',q2=', int(qSpin[2]), ',q3=', int(qSpin[3]), )

print('\n******************************')
print('  Parameter for D-Wave input')
print('     (to be normalized)')
print('******************************\n')
print('Bias b=',b)
print('hi=',hi)
print('Jij=\n',Jij)
