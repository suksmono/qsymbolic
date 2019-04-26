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

# define used symbols
NQ=4
q=symbols('q0:%d'%NQ)

'''
we will calculate the Hamiltonian
H=(Nc- p*q)^2,where p=(1+2*q0+4*q1),q=(1+2*q2);Nc=15, 21 ...
------------------------------------------------------------
    1. USER WRITES THE PROBLEM
------------------------------------------------------------
'''
Nc=25  # a composite number to be factored
# it is known that the prime factors are odd numbers
P=1+2*q[0]+4*q[1]
Q=1+2*q[2]

# Hamiltonian of the problem
Hq=(Nc-P*Q)**2
 
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
vSpin=qsimu.prb_sa(b, hi, Jij, NQ, SUBITER, MAXITER, 0.)

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
#
vp=int(P.subs({q[0]:qSpin[0],q[1]:qSpin[1]}))
vq=int(Q.subs({q[2]:qSpin[2]}))
print('Prime factors of ', Nc, ' are', vp, 'and', vq)

"""
print('\n******************************')
print('  Parameter for D-Wave input')
print('     (to be normalized)')
print('******************************\n')
print('Bias b=',b)
print('hi=',hi)
print('Jij=\n',Jij)
"""