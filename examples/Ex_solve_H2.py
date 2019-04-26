"""
Created on Wed Jul 11 14:41:15 2018
author: Suksmono@{STEI-ITB, MDR Inc.}
Solve the Hadamard problem
MODEL, index translation:
  s11 s12 s13 > 0 1 4  
   0   1   4   
  s21 s22 s23  > 2 3 5
   2   3   5   
Hamiltonian of H2 matrix :: calculated by other script
H2s: 
  4*s11 + 4*s12 - 8*s13 + 4*s21 + 4*s22 - 8*s23 
  + 2*s11*s12 - 4*s11*s13 - 4*s12*s13 + 2*s11*s21 
  + 2*s11*s22 + 2*s12*s21 - 4*s11*s23 + 2*s12*s22 
  - 4*s13*s21 - 4*s12*s23 - 4*s13*s22 + 8*s13*s23 
  + 2*s21*s22 - 4*s21*s23 - 4*s22*s23 + 16
>>>
H2s: 
  4*s0 + 4*s1 - 8*s4 + 4*s2 + 4*s3 - 8*s5 
  + 2*s0*s1 - 4*s0*s4 - 4*s1*s4 + 2*s0*s2 
  + 2*s0*s3 + 2*s1*s2 - 4*s0*s5 + 2*s1*s3 
  - 4*s4*s2 - 4*s1*s5 - 4*s4*s3 + 8*s4*s5 
  + 2*s2*s3 - 4*s2*s5 - 4*s3*s5 + 16 
"""
from qsymbolic import *
from sympy import *
import neal
import numpy as np

'''
***************************************************
1. USER DEFINE THE PROBLEM AS HAMILTONIAN (spin var)
***************************************************
'''
# number of qubits
NQ=6

s=symbols('s0:%d'%NQ)
# Hamiltonian
H2s=  4*s[0] + 4*s[1] - 8*s[4] + 4*s[2] + 4*s[3] - 8*s[5]  \
      + 2*s[0]*s[1] - 4*s[0]*s[4] - 4*s[1]*s[4] + 2*s[0]*s[2]  \
      + 2*s[0]*s[3] + 2*s[1]*s[2] - 4*s[0]*s[5] + 2*s[1]*s[3]  \
      - 4*s[4]*s[2] - 4*s[1]*s[5] - 4*s[4]*s[3] + 8*s[4]*s[5]  \
      + 2*s[2]*s[3] - 4*s[2]*s[5] - 4*s[3]*s[5] + 16 
'''
***************************************************
2. QSYMBOLIC SOLVE THE PROBLEM 
***************************************************
'''

'''
------------------------------------------------------------
2.a. CONVERT THE HAMILTONIAN TO 2-BODY S-domain HAMILTONIAN
------------------------------------------------------------
'''

'''
------------------------------------------------------------
2.b. EXTRACT ISING PARAMETERS FROM SYMBOLIC SOLUTION
------------------------------------------------------------
'''

b, hi, Jij = qsym.isingCoeffs(H2s,NQ)

'''
-----------------------------------------------------------------------------
2.c convert the problem into Ising coefficients
-----------------------------------------------------------------------------
'''
#in dictionary format
h={0:0}
J={(0,1):1}

for m in range(0,len(hi)):
    h[m]=hi[m]
    for n in range (m+1,len(hi)):
        J[m,n]=Jij[m,n]
    
'''
------------------------------------------------------------
2.c. SOLVE THE PROBLEM
------------------------------------------------------------
'''

'''
-----------------------------------------------------------------------------
select a solver
> dimod: ExaxtSolver
> neal:  SimulatedAnnealingSampler
-----------------------------------------------------------------------------
'''
#
solver=neal.SimulatedAnnealingSampler()
response = solver.sample_ising(h, J)
#print('Configurations:\n', response.samples_matrix)
#print('Energy:\n',response.data_vectors['energy'])
vE=response.data_vectors['energy']
aSol=response.samples_matrix
# selecet a solution
idxMinE=np.argmin(vE)
tSol=aSol[idxMinE]
vSol=tSol[0]
print('*** THE SOLUTION ***')
print('Minimum Energy:',vE[idxMinE],'->supposed to be:', -b)
print('Minimum Configurations:',vSol)