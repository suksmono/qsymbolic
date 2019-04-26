"""
prb_simulation 
A SA (simulated annealing) module
"""
import numpy as np
from random import randint
from sympy import *
'''
---------------------------------------------
generate M-length binary vector {-1,+1}
---------------------------------------------
'''
def vRandSpin(M):
    return(np.sign(np.random.randn(M)))
    
'''
---------------------------------------------
generate annealing schedule for N-iteration
---------------------------------------------
'''
def pMCSchedule(N,pStart):
    N025=int(N/4);
    tSched=np.ones(N)
    tSched[3*N025:N-1]=1.0
    tSched[2*N025: 3*N025-1]=0.9999
    tSched[N025:2*N025+1]=0.999  
    dGap=1.0-pStart
    xt1=np.arange(0,N025-1,1)
    st1=np.exp(-xt1/N025)
    st1=pStart+dGap*(1-st1)/max(1-st1)
    tSched[0:N025-1]=st1
    return(tSched)

'''
---------------------------------------------
Energy function 
---------------------------------------------
'''
def vEnergyBx(b, hi, Jij,vSpin):
    E=b+np.dot(hi,vSpin)+np.dot(vSpin,np.dot(Jij,vSpin))
    return(E)

'''
------------------------------------------------------------
SIMULATED ANNEALING
------------------------------------------------------------
'''
def prb_sa(b, hi, Jij, NQ, SUBITER, MAXITER, Etresh):   
    ITR1=SUBITER;
    vStart=0;
    # random vector
    vSpin0=vRandSpin(NQ);
    
    maxIter=MAXITER #1*10*100;
    Ptresh=pMCSchedule(maxIter,0.5)
    #-- Ptresh is the schedule
    k=0;
    vSpin=vSpin0;
    # init energy setting
    Ecurr = vEnergyBx(b, hi, Jij,vSpin);
    #Etresh=-5; #0;
    while (k<maxIter) and (Ecurr>Etresh):
        k=k+1;
        Eprev=Ecurr;
        # try Monte Carlo
        # --flip a row randomly
        RD=randint(0,NQ-1);
        # copy current spin to a template
        tvSpin=vSpin
        tvSpin[RD] = -1*tvSpin[RD]; #flip spin
        Ecurr=vEnergyBx(b, hi, Jij,tvSpin);
        if (Ecurr<Eprev) or (np.random.rand()>Ptresh[k]): 
            vSpin=tvSpin;
        Esys=vEnergyBx(b, hi, Jij,vSpin);
        Ecurr=Esys;
        print('Iter:',k,'of',maxIter,'E=', Esys)
    return(vSpin)