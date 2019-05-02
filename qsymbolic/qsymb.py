'''
qsymb: new version of qsym
'''
from sympy import *
import itertools as itr
import numpy as np
import neal
'''
dependencies: sympy, itrtools, numpy
'''
"""
####################################################
def is_ho_terms(tTerm):
    '''
    check, whether a term's order > 2
    input : symbolic term, eg. q0*q1, s0*s1*s2 ...
    output: True/False
    '''
    if len (tTerm.free_symbols):
        return True
    else:
        return False
"""

####################################################
def is_ho_equation(tEq):
    '''
    check, whether an Equation has a term of order > 2
    input : symbolic Equation, eg. 1+q0*q1, 2+s3**2+s0*s1*s2 ...
    output: True/False
    '''
    oTerms=tEq.as_ordered_terms()
    #print('term', oTerms,'->free symbols:', oTerms[0].free_symbols)
    if len( oTerms[0].free_symbols )>2:    #order: decreasing
        return True
    else:
        return False
####################################################
    
class qsymaqc:
    '''
    object of symbolic quantum computing to solve by 
    AQC(adiabatic quantum computing). Processing flow depend on variable type
        's': the values of variables are {-1,+1}
                    H_string -> Hks-> Hkq -> H2q -> H2s->s_sol
        'q': the values variables are {0,1}
                    H_string -> Hkq->H2q -> H2s-> s_sol -> q_sol
    '''
    def __init__(self):
        # a blank string of problem hamiltonian
        self.EQ=''
        self.Ekx=''
        self.Ekq=''
        self.Hks=''
        self.Hkq=''
        self.H2s=''
        self.H2q=''
        # number of sweeps in SA
        self.NSWEEPS=1*5*10*10*1000  
        # number of rads in SA
        self.NREADS=10
        # output all configurations
        self.configurations=[]
        # output all energies
        self.energies=[]
        # output minimum configurations
        self.min_configuration=[]
        # output minimum energy
        self. min_energy=0
        # b values in Ising equation
        #self.b=0
        # list of original variables
        self.var_original=[]
        # list of standardized variables
        self.var_standard=[]
        self.var_type='s'
        self.err_code=0
        # k-body --> 2-body substitution parameters
        self.delta=100
        self.ising_J={}
        self.ising_h={}
        self.ising_b= 0.0
        self.ising_Jij=[]
        self.ising_hi=[]
        self.ising_maxCoeffs=[]
        #
        self.verbose=True
   
    def extract_ising_coeffs(self):
        '''
        given a 2-body hamiltonian in s-domain, estract Ising parameters
        '''
        print('Obtaining Ising coefficients ...')
        b, hi, Jij = get_ising_coeffs(self.H2s)
        maxCoeff=np.max([np.max(abs(hi)), np.max(abs(Jij))])
        self.ising_maxCoeff=maxCoeff # for renormalization when required
        #
        self.ising_hi=hi /maxCoeff #hi
        self.ising_Jij=Jij /maxCoeff #Jij
        
        # supposed to be the geound state energy
        self.ising_b=b 
        # put hi in the diagonal forms of Jij
        # conform with neal data format ->J
        h, J=embed_ising_coeffs(hi, Jij)
        self.ising_h=h
        self.ising_J=J
        #return h,J

    def extract_qs_coeffs(self):
        '''
        given a 2-body hamiltonian in q/s-domain, estract qubo/ising coefficients
        '''
        print('Obtaining coefficients ...')
        
        if (self.var_type.upper()=='S'):
            b, hi, Jij = get_qs_coeffs(self.H2s,'s')
        else:
            b, hi, Jij = get_qs_coeffs(self.H2q,'q')
        maxCoeff=np.max([np.max(abs(hi)), np.max(abs(Jij))])
        self.ising_maxCoeff=maxCoeff # for renormalization when required
        #
        self.ising_hi=hi /maxCoeff #hi
        self.ising_Jij=Jij /maxCoeff #Jij
        # supposed to be the geound state energy
        self.ising_b=b 
        # put hi in the diagonal forms of Jij
        # conform with neal data format ->J
        h, J=embed_ising_coeffs(hi, Jij)
        self.ising_h=h
        self.ising_J=J       
     
    def solve_s(self):
        '''
        solve problem with variable 's'
        '''            
        print('Problem with a type-s (spin -1,+1) variable');
        varOrg,varStd,Eks = eq_standardized(self.Ekx,'s')
        #
        self.var_original=varOrg
        self.var_standard=varStd
        print('problem in symbolic standard format',Eks)
        # remove squared variables
        #Hks=simplify_squares(expand(Eks))
        Hks=expand(simplify_sq_squares(expand(Eks),'s'))
        self.Hks=Hks
        if self.verbose:            
            print('Hks->',Hks)
        # k-body to 2-body conversion, assume that k-max is 4
        Hkq=Hks2Hkq(Hks)
        self.Hkq=Hkq
        H2q =boole_reduce(Hkq,self.delta)
        self.H2q=H2q
        H2s=Hkq2Hks(H2q)
        self.H2s=H2s
        #
        self.extract_ising_coeffs()
        '''
        -----------------------------------------------------------------------------
        4. SOLVE THE PROBLEM
        -----------------------------------------------------------------------------
        select a solver
        > dimod: ExaxtSolver
        > neal:  SimulatedAnnealingSampler
        '''
        #
        print('Solving the problem using neal  ...')
        solver=neal.SimulatedAnnealingSampler()
        response = solver.sample_ising(self.ising_h, self.ising_J, sweeps=self.NSWEEPS, \
                                       num_reads=self.NREADS)
        #
        vE=response.data_vectors['energy']
        aSol=response.samples_matrix
        
        #
        idxMinE=np.argmin(vE)
        tSol=aSol[idxMinE]
        # copy solution to class global variables
        self.min_configuration=tSol[0]
        self.min_energy=vE[idxMinE]
        self.energies=vE
        self.configurations=aSol
    
    def solve_q(self):
        '''
        ######################################################
        Ekq -> E2q -> E2s ->SOL(si)-> SOL(qi)            
        ######################################################
        '''
        print('Problem with a type-q (boolean 0/1) variable')
        varOrg,varStd,Ekq = eq_standardized(self.Ekx,'q')
        self.Ekq=Ekq
        #
        self.var_original=varOrg
        self.var_standard=varStd
        ##
        if self.verbose:            
            print('Problem in symbolic standard format: ',Ekq)
        # simplifying qi**2 <- 1
        #Hkq=simplify_squares(expand(Ekq))
        Hkq=expand(simplify_sq_squares(expand(Ekq),'q'))
        self.Hkq=Hkq
        if self.verbose:            
            print('Simplified Hkq=', Hkq)
        # Ekq-> E2q
        H2q =boole_reduce(Hkq,self.delta)
        self.H2q=H2q
        if self.verbose:            
            print('H2q=', H2q)
        # ALT-1: H2q->H2s, then solve using ising SA
        # ALT-2: solve H2q using q-doman SA solver
        # ----------------------
        # choose ALT-1: H2s
        # ----------------------
        H2s=Hkq2Hks(H2q)
        self.H2s=H2s
        if self.verbose:            
            print('H2s=', H2s)
        #
        self.extract_ising_coeffs()
        
        '''
        -----------------------------------------------------------------------------
        4. SOLVE THE PROBLEM
        -----------------------------------------------------------------------------
        select a solver
        > dimod: ExaxtSolver
        > neal:  SimulatedAnnealingSampler
        '''
        #
        print('Solving the problem using neal  ...')
        solver=neal.SimulatedAnnealingSampler()
        response = solver.sample_ising(self.ising_h, self.ising_J, \
                                       sweeps=self.NSWEEPS, \
                                       num_reads=self.NREADS)
        vE=response.data_vectors['energy']
        aSol=response.samples_matrix
        idxMinE=np.argmin(vE)
        tSol=aSol[idxMinE]
        # copy solution to global class variables
        self.min_configuration=vs2q(tSol[0])
        self.min_energy=vE[idxMinE]
        self.energies=vE
        self.configurations=vs2q(aSol)
 #
    '''
    to be replaced by extract_qs_coefss_from_string
    '''
    def extract_ising_from_string(self):
        '''
        a problem is given in string format
        '''
        Ekx=eqString2Symbolic(self.EQ)
        self.Ekx=Ekx
        print('problem in symbolic format',Ekx)
        if (self.var_type.upper()=='S'):
            print('Problem with a type-s (spin -1,+1) variable');
            varOrg,varStd,Eks = eq_standardized(self.Ekx,'s')
            #
            self.var_original=varOrg
            self.var_standard=varStd
            print('problem in symbolic standard format',Eks)
            # remove squared variables
            #Hks=simplify_squares(expand(Eks))
            Hks=expand(simplify_sq_squares(expand(Eks),'s'))           
            self.Hks=Hks
            if self.verbose:            
                print('Hks->',Hks)
            # k-body to 2-body conversion, assume that k-max is 4
            Hkq=Hks2Hkq(Hks)
            self.Hkq=Hkq
            H2q =boole_reduce(Hkq,self.delta)
            self.H2q=H2q
            H2s=Hkq2Hks(H2q)
            self.H2s=H2s
            self.extract_ising_coeffs()
        else:
            print('\nIsing coeffs only for s-domain problem')

    def extract_qs_coeffs_from_string(self):
        '''
        a problem is given in string format
        '''
        Ekx=eqString2Symbolic(self.EQ)
        self.Ekx=Ekx
        print('Problem in symbolic format',Ekx)
        if (self.var_type.upper()=='S'):
            print('Problem with a type-s (spin -1,+1) variable');
            varOrg,varStd,Eks = eq_standardized(self.Ekx,'s')
            #
            self.var_original=varOrg
            self.var_standard=varStd
            print('problem in symbolic standard format',Eks)
            # remove squared variables
            # Hks=simplify_squares(expand(Eks))
            Hks=expand(simplify_sq_squares(expand(Eks),'s'))           
            self.Hks=Hks
            if self.verbose:            
                print('Hks->',Hks)
            # k-body to 2-body conversion, assume that k-max is 4
            Hkq=Hks2Hkq(Hks)
            self.Hkq=Hkq
            H2q =boole_reduce(Hkq,self.delta)
            self.H2q=H2q
            H2s=Hkq2Hks(H2q)
            self.H2s=H2s
            self.extract_ising_coeffs()
        else:
            print('Problem with a type-q/qubo (spin 0,1) variable');
            varOrg,varStd,Ekq = eq_standardized(self.Ekx,'q')
            #
            self.var_original=varOrg
            self.var_standard=varStd
            print('problem in symbolic standard format',Ekq)
            # remove squared variables
            #####################################################################
            #Hkq=expand(simplify_squares(expand(Ekq)))
            Hkq=expand(simplify_sq_squares(expand(Ekq),'q'))
            self.Hkq=Hkq
            if self.verbose:            
                print('Hkq->',Hkq)
            # k-body to 2-body conversion, assume that k-max is 4
            H2q =boole_reduce(Hkq,self.delta)
            self.H2q=H2q
            self.extract_qs_coeffs()
        #####
        return self.ising_hi, self.ising_Jij
    ########################################################################
    def solve(self):
        '''
        solve a problem given in string format
        '''
        Ekx=eqString2Symbolic(self.EQ)
        self.Ekx=Ekx
        print('problem in symbolic format',Ekx)
        if (self.var_type.upper()=='S'):
            # type s-variables
            self.solve_s()
            self.err_code=0
        elif (self.var_type.upper()=='Q'):
            # type q-variables
            self.solve_q()
            self.err_code=0
        else: 
            # other types are invalid
            print('Invalid problem type !!! \nShould have been either "s" or "q"');
            self.err_code=1                 
#
def embed_ising_coeffs(hi, Jij):
    '''
    -----------------------------------------------------------------------------
    convert the problem into Ising coefficients
    -----------------------------------------------------------------------------
    '''
    #in dictionary format
    h={0:0}
    J={(0,1):1}
    
    for m in range(0,len(hi)):
        h[m]=hi[m]
        for n in range (m+1,len(hi)):
            J[m,n]=Jij[m,n]
    #
    return h, J

def genSubtitutionPairs(hoT,idxStart):
    '''
    input: hoT-a set of high order variables 
           idxStart -start index of ancillary variables
    output: list of substitution pair, [q1,q2,q3] means:  (q1*q2)<--q3
    '''
    #
    lsubP=[]
    #cp=set(qs.genSubtitutionPairs(hoT,nvTO))
    cp=set(itr.combinations(hoT,2))
    idx=idxStart
    for tz in cp:
        ttz=list(tz)
        # add ancillary variable to var pair->triple
        ttz.append(symbols('q%d'%idx))
        idx=idx+1
        lsubP.append(ttz)
    return lsubP

def H2sub(x1,x2,y,d12):
    '''
    input:  x1, x2 of x1*x2 product
            y is var result x1*x2 <-- y
            d12 compensation factor
    output:two-body polynomial compensation term
    '''
    return(d12*(3*y+x1*x2-2*x1*y-2*x2*y))

'''
There are two binary representation
    s={-1,+1} and q={0,1}
The default in Ising simulation is s-domain
------------------------------------------------------------
symbolic: q-to-s and s-t-q transform
------------------------------------------------------------
'''
def q2s(x):
    return(1/2 -x/2)

def s2q(x):
    return(1 - 2*x)

# define function vq2s
def vq2s(x):
    return(1-2*x)
# define function vs2q
def vs2q(x):
    return(1/2-x/2)

def Hks2Hkq(Hks):
    '''
    ------------------------------------------------------------
     define Hks->Hkq transform as function
    ------------------------------------------------------------
    input: Hks
    output: Hkq
    '''
    # identify free symbols
    fsS=Hks.free_symbols
    Hkq=Hks
    for tx in fsS:
        #print(ts)
        #get symbol index 
        strTx=str(tx)
        sydex=int(strTx[1:])
        syq=symbols('q%d'%sydex)
        Hkq= Hkq.subs({tx:s2q(syq)}) 
    return Hkq.expand()
#
def Hkq2Hks(Hkq):
    '''
    ------------------------------------------------------------
     define Hkq->Hks transform as function
    ------------------------------------------------------------
    input: Hkq
    output: Hks
    '''
    # identify free symbols
    fsQ=Hkq.free_symbols
    Hks=Hkq
    for tx in fsQ:
        #print(ts)
        #get symbol index 
        strTx=str(tx)
        sydex=int(strTx[1:])
        sys=symbols('s%d'%sydex)
        Hks= Hks.subs({tx:q2s(sys)}) 
    return Hks.expand()
#
def boole_reduce(Hkq,delta):
    '''
    input:  k-body hamiltonian in q-domain Hkq
    output: 2-body hamiltonian in q-domain H2s
    stages: Hkq -> H2q
    '''
    # list all involved variables
    # define higher order set >2 
    hoT=set() 
    # collection all variables in Hk
    toT=set()
    # convert Hks->Hkq
    #Hkq=Hks2Hkq(Hks)
    tSet= set(Hkq.args)
    for tt in tSet:
        # get all variables in a term
        ta=tt.free_symbols # obtain var only
        toT=toT.union(ta)
        #print('ta->',ta)
        if( len(ta)>2 ):
            hoT=hoT.union(ta)
    # knowing vadiables in high order term, now construct substitution list
    # the index of ancillary variable should start after number of var in Hk
    nvTO=len(toT) # number of total variables in Hk
    nvHO=len(hoT) # number of variables involved higher order terms in Hk 
    qPair=genSubtitutionPairs(hoT,nvHO)
    # 
    #delta=100
    d=delta; #2*vMax 
    '''
    # do substitution iteratively 
    '''
    H2q=Hkq
    #print(qPair)
    for tx in qPair:
        while (is_ho_equation(H2q)==True):   
            #print(tx)
            # do substitition
            print(tx[0],'*' , tx[1], '->',tx[2] )
            H2q = H2q.subs({tx[0]*tx[1]:tx[2]} ) \
                    + H2sub(tx[0], tx[1],tx[2],d)
    
            H2q=simplify(H2q)
            print('H2q->', H2q)
            print('\n')
    #
    return expand(H2q) 
##
def eq_standardized(Hk,sqChar):
    '''
    change all input variables into a standard variables
    eg. E = x1 + x2*x4  -> E=s1 + s2*s3 
    assume EQU is already symbolic
    '''
    # get a list of all variables in the expression
    # get the symbols ordered  [x0, x1, x2, ...]
    symbSet=tuple(ordered(Hk.free_symbols))
    symOrg=[]
    symStd=[]
    idx=0
    for tx in symbSet:
        symOrg.append(tx)
        if (sqChar.upper()=='S'):
            tsq=symbols('s%d'%idx)
        else:
            tsq=symbols('q%d'%idx)
        symStd.append(tsq)
        # Hamiltonian with standardized variables
        Hk=Hk.subs({tx:tsq})
        '''
        symOrg[x0, x1, ..., xk]
        symStd[s0, s1, ..., sk]
        '''
        #print('idx substitution->',idx)
        idx=idx+1
        #
    return symOrg, symStd, Hk
    

def eqString2Symbolic(EQ):
    '''
    convert string-form expression to symbolic form
    '''
    return(sympify(EQ,evaluate=False))
    
def simplifySquare(Hks):
    '''
    simplify (binary symbolic) function by assigning all
    squared variables with one: si**2 <- 1
    '''
    print('Simplifying si**2 <-1 ')
    # get all free symbols in Hks
    ss=Hks.free_symbols
    # substitute
    for ts in ss:
         Hks= Hks.subs( {ts**2:1}) 
         #print('Processing qi^2->1:',ts)
    return(simplify(Hks))
#
'''
to be replaced by simplify_sq_squares
'''
def simplify_squares(Hkb):
    '''
    simplify (binary symbolic) function by assigning all
    squared variables with one: bi**2 <- 1 
    where bi is a binary variable, either qi={0,1} or si{-1,+1}
    '''
    print('Simplifying bi**2 <-1 ')
    # get all free symbols in Hkb
    bb=Hkb.free_symbols
    # substitute
    for tb in bb:
         Hkb= Hkb.subs( {tb**2:1}) 
         Hkb= Hkb.subs( {tb**3:tb}) 
         #print('Processing qi^2->1:',ts)
    return(simplify(Hkb))

#
def simplify_sq_squares(Hkb, tSq):
    '''
    simplify (binary symbolic) function by assigning all
    squared variables with one: qi**2 <- qi, si**2 <-1  
    where bi is a binary variable, either qi={0,1} or si{-1,+1}
    '''
    print('Simplifying bi**2 <-1 ')
    # get all free symbols in Hkb
    bb=Hkb.free_symbols
    # substitute
    for tb in bb:
        if(tSq.upper()=='S'):
            Hkb= Hkb.subs( {tb**2:1}) 
            Hkb= Hkb.subs( {tb**3:tb}) 
            #print('Processing qi^2->1:',ts)
        else:
            Hkb= Hkb.subs( {tb**2:tb}) 
            Hkb= Hkb.subs( {tb**3:tb}) 
         
    return(simplify(Hkb))
'''
to be replaced by general q/s coefficients
'''
def get_ising_coeffs(H2s):
    '''
    getIsing coefficients
    input:  H2s a 2-body hamiltonian
    output: ising coeffs {hi, Jij}
    since we work in s-domain, assume the symbols are
    s0, s1, ...., s_(NQ-1)
    '''
    NQ=len(H2s.free_symbols)
    hi=np.zeros(NQ)
    Jij=np.zeros((NQ,NQ))
    dc=H2s.as_coefficients_dict()
    # list all symbols: s0, s1, ...
    ss=symbols('s0:%d'%NQ)

    # extract b
    b=dc[1]
    # extract hi
    for m in range(NQ):
        hi[m]=dc[ss[m]];
    # extract Jij
    for m in range(NQ):
        for n in range(m+1,NQ):
            #print('m=',m,'n=',n)
            Jij[m,n]=dc[ss[m]*ss[n]]
    # return the results
    return(b, hi, Jij)
#

def get_qs_coeffs(H2qs,qsType):
    '''
    get Ising or qubo coefficients
    input:  > H2qs  :a 2-body hamiltonian
            > qsType: problem type=> 'q'/0: qubo, 's'/1: ising
    output: qs_coeffs {hi, Jij}
    assume the symbols are
    s0, s1, ...., s_(NQ-1) or q0, q1, ...., q_(NQ-1) 
    '''
    NQ=len(H2qs.free_symbols)
    hi=np.zeros(NQ)
    Jij=np.zeros((NQ,NQ))
    dc=H2qs.as_coefficients_dict()
    # list all symbols: s0, s1, ...
    if(qsType.upper()=='S'):
        qs=symbols('s0:%d'%NQ)
    else:
        qs=symbols('q0:%d'%NQ)

    # extract b
    b=dc[1]
    # extract hi
    for m in range(NQ):
        hi[m]=dc[qs[m]];
    # extract Jij
    for m in range(NQ):
        for n in range(m+1,NQ):
            print('m=',m,'n=',n)
            Jij[m,n]=dc[qs[m]*qs[n]]
    # return the results
    return(b, hi, Jij)
#

