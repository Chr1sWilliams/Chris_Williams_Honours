import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.optimize import least_squares
import Moment_Function as MF

def inter_gen(X,Y,d,inc_a = False):
    
    assert(X.shape == Y.shape)
    
    n = X.shape[0]
    
    N = d.shape[0]
    
    assert(n-1 == N)
    
    A = []
    
    b = []
    
    a = []
    
    for i in range(1,n):
        
        A.append(np.array([[ X[i] - X[i-1], 0  ], [Y[i] - Y[i-1] - d[i-1]*Y[-1], d[i-1]]  ] ,dtype=float))
        
        b.append(np.array([  X[i-1]  , Y[i-1]  ],dtype=float).reshape((2,1)))
        
        a.append(X[i] - X[i-1])
        
    if inc_a:
        
        return A,b,a
        
    return A,b
    

def mom_wrap(A,b,it = 20,n=21,direct = False, moment_calc = False):
    
    assert len(A) == len(b)
    
    N = len(A)
    
    a_l = []
    
    b_l = []
    
    c_l = []
    
    d_l = []
    
    e_l = []
    
    f_l = []
    
    for i in range(N):
        
        Atmp = A[i]
        btmp = b[i]
        
        a_l.append(Atmp[0,0])

        b_l.append(Atmp[0,1])

        c_l.append(Atmp[1,0])

        d_l.append(Atmp[1,1])

        e_l.append(btmp[0])

        f_l.append(btmp[1])
    

    p_l =  a_l
    
    if moment_calc:
        
        return MF.moments(a_l,b_l,c_l,d_l,e_l,f_l,p_l,it = it,n=n,direct = direct)

    return  MF.Phi(a_l,b_l,c_l,d_l,e_l,f_l,p_l,it = 1,n = n,rational = False)


def FIF_wrap(x, maps = False,moment_calc = False):
    
    assert (len(x)-2)%3==0
    
    N = int((len(x)-2)//3)
    
    X = x[0:N+1]
    
    Y = x[N+1:2*N+2]
    
    d = x[2*N+2:]
    
    A,b,a = inter_gen(X,Y,d,inc_a = True)
    
    if maps:
        return A,b,a
    
    return mom_wrap(A,b,moment_calc = moment_calc,n=36), a

def find_nearest( value, array):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return int(idx)


def find_nearest_v( value_v, array):
    
    n = len(value_v)
        
    ans = []
        
    for i in range(n):
        
        ans.append(int(find_nearest( value_v[i], array = array)))
        
    return ans

def chaos(A,b,p = None ,it =1000,burn = 5):
    
    assert(len(A) == len(b))
    
    N = len(A)
    
    if p == None:
        p = np.ones(N)/N
    
    
    z = np.random.rand(2).reshape([2,1])
    
    Z = np.zeros((it,2))
    
    for i in range(it):

        prob = np.random.rand()
        
        P = 0
        
        for q in range(N):
            
            P += p[q]
            
            if prob < P:
                
                
                z = A[q]@z + b[q]
                
                if i>burn:
                    
                    Z[i,:] = z.flatten()
                
                break
    Z = Z[Z[:,0].argsort()]            
    return Z
    


