import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.optimize import least_squares

#function to generate a FIF IFS from data points
def inter_gen(X,Y,d,inc_a = False):
    '''
    Y,X - n by 1 numpy array
    d - N by 1 numpy array
    '''
    
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
    
#function for the chaos game
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
                
    return Z


def tau(n,m):
    return ((n+m)+1)*(n+m)/2 + m

def tau_inv(k):
    
    def fl(k):
        return np.floor((1+np.sqrt(1+8*k))/2) - 1
    
    i = fl(k)- (k - (fl(k)*(fl(k)+1)/2) )
    j = k - (fl(k)*(fl(k)+1)/2)
    
    return i,j

def add_inv(n):

    ans = np.zeros((n+1,2))
    
    for i in range(len(ans)):
        
        ans[i,0] = i
        ans[i,1] = n - i
        
    return ans

def trinomial ( i, j, k ):

#This function was taken from the internet 
     
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    11 April 2015
#
#  Author:
#
#    John Burkardt
#
    
    i = int(i)
    j = int(j)
    k = int(k)

    from sys import exit

    if ( i < 0 or j < 0 or k < 0 ):
        print ( '' )
        print ( 'TRINOMIAL - Fatal error!' )
        print ( '  Negative factor encountered.' )
        exit ( 'TRINOMIAL - Fatal error!' )

    value = 1

    t = 1

    for l in range ( 1, i + 1 ):
        #   value = value * t // l
        t = t + 1

    for l in range ( 1, j + 1 ):
        value = value * t // l
        t = t + 1

    for l in range ( 1, k + 1 ):
        value = value * t // l
        t = t + 1

    return value



def A_r(a,b,c,d,e,f,N=6,rational= False):
    
    A = np.zeros((N,N))
    
    if rational:
        
         A = sp.zeros(N,N)
    
    
    for u in range(N):
        for v in range(N):

            n,m = tau_inv(u)

            i_i, j_j = tau_inv(v)

            i_i = add_inv(int(i_i))

            j_j = add_inv(int(j_j))

            for i in range(i_i.shape[0]):

                for j in range(j_j.shape[0]):

                    if (i_i[i,0] + j_j[j,0]<=n):
                        
                        if (i_i[i,1] + j_j[j,1]<=m):
                            
                            i = int(i)
                            
                            j = int(j)

                            B = trinomial( i, j, n-i-j )*trinomial( int(i_i[i,1]), int(j_j[j,1]), m-int(i_i[i,1])-int(j_j[j,1]) )
                            
                            
                            C = a**i * c**int(i_i[i,1]) * b**j * d**int(j_j[j,1]) * e**(int(n-i-j)) * f**(int(m -i_i[i,1]-j_j[j,1]))
                                                        
                            A[u,v] += int(B)*C
          
    return A


def Phi(a_l,b_l,c_l,d_l,e_l,f_l,p_l,it = 1,n = 6,rational = False):
    
    assert len(a_l) == len(b_l) and len(a_l) == len(p_l)
    
    mat = np.zeros((n,n))
    
    if rational:
        
        mat = sp.zeros(n,n)
    
    N = len(a_l)
    
    for i in range(N):
        
        mat = p_l[i]*A_r(a_l[i],b_l[i],c_l[i],d_l[i],e_l[i],f_l[i],N=n,rational = rational) + mat
    
    if not rational:
    
        mat = np.linalg.matrix_power(mat,it)
        
    return mat

def moments(a_l,b_l,c_l,d_l,e_l,f_l,p_l,it = 20,n=6,direct = False):
    
    if direct:
        
        tmp = Phi(a_l,b_l,c_l,d_l,e_l,f_l,p_l,it = 1,n=n, rational = True) - sp.eye(n)
        
        tmp = tmp.nullspace()
        
        if len(tmp) != 0:
            
            'Fire!! Nullspace larger than expected'
        
        ans = (tmp[0])/tmp[0][0]
        
        return ans
    
    return Phi(a_l,b_l,c_l,d_l,e_l,f_l,p_l,it = it,n=n)@np.ones(n)

#wrapper function for FIF
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
        
        return moments(a_l,b_l,c_l,d_l,e_l,f_l,p_l,it = 20,n=n,direct = direct)

    return  Phi(a_l,b_l,c_l,d_l,e_l,f_l,p_l,it = 1,n = n,rational = False)
    

def FIF_wrap(x, moment_calc = False):
    
    assert (len(x)+1)%3==0
    
    N = int((len(x)+1)//3)
    
    X = np.zeros(N+1)
    
    Y = np.zeros(N+1)
    
    d = np.zeros(N)
    
    X[-1] = 1
    
    X[1:N] = x[0:N-1]
    
    Y[1:] = x[N-1:2*N-1]
    
    Y[-1] = 0
    
    d = x[2*N-1:]

    A,b,a = inter_gen(X,Y,d,inc_a = True)
    
    
    return mom_wrap(A,b,moment_calc = moment_calc)
    




