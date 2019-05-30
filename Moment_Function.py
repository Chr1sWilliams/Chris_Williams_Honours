import numpy as np
import sympy as sp


#tau function as in text
def tau(n,m):
    return ((n+m)+1)*(n+m)/2 + m

#inverse tau function 
def tau_inv(k):
    
    def fl(k):
        return np.floor((1+np.sqrt(1+8*k))/2) - 1
    
    i = fl(k)- (k - (fl(k)*(fl(k)+1)/2) )
    j = k - (fl(k)*(fl(k)+1)/2)
    
    return i,j

#inverse for add function
def add_inv(n):

    ans = np.zeros((n+1,2))
    
    for i in range(len(ans)):
        
        ans[i,0] = i
        ans[i,1] = n - i
        
    return ans


#compute the trinomial cp-efficients quickly 
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


#construct the matrix `phi'
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

#make the operator Phi
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

#function that computes the moments
def moments(a_l,b_l,c_l,d_l,e_l,f_l,p_l,it = 20,n=6,direct = False):
    
    if direct:
        
        tmp = Phi(a_l,b_l,c_l,d_l,e_l,f_l,p_l,it = 1,n=n, rational = True) - sp.eye(n)
        
        tmp = tmp.nullspace()
        
        if len(tmp) != 0:
            
            'Fire!! Nullspace larger than expected'
        
        ans = (tmp[0])/tmp[0][0]
        
        return ans
    
    return Phi(a_l,b_l,c_l,d_l,e_l,f_l,p_l,it = it,n=n)@np.ones(n)

#function for eltons theorem
def elt(X,size = 21):
    
    ans = np.zeros(size)
    
    for k in range(size):
        
        i , j = tau_inv(k)
        
        ans[k] = np.sum(np.power(X[:,0],i)*np.power(X[:,1],j),axis = 0)/len(X)
    
    return ans
