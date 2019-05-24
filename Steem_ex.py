import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
import time
from PIL import Image
import random
from scipy import optimize
import Moment_Function as mf



def mom_wrap(x,D,it = 20,n=5):
    
    assert len(x)%6==0
    
    N = len(x)
    
    a_l = x[:N//6]
    
    b_l = x[1*(N//6):2*(N//6)]
    
    c_l = x[2*(N//6):3*(N//6)]
    
    d_l = x[3*(N//6):4*(N//6)]
    
    e_l = x[4*(N//6):5*(N//6)]
    
    f_l = x[5*(N//6):]
    
    p_l =  np.power(np.abs(a_l),D-1)*np.abs(d_l)/(np.sum(np.power(np.abs(a_l),D-1)*np.abs(d_l)))
    
    return moments(a_l,b_l,c_l,d_l,e_l,f_l,p_l,it = it,n=n)

def dim(a_l):
    
    f = lambda x : np.sum(np.power(np.abs(a_l),x))-1
    
    ans = optimize.root(f,x0 = 0.5)
    
    return ans.x

def convertIFS(IFS_mat):
    N = len(IFS_mat)
    
    ans = np.zeros(N*6)
    
    r = 0
    for j in range(6):
        for i in range(N):
            ans[r]=IFS_mat[i][j]
            r = r+1
            
    return ans


def mom_wrap_Steem(x,D,it = 20,n=5,sim = False):

    
    assert len(x)%6==0
    
    N = len(x)
    
    a_l = x[:N//6]
    
    b_l = x[1*(N//6):2*(N//6)]
    
    c_l = x[2*(N//6):3*(N//6)]
    
    d_l = x[3*(N//6):4*(N//6)]
    
    #project to similitude
    
    if sim:
    
        ad_l = (np.abs(a_l)+np.abs(d_l))/2

        a_l = np.sign(a_l)*ad_l

        d_l = np.sign(d_l)*ad_l

        bc_l = (np.abs(b_l)+np.abs(c_l))/2

        b_l = np.sign(b_l)*bc_l

        c_l = np.sign(c_l)*bc_l
    
    
    e_l = x[4*(N//6):5*(N//6)]
    
    f_l = x[5*(N//6):]
    
    
    prob = np.power(np.sqrt(np.abs(a_l*d_l - b_l*c_l)),D-1)
    
    p_l =  prob/(np.sum(prob))
    
    print('Probabilities:')
    print(p_l)
    

    return mf.moments(a_l,b_l,c_l,d_l,e_l,f_l,p_l,it = it,n=n)


def chaos(IFS,prob,it = 10000,imgxy = 480):
    
    
    N = len(prob)
    
    image = Image.new("RGB", (imgxy, imgxy),"white")
    
    xa = -0.2
    xb = 1.2
    ya = -0.2
    yb = 1.2

    #starting values
    x=np.random.rand()#0.5
    y=np.random.rand()#0.5 
    
    X = np.zeros((it,2))
    
    for i in range(it):
        
        p = random.random()
        
        P = prob[0]
        
        for j in range(N):
            
            if p < P:
                
                
                x0 = IFS[j]*x + IFS[(N) + j]*y + IFS[4*(N) + j]
                
                y = IFS[2*(N) + j]*x + IFS[3*(N) + j]*y + IFS[5*(N) + j]
                
                x = x0
                
                X[i,0] = x
                X[i,1] = y
                
                
                if x>xa and x<xb and y>ya and y<yb:
                    
                    jx = int((x - xa) / (xb - xa) * (imgxy - 1)) 
                    jy = (imgxy - 1) - int((y - ya) / (yb - ya) * (imgxy - 1))
                
                    if j == 0:
                        image.putpixel((jx, jy), (255,0,0,255))
                    elif j == 1:
                        image.putpixel((jx, jy), (255,255,0,255))
                    else: 
                        image.putpixel((jx, jy), (0,0,255,255))
                
                
                break
                
            else: 
                P = P + prob[j+1]

        
    return image, X
                
    
     
def Phi_wrap(x,D,it = 1,n=6, sim = False ):    
    assert len(x)%6==0
    
    N = len(x)
    
    a_l = x[:N//6]
    
    b_l = x[1*(N//6):2*(N//6)]
    
    c_l = x[2*(N//6):3*(N//6)]
    
    d_l = x[3*(N//6):4*(N//6)]
    
    e_l = x[4*(N//6):5*(N//6)]
    
    f_l = x[5*(N//6):]
    
    if sim:
    
        ad_l = (np.abs(a_l)+np.abs(d_l))/2

        a_l = np.sign(a_l)*ad_l

        d_l = np.sign(d_l)*ad_l

        bc_l = (np.abs(b_l)+np.abs(c_l))/2

        b_l = np.sign(b_l)*bc_l

        c_l = np.sign(c_l)*bc_l
    
    prob = np.power(np.sqrt(np.abs(a_l*d_l - b_l*c_l)),D-1)
    
    p_l =  prob/(np.sum(prob))               
    

    global k
    print(k)
    if k % 15 == 0:
        print('Current image:')
        img, _ = chaos(x,p_l,it = 20000)
        plt.imshow(img)
        plt.show()


    k = k + 1
    
    return mf.Phi(a_l,b_l,c_l,d_l,e_l,f_l,p_l,it = it,n=n)


