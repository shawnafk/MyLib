import numpy as np
import matplotlib.pyplot as plt
#等额本息
class equal_loan_payment:
    def __init__(self,P,R,N):
        self.P = P
        self.R = R
        self.N = N

    def BX(self):
        P=self.P
        R=self.R
        N=self.N
        return P*(R*(1+R)**N)/((1+R)**N-1)
    def BandX(self,n):
        P=self.P
        R=self.R
        N=self.N
        bx=self.BX()
        b = P * R*(1+R)**(n-1)/((1+R)**N-1)
        x = bx - b
        return np.array((b,x))
    #t=equal_loan_payment()
    #total = t.BX()*t.N
'''
     def BandX_all(self):
        P=self.P
        R=self.R
        N=self.N
        bx=self.BX()
        b=np.zeros(N)
        x=np.zeros(N)
        q=np.zeros(N)
        b[-1] = bx
        x[-1] = 0
        q[-1] = 0
        for n in np.arange(start=N-1,stop=0,step=-1):
            q[n-1] = b[n] + q[n]
            x[n-1] = q[n-1] * R
            b[n-1] = bx - x[n-1]
        return b,x,q
    
    def BandX_rec(self,n):
        P=self.P
        R=self.R
        N=self.N
        bx=self.BX()
        if(N==n):
            #b and x
            q = 0
            x = q * R
            b = bx - x 
            return np.array((b,x,q))
        else:
            q = BandX_rec(P,R,N,n+1)[0] + BandX_rec(P,R,N,n+1)[-1]
            x = q * R
            b = bx - x 
            return np.array((b,x,q))
'''  
    
   
 
#等额本金
class equal_principal_payment:
    def __init__(self,P,R,N):
        self.P = P
        self.R = R
        self.N = N
    def BandX(self,n):
        P=self.P
        R=self.R
        N=self.N
        
        b = np.ones(n.shape)*P/N
        x = (P - b*(n-1) )*R
        return np.array((b,x))

Ps=183e4
Rs=5.1/100/12

Pg=120e4
Rg=3.25/100/12
    
N=300

#epp = equal_principal_payment(P,R,N)
elp_s = equal_principal_payment(Ps,Rs,N)
elp_g = equal_principal_payment(Pg,Rg,N)

m = np.arange(1,301)
b_s,x_s=elp_s.BandX(m)
b_g,x_g=elp_g.BandX(m)

