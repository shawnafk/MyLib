from functools import wraps
import numpy as np
import cmath
import matplotlib.pyplot as plt
from sympic_io.tools import unit as unit
#notes
#dr need to be classified
#omega k etc parameters should in a unify naming
#using function to do repeated works


class Plasma:
    def __init__(self, m, q, U):
        self.U = U
        self.m = m
        self.q = q
        #ref this works under wci much less than wpi, i.e. high density low magnetic field
        #full  exp
        #wlh=(omega_pi**2/(1+omega_pe**2/omega_ce**2) + omega_ci**2)**0.5

    def gen_w(self, n, B):
        w_ce = unit.Oce(B, self.U)
        w_pe = unit.Ope(n, self.U)
        w_ci = unit.Oci(self.m, self.q, B, self.U)
        w_pi = unit.Opi(self.m, self.q, n, self.U) 
        w_gmg = np.sqrt(w_ci*w_ce)
        #resonance
        w_lh = w_pi/(1+w_pe**2/w_ce**2)**0.5
        w_hh = (w_ce**2 + w_pe**2)**0.5
        #cut-off
        w_p = (w_pi**2 + w_pe**2)**0.5
        #w_Xm
        w_L = (w_ci - w_ce )/2 + ((w_ci+w_ce)**2/4 + w_p**2)**0.5
        #w_Xp
        w_R = (w_ce - w_ci )/2 + ((w_ce+w_ci)**2/4 + w_p**2)**0.5
        return (w_ce, w_pe, w_ci, w_pi, w_gmg, w_lh,w_hh,w_L,w_R)

    def k2w(self, n, B, k, ct):
        (wce, wpe, wci, wpi, wgmg, wlh,whh,wL,wR) = self.gen_w(n, B)
        wp = (wpe**2 + wpi**2)**0.5
        m = self.m
        q = self.q
        c = 1
        c10 = 1
        c8 = 2*c**2 * k**2 + wce**2 + wci**2 + 3*wp**2

        c6 = c**4 * k**4 + (2 * c**2 * k**2 + wp**2) * \
                            (wce**2 + wci**2 + 2*wp**2) + (wp**2 + wce*wci)**2

        c4 = c**4 * k**4 * (wce**2 + wci**2 + wp**2) + 2*c**2 * k**2 * (wp**2 + wci*wce)**2 + \
            c**2 * k**2 * wp**2 * (wce**2 + wci**2 - wci*wce) * \
                                   (1+ct**2) + wp**2 * (wp**2+wci*wce)**2

        c2 = c**4 * k**4 * (wp**2 * (wce**2 + wci**2 - wci*wce)*ct**2 + wce*wci*(wp**2 + wci*wce)) +\
            c**2 * k**2 * wp**2*wci*wce*(wp**2 + wci*wce)*(1+ct**2)

        c0 = c**4 * k**4 * wce**2 * wci**2 * wp**2 * ct**2
        w = np.sqrt(np.roots([c10, -c8, c6, -c4, c2, -c0]))
        return w

    def pomega(self, n, B):
        (wce, wpe, wci, wpi, wgmg, wlh, whh, wL,wR) = self.gen_w(n, B)
        print('w_ce={:f}'.format(wce))
        print('w_pe={:f}'.format(wpe))
        print('w_ci={:f}'.format(wci))
        print('w_pi={:f}'.format(wpi))
        print('w_gmg={:f}'.format(wgmg))
        print('w_lh={:f}'.format(wlh))
        print('w_hh={:f}'.format(whh))
        print('w_l={:f}'.format(wL))
        print('w_r={:f}'.format(wR))
    def omegakv_hot(self,n,B,k_para,omega,vte,vti):
        c=1
        w_ce =  unit.Oce(B,self.U)
        w_pe =  unit.Ope(n,self.U)
        w_pi  = unit.Opi(self.m,self.q,n,self.U)
        w_ci  = unit.Oci(self.m,self.q,B,self.U)
        kv=[]
        n_para=k_para*c/omega
        epsi_perp=1+(w_pe/w_ce)**2 - (w_pi/omega)**2
        epsi_para=1-(w_pe/omega)**2 - (w_pi/omega)**2
        epsi_xy=(w_pe)**2 / (w_ce*omega)
        P0=epsi_para*( (n_para**2-epsi_perp)**2 - epsi_xy**2)
        P2=(epsi_perp + epsi_para)*(n_para**2 - epsi_perp) + epsi_xy**2
        P4=epsi_perp
        #C = - epsi_para/epsi_perp *omega/k_para/1000
        P6i= - (3*w_pi**2/omega**2*vti**2/c**2)
        P6e= - (0.75*(w_pe*omega/w_ce**2)**2*vte**2/c**2)
        P6 = P6i + P6e
        for _ in np.roots([P6,P4,P2,P0]):
            if(np.isreal(_) and np.real(_) > 0):
                kv.append(np.sqrt(_)*omega/c)
        return kv
       #WR WL WUH
    def general_kv(self, n, B, kp, w):
        (wce, wpe, wci, wpi, wgmg, wlh,whh,w_L,w_R) = self.gen_w(n, B)
        kv = []
        c = 1
        R = 1 - (wpe**2/w/(w-wce) + wpi**2/w/(w+wci))
        L = 1 - (wpe**2/w/(w+wce) + wpi**2/w/(w-wci))
        P = 1 - (wpe**2/w**2 + wpi**2/w**2)
        S = 0.5 * (R + L)
        D = 0.5 * (R - L)
        A = S
        B = R*L + P*S - P * (kp*c/w)**2 - S * (kp*c/w)**2
        C = P * (R*L - 2 * S * (kp*c/w)**2 + (kp*c/w)**4)
        #P6=-( 3 * (omega_pi*vti/omega)**2 + 3./4. * (omega_pe*omega*vte/omega_ce/omega_ce)**2)
        #print(P6)
        for _ in np.roots([A, -B, C]):
            if(np.isreal(_) and _ > 0):
                kv.append(np.real(cmath.sqrt(_))*w/c)
        return kv

    def general_k(self, n, B, theta, w):
        (wce, wpe, wci, wpi, wgmg, wlh,whh,w_L,w_R) = self.gen_w(n, B)
        k = []
        c = 1
        R = 1 - (wpe**2/w/(w-wce) + wpi**2/w/(w+wci))
        L = 1 - (wpe**2/w/(w+wce) + wpi**2/w/(w-wci))
        P = 1 - (wpe**2/w**2 + wpi**2/w**2)
        S = 0.5 * (R + L)
        D = 0.5 * (R - L)
        sint=np.sin(theta)
        cost=np.cos(theta)
        A = S *sint**2 + P * cost**2
        B = R*L *sint**2 + P*S *(1+cost**2)
        C = P*R*L
        #P6=-( 3 * (omega_pi*vti/omega)**2 + 3./4. * (omega_pe*omega*vte/omega_ce/omega_ce)**2)
        #print(P6)
        for _ in np.roots([A, -B, C]):
            if(np.isreal(_) and _ > 0):
                k.append(np.real(cmath.sqrt(_))*w/c)
        return k

    def gdr(case,n,B,kp,w_range):
        w=[]
        kv=[]
        for __w in w_range:
            for __kv in case.general_kv(n,B,kp,__w):
                if np.isreal(__kv):
                    w.append(__w)
                    kv.append(__kv)
        #zip and sort
        DP=zip(kv,w)
        res=sorted(DP,key= lambda v:v[0])
        kv,w=zip(*res)
        return np.array((kv,w))

    #scan usually we treat kp fixed
    def scan_density(case,n_range,B,kp,w):
        n=[]
        kv=[]
        for __n in n_range:
            for __kv in case.general_kv(__n,B,kp,w):
                n.append(__n)
                kv.append(np.real(__kv))
        DP=zip(n,kv)
        res=sorted(DP,key= lambda v:v[1])
        n,kv=zip(*res)
        return (n,kv)

    def scan_magnet(case,n,B_range,kp,w):
        b=[]
        kv=[]
        for __b in B_range:
            for __kv in case.general_kv(n,__b,kp,w):
                b.append(__b)
                kv.append(np.real(__kv))
        DP=zip(b,kv)
        res=sorted(DP,key= lambda v:v[1])
        b,kv=zip(*res)
        return (b,kv)
    
    #special case ox lr general k can not distinguish them readily
    #wave parallel to B
    def LR_wave(self, n, B, w):
        (wce, wpe, wci, wpi, wgmg, wlh,whh,w_L,w_R) = self.gen_w(n, B)
        R = 1 - (wpe**2/w/(w-wce) + wpi**2/w/(w+wci))
        L = 1 - (wpe**2/w/(w+wce) + wpi**2/w/(w-wci))
        P = 1 - (wpe**2/w**2 + wpi**2/w**2)
        kr = cmath.sqrt(R)*w
        kl = cmath.sqrt(L)*w
        return (kl, kr)

    def OX_wave(self, n, B, w):
        (wce, wpe, wci, wpi, wgmg, wlh,whh,w_L,w_R) = self.gen_w(n, B)
        R = 1 - (wpe**2/w/(w-wce) + wpi**2/w/(w+wci))
        L = 1 - (wpe**2/w/(w+wce) + wpi**2/w/(w-wci))
        S = 0.5 * (R + L)
        P = 1 - (wpe**2/w**2 + wpi**2/w**2)
        kx = cmath.sqrt(R*L/S)*w
        ko = cmath.sqrt(P)*w
        return (ko, kx)
    
    def oxdr(case, n, B, w_range):
        wo = []
        wx = []
        ko = []
        kx = []
        for __w in w_range:
            __k = case.OX_wave(n, B, __w)
            if np.isreal(__k[0]):
                wo.append(__w)
                ko.append(np.real(__k[0]))
            if np.isreal(__k[1]):
                wx.append(__w)
                kx.append(np.real(__k[1]))
        DPo = np.array((ko, wo))
        DPx = np.array((kx, wx))
        return DPo,DPx
    def lrdr(case, n, B, w_range):
        wl = []
        wr = []
        kl = []
        kr = []
        for __w in w_range:
            __k = case.LR_wave(n, B, __w)
            if np.isreal(__k[0]):
                wl.append(__w)
                kl.append(np.real(__k[0]))
            if np.isreal(__k[1]):
                wr.append(__w)
                kr.append(np.real(__k[1]))
        DPl = np.array((kl, wl))
        DPr = np.array((kr, wr))
        return DPl,DPr

#    def gamma_col(self,k11,kL,wr,nue,nui,wpe,wpi,wce):
#        K=(k11**2 + kL**2)**0.5
#        gamma=( ( (wpe/wr)**2 * (k11/K)**2 + (wpe/wce)**2 * (kL/K)**2 ) *nue + (wpi/wr)**2*nui)/((wpe/wr)**2 * (k11/K)**2 + (wpi/wr)**2)
#        return gamma

    #return vertical wave length with given karallel
    #from plot.py, general functions

#    def wz(self,n,B,kp,omega):
#        w_ce=self.w_ce
#        w_ci=self.w_ci
#        w_pe=self.w_pe
#        w_pi=self.w_pi
#        ne=self.n
#        ni=ne/self.q
#        kxx=1-(wpi/w)**2 + (wpe/wce)**2
#        kzz=1-(wpe/w)**2;
#        return 2*np.pi/(kp*np.sqrt(abs(kzz/kxx)))
