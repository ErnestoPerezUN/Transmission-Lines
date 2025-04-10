# Definition of classes for conductors
import numpy as np
import pandas as pd


class Tower:
    n_cond=1
    name=''
    phases=[]  
    
    def __init__(self, name, cables):
        self.name = name
        self.phases = cables
        self.n_cond=len(self.phases)
class cable:
    name=''
    r:float #[m]
    Rdc:float
    Rac:float
    Tac:float #[degC]
    Imax:float #[A]
    def __init__(self, name='', r=.05, Rdc=0, Rac=0, Tac=75, Imax=100):
        self.name = name
        self.r = r
        self.Rdc = Rdc
        self.Rac = Rac
        self.Tac = Tac
        self.Imax = Imax
        
        
    def __str__(self):
        return self.name
class conductor:
    c:cable
    pos_x:float
    pos_y:float
    def __init__(self, c, phase='', pos_x=0, pos_y=10):
        self.c = c  # Instance of the cable class
        self.phase = phase
        self.pos_x = pos_x
        self.pos_y = pos_y
class impedance:
    def __init__(self,Linf, Z, Z012, Cap, Zinf, Zeq,zp,zm,z0,z1):
        self.Linf = Linf
        self.Z012 = Z012
        self.Cap = Cap
        self.Zinf = Zinf
        self.Zeq = Zeq
        self.zp = zp
        self.zm = zm
        self.z0 = z0
        self.z1 = z1
def calc_param(T1:Tower, freq, sigma_g, mu_gr, eps_gr) ->impedance:
    
    # Physical constants
    C0 = 299_792_458           # Speed of light in vacuum (m/s)
    mu_0 = 4 * np.pi * 1e-7     # Permeability of free space (H/m)
    eps_0 = 1 / (mu_0 * C0**2)   # Permittivity of free space (F/m)
    Pos_cond=np.zeros((T1.n_cond,2))
    rw=np.zeros((T1.n_cond,1))
    for i in range(T1.n_cond):
        Pos_cond[i,0]=T1.phases[i].pos_x
        Pos_cond[i,1]=T1.phases[i].pos_y
        rw[i]=T1.phases[i].c.r
    #Defining Transformation Matrix
    a120=120*np.pi/180
    a=1*np.cos(a120)+1j*np.sin(a120)
    T=np.array([[1, 1, 1],
            [1, a**2, a],
              [1, a, a**2]])
    Ti=np.linalg.inv(T)

    n_cond=T1.n_cond

    if n_cond==len(rw):

        w=2* np.pi*freq 
        #CÁLCULO DE LOS VALORES DE PERMITIVIDAD Y PERMEABILIDAD
        mu_g=mu_0*mu_gr 
        eps_g=eps_0*eps_gr 

        p=np.sqrt(1/(1j*w*mu_g*(sigma_g+1j*w*eps_g))) #Calculo de Distancia Compleja
        I_vec=np.ones((n_cond,1))
        M_X=Pos_cond[:,0]*I_vec
        M_Y=Pos_cond[:,1]*I_vec
        M_Xt=np.transpose(M_X)
        M_Yt=np.transpose(M_Y)
        dij=(((M_X-M_Xt)**2+(M_Y-M_Yt)**2)**0.5)
        dij=dij+np.diag(rw) 
        Dij=(((M_X-M_Xt)**2+(M_Y+M_Yt)**2)**0.5)
        Dij_p=(((M_X-M_Xt)**2+(M_Y+M_Yt+2*p)**2)**0.5)
        #Cálculo de distancias
        Mat_P=np.log(Dij/dij)
        Mat_Pp=np.log(Dij_p/dij)
        Linf=mu_0/(2* np.pi)*Mat_P 
        Z=1j*2* np.pi*freq*mu_0/(2* np.pi)*Mat_Pp 
        Zinf=1j*2* np.pi*freq*Linf 
        Cap=2* np.pi*eps_0*np.linalg.inv(Mat_P) 
    

        #c=2* np.pi*eps_0/np.log((2*h)/rw)
        #z=1j*w*mu_0*(2* np.pi)*np.log((2*h+2*p)/rw)
        A=Z[0:3,0:3]
        B=Z[0:3,3]
        C=Z[3,0:3]
        D=Z[3,3]
        D=np.array(D)
        Dinv=1/D
        DC=Dinv*C
        Zeq=A-B@DC

        Z012=Ti @ Zeq @ Ti.T
        zm=np.average(np.triu(Zeq))
        zp=np.average(np.diag(Zeq))
        z0=zp+2*zm
        z1=zp-zm
        impe=impedance(Linf, Z, Z012, Cap, Zinf, Zeq,zp,zm,z0,z1)
        return impe
        #return Linf, Z, Z012, Cap, Zinf, Zeq,zp,zm,z0,z1
        
    else: 
        print('La matriz de posición debe coincidir con la matriz de conductores')