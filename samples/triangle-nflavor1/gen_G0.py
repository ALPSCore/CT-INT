import sys
import numpy as np
from scipy.integrate import simps
from math import pi
from SOI_lib import *
from itertools import *

def ft_to_tau_hyb(ndiv_tau, beta, matsubara_freq, tau, Vek, data_n, data_tau, cutoff):
 for it in range(ndiv_tau+1):
     tau_tmp=tau[it]
     if it==0:
         tau_tmp=1E-4*(beta/ndiv_tau)
     if it==ndiv_tau:
         tau_tmp=beta-1E-4*(beta/ndiv_tau)
     ztmp=0.0
     for im in range(cutoff):
         ztmp+=(data_n[im]+1J*Vek/matsubara_freq[im])*np.exp(-1J*matsubara_freq[im]*tau_tmp)
     ztmp=ztmp/beta
     data_tau[it]=2.0*ztmp.real-0.5*Vek

ndiv_dos = 10000
W=2.0
e_smp=np.linspace(-W,W,ndiv_dos)
dos_smp=np.sqrt(W**2-e_smp**2)/(0.5*pi*W**2)
ek_var = simps(dos_smp*(e_smp**2),e_smp)
print "ek_var", ek_var

vbeta=10.
mu = 0.0 #please set mu=0.0 for half filling (U/2)

ndiv_tau=1000
cutoff = ndiv_tau/2
nf=6
norb=nf/2

matsubara_freq=np.zeros((ndiv_tau,),dtype=float)
for im in range(ndiv_tau):
    matsubara_freq[im]=((2*im+1)*np.pi)/vbeta

tau=np.zeros((ndiv_tau+1,),dtype=float)
for it in range(ndiv_tau+1):
    tau[it]=(vbeta/ndiv_tau)*it

Himp = np.zeros((norb,2,norb,2),dtype=float)
for i in xrange(norb):
    for j in xrange(norb):
        if i!=j:
            Himp[i,0,j,0] = -1
            Himp[i,1,j,1] = -1
Himp = Himp.reshape((nf,nf))

#print(Himp)
#evals,evecs = np.linalg.eigh(Himp)
#print evals

G0_iwn = np.zeros((ndiv_tau,nf,nf),dtype=complex)
G0_tau =np.zeros((ndiv_tau+1,nf,nf),dtype=complex)

#Compute G0
for im in range(ndiv_tau):
    f_tmp = dos_smp/(1J*matsubara_freq[im]-e_smp)
    delta_im = simps(f_tmp, e_smp)
    G0_iwn[im, :, :] = np.linalg.inv((1J*matsubara_freq[im]+mu)*np.identity(nf,dtype=complex) - Himp - delta_im * np.identity(nf))

ft_to_tau_hyb(ndiv_tau, vbeta, matsubara_freq, tau, np.identity(nf), G0_iwn, G0_tau, ndiv_tau)

def cut_small_value(v, eps=1e-12):
    if np.abs(v)>eps:
        return v
    else:
        return 0.0

G0_tau = G0_tau.reshape((ndiv_tau+1,norb,2,norb,2))
G0_iwn = G0_iwn.reshape((ndiv_tau,norb,2,norb,2))

f = open('G0_TAU.txt','w')
print>>f, 1, 2*norb, ndiv_tau+1, vbeta

for iorb, ispin, jorb, jspin, itau in product(range(norb), range(2), range(norb), range(2), range(ndiv_tau+1)):
    if ispin == jspin:
         print>>f, 0, 2*iorb+ispin, 2*jorb+jspin, itau, cut_small_value(G0_tau[itau,iorb,0,jorb,0].real), cut_small_value(G0_tau[itau,iorb,0,jorb,0].imag)
    else:
         print>>f, 0, 2*iorb+ispin, 2*jorb+jspin, itau, 0.0, 0.0
f.close()


np.save('G0_iwn.npy', G0_iwn)
