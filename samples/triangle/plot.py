import numpy
import h5py

nf = 6
norb = int(nf/2)

mu = 0.0

hf = h5py.File('input.out.h5','r')

#Dataset {1000, 3, 3, 2, 2} (iwn, site, site, flavor=spin, re/imag)
SigmaG = hf['/SigmaG_omega'].value[:,:,:,:,0] + 1J * hf['/SigmaG_omega'].value[:,:,:,:,1]
SigmaG_l = hf['/SigmaG_legendre'].value[:,:,:,:,0] + 1J * hf['/SigmaG_legendre'].value[:,:,:,:,1]

niw = SigmaG.shape[0]
nsite = SigmaG.shape[1]
nf = SigmaG.shape[3]
nl = SigmaG_l.shape[0]

G0_iwn = numpy.load('G0_iwn.npy')

Giw = numpy.zeros((2,niw,nsite,nsite), dtype=complex)
up=0
down=1

for spin in range(2):
    for n in range(niw):
        Giw[spin,n,:,:] = G0_iwn[n,:,spin,:,spin] - numpy.dot(G0_iwn[n,:,spin,:,spin], SigmaG[n,:,:,spin])

Sigma = numpy.zeros_like(Giw)
for spin in range(2):
    for iw in xrange(niw):
        #Sigma[spin, iw, :, :] = - numpy.dot(SigmaG[iw, :, :, spin], numpy.linalg.inv(Giw[spin, iw, :, :]))
        Sigma[spin, iw, :, :] = numpy.linalg.inv(G0_iwn[iw, :, spin, :, spin]) -  numpy.linalg.inv(Giw[spin, iw, :, :])

with open("G_omega.txt", "w") as f:
    for iw in xrange(niw):
        print>>f, iw, Giw[up,iw,0,0].imag, Giw[down,iw,0,0].imag
    f.close()

with open("G_omega_all.txt", "w") as f:
    for spin in xrange(2):
        for i in xrange(3):
            for j in xrange(3):
                for iw in xrange(niw):
                    print>>f, iw, Giw[spin,iw,i,j].real, Giw[spin,iw,i,j].imag, G0_iwn[iw,i,spin,j,spin].real, G0_iwn[iw,i,spin,j,spin].imag
                print >>f, ""
                print >>f, ""

with open("Sigma_omega.txt", "w") as f:
    for spin in xrange(2):
        for i in xrange(3):
            for j in xrange(3):
                for iw in xrange(niw):
                    print>>f, iw, Sigma[spin,iw,i,j].real, Sigma[spin,iw,i,j].imag
                print >>f, ""
                print >>f, ""

with open("SigmaG_omega.txt", "w") as f:
    for i in xrange(3):
        for j in xrange(3):
            print>>f, "#", i, j
            for iw in xrange(niw):
                print>>f, iw, SigmaG[iw,i,j, 0].real, SigmaG[iw,i,j, 0].imag
            print>>f, ""

with open("SigmaG_l.txt", "w") as f:
    for i in xrange(3):
        for j in xrange(3):
            print>>f, "#", i, j
            for l in xrange(nl):
                print>>f, l, SigmaG_l[l,i,j, 0].real, SigmaG_l[l,i,j, 0].imag, SigmaG_l[l,i,j, 1].real, SigmaG_l[l,i,j, 1].imag
            print>>f, ""


import matplotlib.pyplot as plt
import pylab
 
params = {'backend': 'pdf',
          'axes.labelsize': 24,
          'axes.titlesize': 24,
          'legend.fontsize': 24,
          'xtick.labelsize': 24,
          'ytick.labelsize': 24,
          'text.usetex': True,
          'font.family': 'serif',
          }
pylab.rcParams.update(params)

beta = 10.0
wn = numpy.array([(2*iw+1)*numpy.pi/beta for iw in xrange(niw)])

plt.figure(1)
plt.xlim([0, 10])
plt.xlabel("$\omega_n$")
plt.ylabel("$\mathrm{Re}~\Sigma_{\sigma,ij}(i\omega_n)$")
for spin in xrange(2):
   for i in xrange(3):
      for j in xrange(3):
          if i == j:
              plt.plot(wn, Sigma[spin, :, i, j].real, ls='-', marker='o', color='r', markersize=10)
          else:
              plt.plot(wn, Sigma[spin, :, i, j].real, ls='-', marker='x', color='b', markersize=10)
plt.tight_layout()
plt.savefig("Sigma-Re.pdf")
plt.close()

plt.figure(1)
plt.xlim([0, 10])
plt.xlabel("$\omega_n$")
plt.ylabel("$\mathrm{Im}~\Sigma_{\sigma,ij}(i\omega_n)$")
for spin in xrange(2):
   for i in xrange(3):
      for j in xrange(3):
          if i == j:
              plt.plot(wn, Sigma[spin, :, i, j].imag, ls='-', marker='o', color='r', markersize=10)
          else:
              plt.plot(wn, Sigma[spin, :, i, j].imag, ls='-', marker='x', color='b', markersize=10)
plt.tight_layout()
plt.savefig("Sigma-Im.pdf")
plt.close()
