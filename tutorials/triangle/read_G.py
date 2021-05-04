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
    for iw in range(niw):
        #Sigma[spin, iw, :, :] = - numpy.dot(SigmaG[iw, :, :, spin], numpy.linalg.inv(Giw[spin, iw, :, :]))
        Sigma[spin, iw, :, :] = numpy.linalg.inv(G0_iwn[iw, :, spin, :, spin]) -  numpy.linalg.inv(Giw[spin, iw, :, :])

with open("G_omega.txt", "w") as f:
    for iw in range(niw):
        print(iw, Giw[up,iw,0,0].imag, Giw[down,iw,0,0].imag, file=f)
    f.close()

with open("G_omega_all.txt", "w") as f:
    for spin in range(2):
        for i in range(3):
            for j in range(3):
                for iw in range(niw):
                    print(iw, Giw[spin,iw,i,j].real, Giw[spin,iw,i,j].imag, G0_iwn[iw,i,spin,j,spin].real, G0_iwn[iw,i,spin,j,spin].imag, file=f)
                print("", file=f)
                print("", file=f)

with open("Sigma_omega.txt", "w") as f:
    for spin in range(2):
        for i in range(3):
            for j in range(3):
                for iw in range(niw):
                    print(iw, Sigma[spin,iw,i,j].real, Sigma[spin,iw,i,j].imag, file=f)
                print("", file=f)
                print("", file=f)

with open("SigmaG_omega.txt", "w") as f:
    for i in range(3):
        for j in range(3):
            print("#", i, j, file=f)
            for iw in range(niw):
                print(iw, SigmaG[iw,i,j, 0].real, SigmaG[iw,i,j, 0].imag, file=f)
            print("", file=f)

with open("SigmaG_l.txt", "w") as f:
    for i in range(3):
        for j in range(3):
            print("#", i, j, file=f)
            for l in range(nl):
                print(l, SigmaG_l[l,i,j, 0].real, SigmaG_l[l,i,j, 0].imag, SigmaG_l[l,i,j, 1].real, SigmaG_l[l,i,j, 1].imag, file=f)
            print("", file=f)
