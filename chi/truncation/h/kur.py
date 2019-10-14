#!/usr/bin/env python
# -*- coding: utf-8 -*-
# sphinx_gallery_thumbnail_number = 3

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullFormatter  # useful for `logit` scale
import matplotlib.ticker as ticker
import matplotlib as mpl

mpl.style.use('classic')


# Data for plotting
#chi2lam700=np.loadtxt(r'Lam700/buffer/chi2.dat')
#chi4lam700=np.loadtxt(r'Lam700/buffer/chi4.dat')
#r42lam700=chi4lam700 / chi2lam700

#chi2lam600=np.loadtxt(r'Lam600/buffer/chi2.dat')
#chi4lam600=np.loadtxt(r'Lam600/buffer/chi4.dat')
#r42lam600=chi4lam600 / chi2lam600

#chi2lam500=np.loadtxt(r'Lam500/buffer/chi2.dat')
#chi4lam500=np.loadtxt(r'Lam500/buffer/chi4.dat')
#r42lam500=chi4lam500 / chi2lam500

#chi2lam400=np.loadtxt(r'Lam400/buffer/chi2.dat')
#chi4lam400=np.loadtxt(r'Lam400/buffer/chi4.dat')
#r42lam400=chi4lam400 / chi2lam400

#chi2lam300=np.loadtxt(r'Lam300/buffer/chi2.dat')
#chi4lam300=np.loadtxt(r'Lam300/buffer/chi4.dat')
#r42lam300=chi4lam300 / chi2lam300

chi2lam0=np.loadtxt(r'Lam0/buffer/chi2.dat')
chi4lam0=np.loadtxt(r'Lam0/buffer/chi4.dat')
r42lam0=chi4lam0 / chi2lam0


# Create figure
fig=plt.figure(figsize=(4.5,3.5))
#fig=plt.figure()
ax1=fig.add_subplot(111)

ax1.plot(r42lam0,'olive',linewidth=1,markersize=5,label=r'$0$')


ax1.axis([100,250,-70,70])
#ax1.set_xscale('log')

ax1.set_xlabel('$T[\mathrm{MeV}]$', fontsize=14, color='black')
ax1.set_ylabel(r'$\chi^B_4/\chi^B_2$', fontsize=14, color='black')

ax1.legend(loc=0,fontsize='x-small',frameon=False,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1)

for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)


fig.subplots_adjust(top=0.9, bottom=0.15, left=0.15, right=0.95, hspace=0.35,
                    wspace=0.35)


fig.savefig("kur.pdf")

#plt.show()
