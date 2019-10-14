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


chi2=np.loadtxt(r'Lam0/buffer/chi2.dat')
chi4=np.loadtxt(r'Lam0/buffer/chi4.dat')


# Create figure
fig=plt.figure(figsize=(9, 3.5))
ax1=fig.add_subplot(121)

ax1.plot(chi2,color='r',linestyle='--',linewidth=2,markersize=5,label=r'$\chi^B_2$')


ax1.axis([0,300,-0.05,0.2])

ax1.set_xlabel('$T\,[\mathrm{MeV}]$', fontsize=15, color='black')
ax1.set_ylabel(r'$\chi_2$', fontsize=15, color='black')



for label in ax1.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax1.yaxis.get_ticklabels():
    label.set_fontsize(10)


# Plot two
ax2=fig.add_subplot(122)

ax2.plot(chi4,color='k',linestyle='-',linewidth=2,markersize=5,label=r'$\chi^B_4$')

ax2.axis([0,300,-0.15,0.2])

ax2.set_xlabel('$T\,[\mathrm{MeV}]$', fontsize=15, color='black')
ax2.set_ylabel(r'$\chi_4$', fontsize=15, color='black')
ax2.legend(loc=0,fontsize=7.3,frameon=False,shadow=True,handlelength=3.,borderpad=0.5,borderaxespad=1)

for label in ax2.xaxis.get_ticklabels():
    label.set_fontsize(10)
for label in ax2.yaxis.get_ticklabels():
    label.set_fontsize(10)



fig.subplots_adjust(top=0.9, bottom=0.15, left=0.1, right=0.95, hspace=0.35,
                    wspace=0.2)
                   

fig.savefig("chi.pdf")

#plt.show()
