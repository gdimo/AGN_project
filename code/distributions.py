#This file will generate distribution figures
#Bins for the whole sample has been extracted through
#Kaplan Meier Estimator using ASURV package


#For the KM estimator
#7 bins have been used
#       parameter   |   initial value   |   Range   |   max value   |   
#       ______________________________________________________________
#           KT      |       8           | 8 +(7*24) |       170     |
#       ____________|_______________________________|_______________|_
#           tau     |       0.20        |0.2+(7*1.2)|       8       |
#       ______________________________________________________________
#
#INSERT LIBRARIES AND TABLES
from astropy.table import Table, Column, vstack, MaskedColumn

import matplotlib
import matplotlib.ticker as mtick
import matplotlib.pyplot as plt
import numpy as np

all_tbl = Table.read('./data/all_table.csv', format='csv')

all_tbl['KT_high'].fill_value = -99.9 
all_tbl['KT_low'].fill_value = -99.9

all_tbl = all_tbl.filled()

msk_kt = all_tbl['KT_high'] < 0
lo_tbl = all_tbl[msk_kt]
go_tbl = all_tbl[~msk_kt]

#Bins from KM estimator
# KT bins

bval = [14.728,30.599,14.585,6.336,5.188,3.113,12.450]
bval_n = [0.17,0.35,0.17,0.07,0.06,0.036,0.14] #normalized values
bcenter = [20.0,44.0,68.0,92.0,116.0,140.0,164.0]

#tau bins


tval = [30.783,24.011,15.945,3.512,4.264,0.0,3.0]
tval_n = [0.35,0.28,0.18,0.04,0.05,0.0,0.035] #normalized values
tcenter = [0.76,1.88,3.00,4.12,5.24,6.36,7.48]


#------------------------------------------------------------------------------
#---------------------------------------------------------------------------------
#--------------TWO SUBPLOTS--------------------------------------------------------
#-----------------------------------------------------------------------------

fig, [ax1, ax2] = plt.subplots(2,1,figsize=(6,10))
#plt.rcParams['xtick.top'] = True
#plt.rcParams['ytick.right'] = True

ax1.hist(go_tbl['KT'],bins=7,weights=np.ones_like(go_tbl['KT']/len(go_tbl['KT'])),
        density=True,histtype='step',color='b',
        alpha=1,linewidth=2, edgecolor='b', label=r'Constrained $kT_e$')

ax1.hist(lo_tbl['KT'],bins=7,weights=np.ones_like(lo_tbl['KT']/len(lo_tbl['KT'])),histtype='step',color='k',
        density=True,
        alpha=1,linewidth=2,linestyle='dotted', edgecolor='black', label='Low limit fixed values')

ax1.tick_params(axis='both', which='major', labelsize=10, direction='in')
#ax1.set_xticks(fontsize=15)
#ax1.set_yticks([5.0,10.0,15.0,20.0,25.0,30.0],[5,10,15,20,25,30])
#ax1.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1,decimals=0))



ax1.legend()

ax2.bar(bcenter,bval,width=24,
            color='none',alpha=1,edgecolor='black',linewidth=2,linestyle='dashed',
            label='KM estimator')

ax2.hist(go_tbl['KT'],bins=7,histtype='step',color='b',
        alpha=1,linewidth=2, edgecolor='b', label=r'Constrained $kT_e$')

ax2.legend()
#plt.grid(b=True, which='major', ls=':')
plt.xlabel(r'Plasma Temperature $kT_e$',fontsize=15)
ax2.tick_params(axis='both', which='major', direction='in')
#ax2.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1,decimals=0))

plt.xticks(fontsize=10, rotation=0)
#plt.yticks([5.0,10.0,15.0,20.0,25.0,30.0],[5,10,15,20,25,30],fontsize=10, rotation=0)
plt.tight_layout()


plt.savefig('./figures/kt_dist.pdf',dpi=400)

#----------------OPTICAL DEPTH TWO FIGS-----------------------------------------------



#-------------------------------------------------------
fig, [ax1, ax2] = plt.subplots(2,1,figsize=(6,10))
#plt.rcParams['xtick.top'] = True
#plt.rcParams['ytick.right'] = True

ax1.hist(go_tbl['Tau'],bins=7,histtype='step',color='b',#density=1,
        alpha=1,linewidth=2, edgecolor='b', label=r'Constrained $kT_e$')

ax1.hist(lo_tbl['Tau'],bins=7,histtype='step',color='k',#density=1,
        alpha=1,linewidth=2, linestyle='dotted' ,edgecolor='black', label='Low limit fixed values')

ax1.tick_params(axis='both', which='major', labelsize=10, direction='in')
#ax1.set_xticks(fontsize=15)
ax1.set_yticks([0.0,0.05,0.10,0.15,0.20,0.25,0.30,0.35],[0,5,10,15,20,25,30,35])
#ax1.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1,decimals=0))


ax1.legend()

ax2.bar(tcenter,tval,width=1.12,
            color='none',alpha=1,edgecolor='black',linewidth=2,linestyle='dashed',
            label='KM estimator')

ax2.hist(go_tbl['Tau'],bins=7,histtype='step',color='b',#density=1,
        alpha=1,linewidth=2, edgecolor='b', label=r'Constrained $kT_e$')

ax2.legend()
#plt.grid(b=True, which='major', ls=':')
plt.xlabel(r'Optical depth $\tau$',fontsize=15)
ax2.tick_params(axis='both', which='major', labelsize=10, direction='in')
#ax2.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1,decimals=0))


plt.xticks(fontsize=10, rotation=0)
#plt.yticks([5.0,10.0,15.0,20.0,25.0,30.0],[5,10,15,20,25,30],fontsize=15, rotation=0)
plt.tight_layout()


plt.savefig('./figures/tau_dist.pdf',dpi=400)
