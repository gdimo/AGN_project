from astropy.table import Table, Column, vstack, MaskedColumn

import matplotlib
import matplotlib.ticker as mtick
import matplotlib.pyplot as plt
import numpy as np

all_tbl = Table.read('./data/complete_table.csv', format='csv')

all_tbl['KT_high'].fill_value = -99.9 
all_tbl['KT_low'].fill_value = -99.9

all_tbl = all_tbl.filled()

msk_kt = all_tbl['KT_high'] < 0
lo_tbl = all_tbl[msk_kt]
go_tbl = all_tbl[~msk_kt]

fig, ax1 = plt.subplots(1,1)
#plt.rcParams['xtick.top'] = True
#plt.rcParams['ytick.right'] = True


logbins_g = np.logspace(2,6,20)
logbins_l = np.logspace(2,5,15)
ax1.hist(go_tbl['rounded_photons'],bins=logbins_g,histtype='step',color='b',
        alpha=1,linewidth=2, edgecolor='b', label=r'Constrained $kT_e$')
ax1.hist(lo_tbl['rounded_photons'],bins=logbins_l,histtype='step',color='k',linestyle='dotted',
        alpha=1,linewidth=2, edgecolor='black', label=r'Low limit $kT_e$')

ax1.legend(loc=2)
ax1.tick_params(axis='both', which='major', labelsize=10, direction='in')
ax1.set_xlabel('Total counts', fontsize=15)
plt.xscale('log')
plt.tight_layout()
#plt.show()
plt.savefig('./figures/photons.pdf',dpi=400)