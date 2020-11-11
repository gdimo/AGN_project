#### Create Histograms using Numpy library
#import the necessary libraries and load data files

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

#create KT histogram

hist, edges = np.histogram(go_tbl['KT'],bins=7, normed=True)
center = np.array([edges+(12)])
height = hist*21.5
print(type(edges))
print(edges)
print(center)
print(height)
center_2 = np.array([20.57, 42.19, 63.80, 85.42, 128.66, 150.28, 171.9])


plt.bar(center_2,height,width=21,
            color='none',alpha=1,edgecolor='black',linewidth=2,linestyle='dashed',
            label='good')
#
plt.show()