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

hist, edges