#Script to bin and plot Theta and Compactness
#import libraries
from astropy.table import Table, Column, vstack, MaskedColumn
from astropy.constants import G, c, m_p, sigma_T, M_sun, L_sun, m_e #import constants
from astropy import units as u
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from scipy import stats
#Load data files

all_tbl = Table.read("./data/all_table.csv", format='csv')

#-----------------------------------------------------------------------------
#Define useful formulas
def EddLumandRatio(mass, L2_to_10):
    mass = (10**(mass))*M_sun #Make mass quantity
    Ledd = (4.*np.pi*c*G*m_p*mass)/sigma_T #Calculate Eddington Luminosity
    
    Ledd = Ledd/10**(44)
    Ledd = Ledd.cgs
    k2_10 = 20
    L_bol = k2_10 * L2_to_10#Calculate the Bolometric Luminosity
    L_bol = L_bol*(u.erg)*(u.s)**-1
    lamda_ratio = (L_bol)/(Ledd)
    return lamda_ratio

def compactness(mass, lamda):
    #Ledd: Eddington Luminosity; Lagn: bolometric luminosity
    mass = mass = 10**(mass)*M_sun
    R_grav = (G*mass)/c**2
    R_source = 10*R_grav
    #Kamraj_2018 eq.1. ; Ricci+2018 (s4)
    kappa_x = 3.87
    l_compactness = 4*np.pi*(m_p/m_e)*(R_grav/R_source)*(lamda/kappa_x)

    return l_compactness   

def brem(theta):
    af = 1./137.
    l_brem = 3*af*theta**(-0.5)
    return l_brem

def el_pr(theta):
    l_el_pr = 0.04*theta**(-3/2)
    return l_el_pr

def el_el(theta):
    l_el_el = 80*theta**(-3/2)
    return l_el_el

def Svenson(theta):
    l_svenson = 10*(theta**(5/2))*(np.exp(theta**(-1)))
    return l_svenson
#-----------------------------------------------------------------------------
#Create tables and calculate quantities
all_tbl['KT_high'].fill_value = -99.9 
all_tbl['KT_low'].fill_value = -99.9

all_tbl = all_tbl.filled()

all_tbl['lamda'] = EddLumandRatio(all_tbl['mass'],all_tbl['lum210'])
all_tbl['lcom'] = compactness(all_tbl['mass'],all_tbl['lamda'])

all_tbl['lcom'] = all_tbl['lcom'].astype(np.float64)
print(np.log10(all_tbl['lcom']))

msk_kt = all_tbl['KT_high'] < 0
lo_tbl = all_tbl[msk_kt]
go_tbl = all_tbl[~msk_kt]

go_tbl['theta_l'] = go_tbl['KT_low']/511
go_tbl['theta_h'] = go_tbl['KT_high']/511

#--------------------------------------------------------------------------
#Calculate and load theoretical lines
#initial conditions

trialtheta = np.logspace(-2,0,500)
# 
# # #Bremsstrahlung
brem = brem(trialtheta)
# Electron-proton
elpr = el_pr(trialtheta)
# Electron-electron
elel = el_el(trialtheta)
# #Svenson
lsven = Svenson(trialtheta)
# # 
# #Slab
slab = Table.read('./data/slabfabian.csv', format='csv')
# #Hemishpere
hem = Table.read('./data/Hem_kamraj.csv', format='csv')
ric = Table.read('./data/ricci.csv', format='csv')
ric['theta'] = np.log10(ric['theta'])
ric['theta_l'] = np.log10(ric['theta_l'])
ric['theta_u'] = np.log10(ric['theta_u'])
ric['lc'] = np.log10(ric['lc'])
ric['lc_l'] = np.log10(ric['lc_l'])
ric['lc_u'] = np.log10(ric['lc_u'])
# 
# tmpg = go_tbl
# tmpl = lo_tbl
# 
# 
# tmpg['ktnew'] = np.zeros(52)
# for j in range(len(tmpg)):
#     tmpg['ktnew'][j] = go_tbl['KT'][j]
# tmpl['ktnew'] = np.zeros(35)
# for i in range(len(tmpl)):
#     tmpl['ktnew'][i]=169.
# 
# 
# 
# tmp = vstack([tmpg,tmpl])
# ali = vstack([go_tbl,lo_tbl])
# 
# 
# ali['th'] = ali['KT']/511.
# tmp['th'] = tmp['ktnew']/511.
# 
# 
# print(ali['th'] - tmp['th'])
# '''
# 
# '''
#Theta - Compactness
#plt.figure(figsize=(10,8))
#plt.rcParams['xtick.top']=True
#plt.rcParams['ytick.right']=True
#
#plt.scatter(go_tbl['theta'],go_tbl['lcom'],c=(go_tbl['lamda']),s=70,
#    alpha=0.8,cmap='seismic',marker='s',label='Constrained')
#
#plt.errorbar(go_tbl['theta'],go_tbl['lcom'],
#    xerr=[(go_tbl['theta']-go_tbl['theta_l']),(go_tbl['theta_h']-go_tbl['theta'])],
#    yerr=False,
#    fmt='k,',
#    capsize=3,
#    ecolor='lightgrey',
#    elinewidth=0.5,
#    )
#
#
#plt.scatter(lo_tbl['theta'],lo_tbl['lcom'],c=(lo_tbl['lamda']),s=70,
#    alpha=0.8,cmap='seismic',marker='>',label='Unconstrained')
#
#plt.errorbar(lo_tbl['theta'],lo_tbl['lcom'],
#    xerr=lo_tbl['theta']/2,
#    xlolims=True,
#    yerr=False,
#    fmt='k,',
#    capsize=3,
#    ecolor='lightgrey',
#    elinewidth=0.5,
#    )
#
#plt.colorbar(label=r'$\lambda_{Edd}$')
#
#plt.plot(hem['theta'],hem['lc'],'g-',linewidth=2.0,label='Hemishpere')
##plt.plot(trialtheta,lc_th,'k-',linewidth=2.0,label='Svensson 1984')
##plt.plot(trialtheta,brem,'y-',linewidth=2.0,label='Bremsstralung')
##plt.text(0.07,3*1e4,'Pair production',rotation=-72,fontsize=15,color='k')
##plt.text(0.012,2.5e-1,'Bremsstrahlung',rotation=-8,fontsize=15,color='b')
#plt.legend(fontsize=15)
#plt.xscale('log')
#plt.yscale('log')
#
#plt.xlim(0.01,2)
#plt.ylim(0.1,100000)
#plt.xlabel(r'$\Theta$ $(kT_{e}/mc^2)$',fontsize=15)
#plt.ylabel(r'Compactness $- \mathit{l}$',fontsize=15)
#
#plt.tight_layout()
#plt.show()

#Bin data

#LOW LIMIT 169
x1 = all_tbl['theta']
x1 = x1.astype(np.float64)
x1 = np.log10(x1)

print('FOR TMP THETA DONE')

y1 = all_tbl['lcom']
#y1 = y1.astype(np.float64)
y1 = np.log10(y1)

print('FOR TMP COMPACTNESS DONE')

bin_meany1, bin_edgesy1, binnumbery1 = stats.binned_statistic(x1,y1,'mean',bins=5)
bin_stdy1, bin_edgesy1, binnumbery1 = stats.binned_statistic(x1,y1,'std',bins=5)

bin_widthy1 = (bin_edgesy1[1] - bin_edgesy1[0])/2
bin_centersy1 = bin_edgesy1[1:] - bin_widthy1

#Low limit = low limit
x2 = go_tbl['theta']
x2 = x2.astype(np.float64)
x2 = np.log10(x2)

print('FOR Detections DATA THETA DONE')

y2 = go_tbl['lcom']
y2 = y2.astype(np.float64)
y2 = np.log10(y2)

print('FOR Detections DATA COMPACTNESS DONE')

bin_meany2, bin_edgesy2, binnumbery2 = stats.binned_statistic(x2,y2,'mean',bins=5)
bin_stdy2, bin_edgesy2, binnumbery2 = stats.binned_statistic(x2,y2,'std',bins=5)

bin_widthy2 = (bin_edgesy2[1] - bin_edgesy2[0])/2
bin_centersy2 = bin_edgesy2[1:] - bin_widthy2


#NEW UPPER LIMIT TABLE
newl = lo_tbl

for i in range(len(newl)):
    newl['KT'][i] = 169

new = vstack([go_tbl,newl])
new['th'] = new['KT']/511

x3 = new['th']
x3 = x3.astype(np.float64)
x3 = np.log10(x3)

print('FOR Detections DATA THETA DONE')

y3 = new['lcom']
y3 = y3.astype(np.float64)
y3 = np.log10(y3)

print('FOR Detections DATA COMPACTNESS DONE')

bin_meany3, bin_edgesy3, binnumbery3 = stats.binned_statistic(x3,y3,'mean',bins=5)
bin_stdy3, bin_edgesy3, binnumbery3 = stats.binned_statistic(x3,y3,'std',bins=5)

bin_widthy3 = (bin_edgesy3[1] - bin_edgesy3[0])/2
bin_centersy3 = bin_edgesy3[1:] - bin_widthy3

print('BINNING COMPLETED SUCCESSFULLY')

x4 = lo_tbl['theta']
x4 = x4.astype(np.float64)
x4 = np.log10(x4)

print('FOR Detections DATA THETA DONE')

y4 = lo_tbl['lcom']
y4 = y4.astype(np.float64)
y4 = np.log10(y4)

print('FOR Detections DATA COMPACTNESS DONE')

bin_meany4, bin_edgesy4, binnumbery4 = stats.binned_statistic(x4,y4,'mean',bins=5)
bin_stdy4, bin_edgesy4, binnumbery4 = stats.binned_statistic(x4,y4,'std',bins=5)

print(bin_edgesy1)

bin_widthy4 = (bin_edgesy4[1] - bin_edgesy4[0])/2
bin_centersy4 = bin_edgesy4[1:] - bin_widthy4

#PLOT
#---------------------------------------------------
##LOW LIMIT = 169
#plt.figure(figsize=(10,8))
#
#plt.plot(bin_centersy1,bin_meany1,'ko',ms=10,label='Low + Detections')
#plt.errorbar(bin_centersy1,bin_meany1,
#        xerr=bin_widthy1,
#         yerr=bin_stdy1,
#         fmt='k,',
#         capsize=3,
#         ecolor='grey',
#         elinewidth=0.5,
#                
#         )
##-----------------------------------------------
##Low limit = low limit
#
#plt.plot(bin_centersy2,bin_meany2,'ro',ms=10,mfc='r',label='Only Detections')
#plt.errorbar(bin_centersy2,bin_meany2,
#         xerr=bin_widthy2,
#         yerr=bin_stdy2,
#         fmt='r,',
#         capsize=3,
#         ecolor='coral',
#         elinewidth=0.5,
#                 
#         )
#
#
#plt.plot(bin_centersy3,bin_meany3,'bs',ms=10,mfc='b',label='Upper Limits')
#plt.errorbar(bin_centersy3,bin_meany3,
#         xerr=bin_widthy3,
#         yerr=bin_stdy3,
#         fmt='b,',
#         capsize=3,
#         ecolor='cyan',
#         elinewidth=0.5,
#                 
#         )
#
#plt.plot(bin_centersy4,bin_meany4,'gs',ms=10,mfc='g',label='Low Limits')
#plt.errorbar(bin_centersy4,bin_meany4,
#         xerr=bin_widthy4,
#         yerr=bin_stdy4,
#         fmt='g,',
#         capsize=3,
#         ecolor='olive',
#         elinewidth=0.5,
#                 
#         )
#
#
#
#
#plt.plot(x2,y2,marker='o',mec='grey',mfc='none',linestyle='none',label='Detection points')
#plt.plot(x4,y4,marker='^',mec='grey',mfc='none',linestyle='none',label='Low Limit points')
##plt.plot(x3,y3,marker='v',mec='cyan',mfc='none',linestyle='none',label='169 keV')
#
#plt.plot(np.log10(hem['theta']),np.log10(hem['lc']),'g-',linewidth=2.0,label='Hemishpere')
## plt.plot(np.log10(trialtheta),np.log10(lc_th),'k-',linewidth=2.0,label='Svensson 1984')
## #plt.plot(np.log10(trialtheta),np.log10(brem),'y-',linewidth=2.0,label='Bremsstralung')
## 
## plt.ylim(0.5,3)
## plt.xlim(-2,1.5)
## 
##
#plt.legend(fontsize=15) 
#plt.ylabel('Compactness')
#plt.xlabel(r'$\Theta$ $keV/mc^2$')
#plt.tight_layout()
#plt.show()
## '''
#---------------------------------------------------------------------------------------------
# ------------------------------------INVERT BINNING------------------------------------------
#---------------------------------------------------------------------------------------------
#LOW LIMIT 169
w1 = all_tbl['theta']
w1 = w1.astype(np.float64)
w1 = np.log10(w1)

print('FOR TMP THETA DONE')

u1 = all_tbl['lcom']
#y1 = y1.astype(np.float64)
u1 = np.log10(u1)

print('FOR TMP COMPACTNESS DONE')

bin_meanu1, bin_edgesu1, binnumberu1 = stats.binned_statistic(u1,w1,'mean',bins=5)
bin_stdu1, bin_edgesu1, binnumberu1 = stats.binned_statistic(u1,w1,'std',bins=5)

bin_widthu1 = (bin_edgesu1[1] - bin_edgesu1[0])/2
bin_centersu1 = bin_edgesu1[1:] - bin_widthu1

      

#Low limit = low limit
w2 = go_tbl['theta']
w2 = w2.astype(np.float64)
w2 = np.log10(w2)

print('FOR Detections DATA THETA DONE')

u2 = go_tbl['lcom']
u2 = u2.astype(np.float64)
u2 = np.log10(u2)

print('FOR Detections DATA COMPACTNESS DONE')

bin_meanu2, bin_edgesu2, binnumberu2 = stats.binned_statistic(u2,w2,'mean',bins=5)
bin_stdu2, bin_edgesu2, binnumberu2 = stats.binned_statistic(u2,w2,'std',bins=5)

bin_widthu2 = (bin_edgesu2[1] - bin_edgesu2[0])/2
bin_centersu2 = bin_edgesu2[1:] - bin_widthu2


#NEW UPPER LIMIT TABLE
newl = lo_tbl

for i in range(len(newl)):
    newl['KT'][i] = 169

new = vstack([go_tbl,newl])
new['th'] = new['KT']/511

w3 = new['th']
w3 = w3.astype(np.float64)
w3 = np.log10(w3)

print('FOR Detections DATA THETA DONE')

u3 = new['lcom']
u3 = u3.astype(np.float64)
u3 = np.log10(u3)

print('FOR Detections DATA COMPACTNESS DONE')

bin_meanu3, bin_edgesu3, binnumberu3 = stats.binned_statistic(u3,w3,'mean',bins=5)
bin_stdu3, bin_edgesu3, binnumberu3 = stats.binned_statistic(u3,w3,'std',bins=5)


bin_widthu3 = (bin_edgesu3[1] - bin_edgesu3[0])/2
bin_centersu3 = bin_edgesu3[1:] - bin_widthu3

print('BINNING COMPLETED SUCCESSFULLY')

w4 = lo_tbl['theta']
w4 = w4.astype(np.float64)
w4 = np.log10(w4)

print('FOR Detections DATA THETA DONE')

u4 = lo_tbl['lcom']
u4 = u4.astype(np.float64)
u4 = np.log10(u4)

print('FOR Detections DATA COMPACTNESS DONE')

bin_meanu4, bin_edgesu4, binnumberu4 = stats.binned_statistic(u4,w4,'mean',bins=5)
bin_stdu4, bin_edgesu4, binnumberu4 = stats.binned_statistic(u4,w4,'std',bins=5)



bin_widthu4 = (bin_edgesu4[1] - bin_edgesu4[0])/2
bin_centersu4 = bin_edgesu4[1:] - bin_widthu4


#
#plt.figure(figsize=(10,8))
##plt.plot(bin_meanu1,bin_centersu1,'mo',ms=10,label='Low + Detections')
##plt.errorbar(bin_meanu1,bin_centersu1,
##         xerr=bin_stdu1,
##         yerr=bin_widthu1,
##         fmt='m,',
##         capsize=3,
##         ecolor='grey',
##         elinewidth=0.5,
##                
##         )
#
#plt.plot(bin_meanu2,bin_centersu2,'ro-',ms=10,mfc='r',label='Only Detections')
#plt.errorbar(bin_meanu2,bin_centersu2,
#         xerr=bin_stdu2,
#         yerr=bin_widthu2,
#         fmt='r,',
#         capsize=3,
#         ecolor='coral',
#         elinewidth=0.5,
#                 
#         )
#
##plt.plot(bin_meanu3,bin_centersu3,'bs',ms=10,mfc='b',label='Upper Limits')
##plt.errorbar(bin_meanu3,bin_centersu3,
##         xerr=bin_stdu3,
##         yerr=bin_widthu3,
##         fmt='b,',
##         capsize=3,
##         ecolor='b',
##         elinewidth=0.5,
##                 
##         )
#
##plt.plot(bin_meanu4,bin_centersu4,'gs',ms=10,mfc='g',label='Low Limits')
##plt.errorbar(bin_meanu4,bin_centersu4,
##         xerr=bin_stdu4,
##         yerr=bin_widthu4,
##         fmt='g,',
##         capsize=3,
##         ecolor='olive',
##         elinewidth=0.5,
##                 
##         )
#
#plt.plot(ric['theta'],ric['lc'],'kD:',ms=10,mfc='k',label='Ricci+18')
#plt.errorbar(ric['theta'],ric['lc'],
#         xerr=[(ric['theta']-ric['theta_l']),(ric['theta_u']-ric['theta'])],
#         yerr=[(ric['lc']-ric['lc_l']),(ric['lc_u']-ric['lc'])],
##         fmt='k,',
##         capsize=3,
##         ecolor='k',
##         elinewidth=0.5,
##         )
##
###plt.plot(x2,y2,marker='o',mec='grey',mfc='none',linestyle='none',label='Detection points')
###plt.plot(x4,y4,marker='^',mec='grey',mfc='none',linestyle='none',label='Low Limit points')
##
##plt.plot(np.log10(hem['theta']),np.log10(hem['lc']),'g-',linewidth=2.0,label='Hemishpere')
##
##
##
##plt.legend(fontsize=15) 
##plt.ylabel('Compactness')
##plt.xlabel(r'$\Theta$ $keV/mc^2$')
##plt.tight_layout()
##
##plt.show()   
##
#
##---------------------------------------------------------------
##----------DUAL AXES plot---------------------------------------
##---------------------------------------------------------------
##fig, [ax1, ax2] = plt.subplots(1,2,figsize=(9,3))
#
#plt.figure(figsize=(10,8))
#plt.rcParams['xtick.top'] = True
#plt.rcParams['ytick.right'] = True
#
#plt.plot(go_tbl['theta'],go_tbl['lcom'],'ks',ms=6,mec='k',mfc='k',label=r'Constrained $kT_e$')
#plt.errorbar(go_tbl['theta'],go_tbl['lcom'],
#                xerr=[(go_tbl['theta']-go_tbl['theta_l']),(go_tbl['theta_h']-go_tbl['theta'])],
#                yerr=False,
#                fmt='none',
#                capsize=3,
#                ecolor='k',
#                elinewidth=0.5
#                )
#
#plt.plot(lo_tbl['theta'],lo_tbl['lcom'],'ko',ms=7,mec='k',mfc='none',label=r'Unconstrained $kT_e$')
#plt.errorbar(lo_tbl['theta'],lo_tbl['lcom'],
#                xerr=lo_tbl['theta']/2,
#                xlolims=True,
#                yerr=False,
#                fmt='none',
#                capsize=3,
#                ecolor='k',
#                elinewidth=0.5
#                )
#
##ax1.plot(trialtheta,brem,'b--')
##ax1.plot(trialtheta,elpr,'r--',label=r'$e^--p^+$')
##ax1.plot(trialtheta,elel,'y--',label=r'$e^--e^-$')
##plt.plot(trialtheta,lsven,'b:',label='Svensson')
#plt.plot(hem['theta'],hem['lc'],'g--',linewidth=2.0,label='Hemisphere')
#plt.plot(slab['theta'],slab['lc'],'r--',linewidth=2.0,label='Slab')
#
#plt.xticks(fontsize=15, rotation=0)
#plt.yticks(fontsize=15, rotation=0)
#plt.ylabel('Compactness (l)',fontsize=14)
#plt.xlabel(r'$\Theta$ ($kT_e/mc^2$)',fontsize=14)
#plt.xscale('log')
#plt.yscale('log')
#plt.xlim(0.01,1.1)
#plt.ylim(0.1,1000)
#
#plt.legend(loc=4, fontsize=14)
#plt.tight_layout()
#plt.savefig('../figures/paper/compac_scatter.pdf', dpi=400)
##-------------------------------------------------------------------
##       BIN
##--------------------------------------------------------------------------
#plt.figure(figsize=(10,8))
#plt.rcParams['xtick.top'] = True
#plt.rcParams['ytick.right'] = True
#plt.plot(bin_meanu1,bin_centersu1,'ks',ms=10,label='All sources')
#plt.errorbar(bin_meanu1,bin_centersu1,
#         xerr=bin_stdu1,
#         yerr=bin_widthu1,
#         fmt='k,',
#         capsize=3,
#         ecolor='k',
#         elinewidth=0.5
#         )
#
#plt.plot(bin_meanu2,bin_centersu2,'ko',ms=10,mfc='none',label='Only Detections')
#plt.errorbar(bin_meanu2,bin_centersu2,
#         xerr=bin_stdu2,
#         yerr=bin_widthu2,
#         fmt='k,',
#         capsize=3,
#         ecolor='k',
#         elinewidth=0.5
#         )
#
#plt.plot(ric['theta'],ric['lc'],'kD',ms=10,mfc='k',label='Ricci+18')
#plt.errorbar(ric['theta'],ric['lc'],
#         xerr=[(ric['theta']-ric['theta_l']),(ric['theta_u']-ric['theta'])],
#         yerr=[(ric['lc']-ric['lc_l']),(ric['lc_u']-ric['lc'])],
#         fmt='k,',
#         capsize=3,
#         ecolor='k',
#         elinewidth=0.5,
#         )
#
#plt.plot(np.log10(hem['theta']),np.log10(hem['lc']),'g--',linewidth=2.0,label='Hemisphere')
#plt.plot(np.log10(slab['theta']),np.log10(slab['lc']),'r--',linewidth=2.0,label='Slab')
#
#plt.xlim(-1.75,0.25)
#plt.ylim(-0.5,3.5)
#plt.xticks(fontsize=15, rotation=0)
#plt.yticks(fontsize=15, rotation=0)
#plt.ylabel('Compactness (l)',fontsize=14)
#plt.xlabel(r'$\Theta$ ($kT_e/mc^2$)',fontsize=14)
#plt.legend(loc=1, fontsize=14)
#
#plt.tight_layout()
#plt.savefig('../figures/paper/compac_BIN.pdf', dpi=400)
#
#
##----------------------------------------------------------------------
##----ONE FIGURE TWO PANELS---------------------------------------------
##----------------------------------------------------------------------
#
fig, [ax1, ax2] = plt.subplots(1,2,figsize=(14,6))
plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True

ax1.plot(go_tbl['theta'],go_tbl['lcom'],'ks',ms=6,mec='k',mfc='k',label=r'Constrained $kT_e$')
ax1.errorbar(go_tbl['theta'],go_tbl['lcom'],
                xerr=[(go_tbl['theta']-go_tbl['theta_l']),(go_tbl['theta_h']-go_tbl['theta'])],
                yerr=False,
                fmt='none',
                capsize=3,
                ecolor='k',
                elinewidth=0.5
                )

ax1.plot(lo_tbl['theta'],lo_tbl['lcom'],'ko',ms=7,mec='k',mfc='none',label=r'Unconstrained $kT_e$')
ax1.errorbar(lo_tbl['theta'],lo_tbl['lcom'],
                xerr=lo_tbl['theta']/2,
                xlolims=True,
                yerr=False,
                fmt='none',
                capsize=3,
                ecolor='k',
                elinewidth=0.5
                )

#ax1.plot(trialtheta,brem,'b--')
#ax1.plot(trialtheta,elpr,'r--',label=r'$e^--p^+$')
#ax1.plot(trialtheta,elel,'y--',label=r'$e^--e^-$')
#plt.plot(trialtheta,lsven,'b:',label='Svensson')
ax1.plot(hem['theta'],hem['lc'],'g--',linewidth=2.0,label='Hemisphere')
ax1.plot(slab['theta'],slab['lc'],'r--',linewidth=2.0,label='Slab')


ax1.set_ylabel('Compactness (l)',fontsize=14)
ax1.set_xlabel(r'$\Theta$ ($kT_e/mc^2$)',fontsize=14)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(0.01,1.1)
ax1.set_ylim(0.1,1000)

ax1.legend(loc=4, fontsize=14)

#-------------------------------------------------------------------
#       BIN
#--------------------------------------------------------------------------

plt.rcParams['xtick.top'] = True
plt.rcParams['ytick.right'] = True

#------------WE SHOULD CREATE KM POINTS FOR THESE POINTS
#We will create a table of the values that will be used and then assign them to the varieables plotted
#   _______________________________________________
#   | x(theta)    |   y(compact)            |   range          |
#   |   -1.032 +/0.116   |      0.83         |     0.31 -- 1.22        |
#   |   -0.963 +/- 0.008   |     1.39          |    1.22 -- 1.48         |
#   |   -0.740 +/- 0.061      |   1.66    |  1.50 -- 1.80  |
#   |   -1.014 +/- 0.067    |    1.93    |   1.81 -- 2.1   |
#   |    -1.037 +/- 0.052   |  2.42   |   2.1 -- 3.38   |
#  
# Build the tables
bin_meanu1 = np.array([-1.032, -0.963, -0.740, -1.014, -1.037])
bin_centersu1 = np.array([0.83, 1.39, 1.66, 1.93, 2.42])
#bin_widthu1 = np.array([0.24, 0.3, 0.3, 0.3, (0.3,0.5)])
ss_low = [0.12, 0.15, 0.15, 0.15, 0.3]
ss_up = [0.12, 0.15, 0.15, 0.15, 1.0]
bin_stdu1 = np.array([0.116, 0.08, 0.061, 0.067, 0.052])   
ax2.plot(bin_meanu1,bin_centersu1,'ks',ms=10,label='KM sources')
ax2.errorbar(bin_meanu1,bin_centersu1,
         xerr=bin_stdu1,
         yerr=[ss_low, ss_up],
         fmt='k,',
         capsize=3,
         ecolor='k',
         elinewidth=0.5
         )
#------ END OF KM POINTS-------------------------
ax2.plot(bin_meanu2,bin_centersu2,'ko',ms=10,mfc='none',label='Only Detections')
ax2.errorbar(bin_meanu2,bin_centersu2,
         xerr=bin_stdu2,
         yerr=bin_widthu2,
         fmt='k,',
         capsize=3,
         ecolor='k',
         elinewidth=0.5
         )

ax2.plot(ric['theta'],ric['lc'],'kD',ms=10,mfc='k',label='Ricci+18')
ax2.errorbar(ric['theta'],ric['lc'],
         xerr=[(ric['theta']-ric['theta_l']),(ric['theta_u']-ric['theta'])],
         yerr=[(ric['lc']-ric['lc_l']),(ric['lc_u']-ric['lc'])],
         fmt='k,',
         capsize=3,
         ecolor='k',
         elinewidth=0.5,
         )

ax2.plot(np.log10(hem['theta']),np.log10(hem['lc']),'g--',linewidth=2.0,label='Hemisphere')
ax2.plot(np.log10(slab['theta']),np.log10(slab['lc']),'r--',linewidth=2.0,label='Slab')

ax2.set_xlim(-1.75,0.25)
ax2.set_ylim(-0.5,3.5)

ax2.set_ylabel('Compactness (l)',fontsize=14)
ax2.set_xlabel(r'$\Theta$ ($kT_e/mc^2$)',fontsize=14)
ax2.legend(loc=1, fontsize=14)

fig.tight_layout()
fig.savefig('./figures/compac_all.pdf', dpi=400)