l#python 3Dplot.py
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rc
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)


##############################Reading in & assign data##############################
#positionVec = np.loadtxt('/mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/output.00399.Xi')
##xi = positionVec[:,1]
##yi = positionVec[j,2]
##zi = positionVec[j,3]

positionStar = np.loadtxt('/mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/output.00399.stellar_rp.Xi_star')


positionCenter = np.loadtxt('/mn/stornext/d17/extragalactic/personal/gabrierg/Eris_AHF/L90Mpc8000_hithres.00400.amiga.stt')
xc = positionCenter[0,13]
yc = positionCenter[0,14]
zc = positionCenter[0,15]

#Eumass = np.loadtxt('/mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/output.00399.eu_out')        #1st run
##Eumass = np.loadtxt('/mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/output.00399.ord_sml.3dplot') #2nd run
##Eumass = np.loadtxt('/mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/output.00400.ord_sml.3dplot') #2nd run
##euid = Eumass[:,0]#12936000
#eumass = Eumass[:,1]
#eu = Eumass[:,2]
##eugrp = Eumass[:,3]

#Femass = np.loadtxt('/mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/iron.fe_smooth_128')
#fe = Femass[:,:]

#Starmass = np.loadtxt('/mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/stellar.seu_sm_0p01_128')
#solar = Starmass[:,:]

#Dwarfmass = np.loadtxt('/mn/stornext/d17/extragalactic/personal/gabrierg/rprocess_halos/dwarf.seu_dwarf_0p01_128')
#dw = Dwarfmass[:,:]

#print(np.shape(Eumass2))
#print(np.shape(Eumass))
#print(np.shape(Starmass))
#print(np.shape(Dwarfmass))
#print(np.shape(positionCenter))


####################################Shorten arrays###################################
#EU
#size = len(eu)
#print(size) #12936000

#l=0
#Eu_massarry = np.zeros((9611, 2))
#for i in range(0,8649900,900):
#    Eu_massarry[l,0] = euid[i]
#    Eu_massarry[l,1] = eumass[i]
#    l=l+1


#STAR POSITION
n=0
pos_istar = np.zeros((9611, 3))
for j in range(0,7491000,1000):
    pos_istar[n,0] = positionStar[j,1]
    pos_istar[n,1] = positionStar[j,2]
    pos_istar[n,2] = positionStar[j,3]
    n=n+1


#GAS POSITION
#m=0
#pos_i = np.zeros((9611, 3))
#for k in range(0,8649900,900):
##    if Eu_massarry[m,0] == positionVec[j,0]:
#    pos_i[m,0] = positionVec[k,1]
#    pos_i[m,1] = positionVec[k,2]
#    pos_i[m,2] = positionVec[k,3]
#    m=m+1

#for r in range(20):
#    print(Euarry[r,:])
#for w in range(20):
#    print(pos_i[w,:])


#IRON
#p=0
#Fe=np.zeros(8649900)
#for i in range(len(fe)):
#    for j in range(6):
#        Fe[p]=fe[i,j]
#        p=p+1

#print(len(fe))
#print(np.shape(fe))
#print(len(Fe))
#q=0
#iron = np.zeros(9611)
#for j in range(0,8649900,900):
#    iron[q] = Fe[j]
#    q=q+1


#STAR
#pp=0
#StellR=np.zeros(8649900)
#for i in range(len(solar)):
#    for j in range(6):
#        StellR[pp]=solar[i,j]
#        pp=pp+1

#qq=0
#star = np.zeros(9611)
#for j in range(0,8649900,900):
#    star[qq] = StellR[j]
#    qq=qq+1


#DWARF
#ppp=0
#Df=np.zeros(8649900)
#for i in range(len(dw)):
#    for j in range(6):
#        Df[ppp]=dw[i,j]
#        ppp=ppp+1

#qqq=0
#dwarf = np.zeros(9611)
#for j in range(0,8649900,900):
#    dwarf[qqq] = Df[j]
#    qqq=qqq+1


####################################Change units####################################
kpcunit = 1000
Mpc = 1000000
lunit = 90#*Mpc
Xc_codeunit=(xc/lunit)-0.5
Yc_codeunit=(yc/lunit)-0.5
Zc_codeunit=(zc/lunit)-0.5
#print(Xc_codeunit)
#print(Yc_codeunit)
#print(Zc_codeunit)

#DeltaPos = np.zeros((9611, 3))
#for k in range(9611):
#    DeltaPos[k,0] = (pos_i[k,0] - Xc_codeunit) * lunit * kpcunit
#    DeltaPos[k,1] = (pos_i[k,1] - Yc_codeunit) * lunit * kpcunit
#    DeltaPos[k,2] = (pos_i[k,2] - Zc_codeunit) * lunit * kpcunit

DeltaPos_star = np.zeros((9611, 3))
for k in range(9611):
    DeltaPos_star[k,0] = (pos_istar[k,0] - Xc_codeunit) * lunit * kpcunit
    DeltaPos_star[k,1] = (pos_istar[k,1] - Yc_codeunit) * lunit * kpcunit
    DeltaPos_star[k,2] = (pos_istar[k,2] - Zc_codeunit) * lunit * kpcunit


#Eu_Fe = np.zeros(9611)
#for t in range(9611):
#    if iron[t] == 0:
#        Eu_Fe[t] = 0
#    else:
#        Eu_Fe[t] = Eu_massarry[t,1]/iron[t]

#dwarfEu_Fe = np.zeros(9611)
#for v in range(9611):
#    if iron[v] == 0:
#        dwarfEu_Fe[v] = 0
#    else:
#        dwarfEu_Fe[v] = dwarf[v]/iron[v]


######################################Plot data######################################
#np.savetxt("Position.txt", X)
#X = np.zeros(9611)
#Y = np.zeros(9611)
#Z = np.zeros(9611)

#b=0
#for n in range(9611):
#    X[n] = DeltaPos[n,0]
#    Y[n] = DeltaPos[n,1]
#    Z[n] = DeltaPos[n,2]


Xstar = np.zeros(9611)
Ystar = np.zeros(9611)
Zstar = np.zeros(9611)

b=0
for n in range(9611):
    Xstar[n] = DeltaPos_star[n,0]
    Ystar[n] = DeltaPos_star[n,1]
    Zstar[n] = DeltaPos_star[n,2]


#    if DeltaPos[n,0] >= 10000000 or DeltaPos[n,0] <= -10000000:
#        X[n] = b
#    if DeltaPos[n,1] >= 10000000 or DeltaPos[n,1] <= -10000000:
#        Y[n] = b
#    if DeltaPos[n,2] >= 10000000 or DeltaPos[n,2] <= -10000000:
#        Z[n] = b


########################################Plot#########################################
print("go for plot")

####################plot of [Eu/Fe] vs [Fe/H]
#fig=plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.scatter(iron, Eu_Fe, s=5, c='r', marker="o", label='second')
#ax1.set_ylabel('Eu/Fe')                                                   
#ax1.set_xlabel('Fe/H')                                                   
#ax1.set_title('Eu vs Fe abondance')                            
#plt.xlim(-2500,2500)
#plt.ylim(-1500,3000)
#plt.show()

####################plot of dwarf [Eu/Fe] vs [Fe/H]
#fig=plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.scatter(iron, dwarfEu_Fe, s=5, c='r', marker="o", label='second')
#ax1.set_ylabel('Dwarf Eu/Fe')
#ax1.set_xlabel('Fe/H')
#ax1.set_title('Dwarf Eu vs Fe abondance')
#plt.xlim(-2500,2500)
#plt.ylim(-1500,3000)
#plt.show()

#
#fig=plt.figure()
#ax1 = fig.add_subplot(111)        
#ax1.scatter(iron, Euarry[:,1], s=5, c='r', marker="o", label='second')
#ax1.set_ylabel('Eu')                                                    
#ax1.set_xlabel('Fe')                                                    
#ax1.set_title('Eu vs Fe abondance')
#plt.xlim(-2500,2500)
#plt.ylim(-1500,3000)
#plt.show()

#fig=plt.figure()
#ax=fig.add_axes([0,0,1,1])
#ax.scatter(X, Y, color='r')
#ax.scatter(X, Z, color='b')
#ax.set_xlabel('$\Delta X$')
#ax.set_ylabel('$\Delta Y$')
#ax.set_title('scatter plot')
#plt.xlim(-2500,2500)
#plt.ylim(-1500,3000)
#plt.show()


########################################GAS##########################################
####################plot XvsY####################
#fig=plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.scatter(X, Y, s=3, c='r', marker="o", label='second')
#ax1.set_xlabel('$\Delta X$')#kpc
#ax1.set_ylabel('$\Delta Y$')#kpc
#ax1.set_title('Eu particles from center of host galaxy')
#plt.xlim(-2500,2500)
#plt.ylim(-1500,3000)
#plt.show()

#fig=plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.scatter(X, Y, s=3, c='r', marker="o", label='second')
#ax1.set_xlabel('$\Delta X$')#kpc
#ax1.set_ylabel('$\Delta Y$')#kpc
#ax1.set_title('Stars from center of host galaxy')
#plt.xlim(-2000,2000)
#plt.ylim(-1500,2000)
#plt.show()

#fig=plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.scatter(X, Y, s=3, c='r', marker="o", label='second')
#ax1.set_xlabel('$\Delta X$')#kpc
#ax1.set_ylabel('$\Delta Y$')#kpc
#ax1.set_title('Stars from center of host galaxy')
#plt.xlim(-1000,1000)
#plt.ylim(-1000,1500)
#plt.show()

#fig=plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.scatter(X, Y, s=3, c='r', marker="o", label='second')
#ax1.set_xlabel('$\Delta X$')#kpc
#ax1.set_ylabel('$\Delta Y$')#kpc
#ax1.set_title('Stars from center of host galaxy')
#plt.xlim(-500,500)
#plt.ylim(-1000,1000)
#plt.show()

#fig=plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.scatter(X, Y, s=3, c='r', marker="o", label='second')
#ax1.set_xlabel('$\Delta X$')#kpc
#ax1.set_ylabel('$\Delta Y$')#kpc
#ax1.set_title('Stars from center of host galaxy')
#plt.xlim(-200,200)
#plt.ylim(-500,500)
#plt.show()

#fig=plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.scatter(X, Y, s=3, c='r', marker="o", label='second')
#ax1.set_xlabel('$\Delta X$')#kpc
#ax1.set_ylabel('$\Delta Y$')#kpc
#ax1.set_title('Stars from center of host galaxy')
#plt.xlim(-100,100)
#plt.ylim(-500,500)
#plt.show()



#fig=plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.scatter(X, Y, s=10, c='r', marker="o", label='second')
#ax1.set_xlabel('$\Delta X$')#kpc
#ax1.set_ylabel('$\Delta Y$')#kpc
#ax1.set_title('Stars from center of host galaxy')
##plt.xlim(-2500,2500)
##plt.ylim(1.1e292,1.2e292)
###plt.xlim(-1.3e210,-1.4e210)
###plt.ylim(-1500,3000)
#plt.show()

####################plot XvsZ####################
#fig=plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.scatter(X, Z, s=10, c='r', marker="o", label='second')
#ax1.set_xlabel('$\Delta X$')#kpc
#ax1.set_ylabel('$\Delta Z$')#kpc
#ax1.set_title('Eu particles from center of host galaxy')
#plt.xlim(-2500,2500)
#plt.ylim(-1500,3000)
#plt.show()

#fig=plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.scatter(X, Y, s=10, c='r', marker="o", label='second')
#ax1.set_xlabel('$\Delta X$')#kpc
#ax1.set_ylabel('$\Delta Z$')#kpc
#ax1.set_title('Stars from center of host galaxy')
#plt.xlim(-2000,2000)
#plt.ylim(-1500,2000)
#plt.show()

#fig=plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.scatter(X, Y, s=10, c='r', marker="o", label='second')
#ax1.set_xlabel('$\Delta X$')#kpc
#ax1.set_ylabel('$\Delta Z$')#kpc
#ax1.set_title('Stars from center of host galaxy')
#plt.xlim(-1000,1000)
#plt.ylim(-1000,1500)
#plt.show()

#fig=plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.scatter(X, Y, s=10, c='r', marker="o", label='second')
#ax1.set_xlabel('$\Delta X$')#kpc
#ax1.set_ylabel('$\Delta Z$')#kpc
#ax1.set_title('Stars from center of host galaxy')
#plt.xlim(-500,500)
#plt.ylim(-1000,1000)
#plt.show()

#fig=plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.scatter(X, Y, s=10, c='r', marker="o", label='second')
#ax1.set_xlabel('$\Delta X$')#kpc
#ax1.set_ylabel('$\Delta Z$')#kpc
#ax1.set_title('Stars from center of host galaxy')
#plt.xlim(-200,200)
#plt.ylim(-500,500)
#plt.show()

#fig=plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.scatter(X, Y, s=10, c='r', marker="o", label='second')
#ax1.set_xlabel('$\Delta X$')#kpc
#ax1.set_ylabel('$\Delta Z$')#kpc
#ax1.set_title('Stars from center of host galaxy')
#plt.xlim(-100,100)
#plt.ylim(-500,500)
#plt.show()



#fig=plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.scatter(X, Z, s=10, c='r', marker="o", label='second')
#ax1.set_xlabel('$\Delta X$')#kpc
#ax1.set_ylabel('$\Delta Z$')#kpc
#ax1.set_title('Stars from center of host galaxy')
##plt.xlim(-2500,2500)
##plt.ylim(1e292,2e292)
###plt.xlim(-1.3e210,-1.4e210)
###plt.ylim(-1500,3000)
#plt.show()

########################################STAR#########################################
####################plot XvsY####################
fig=plt.figure()
ax1 = fig.add_subplot(111)
ax1.scatter(Xstar, Ystar, s=5, c='r', marker="o", label='second')
ax1.set_xlabel('$\Delta X$')#kpc
ax1.set_ylabel('$\Delta Y$')#kpc
ax1.set_title('Eu particles from center of host galaxy')
plt.xlim(-2500,2500)
plt.ylim(-1500,3000)
#plt.show()

#fig=plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.scatter(Xstar, Ystar, s=2, c='r', marker="o", label='second')
#ax1.set_xlabel('$\Delta X$')#kpc
#ax1.set_ylabel('$\Delta Y$')#kpc
#ax1.set_title('Stars from center of host galaxy')
#plt.xlim(-2000,2000)
#plt.ylim(-1500,2000)
#plt.show()

#fig=plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.scatter(Xstar, Ystar, s=2, c='r', marker="o", label='second')
#ax1.set_xlabel('$\Delta X$')#kpc
#ax1.set_ylabel('$\Delta Y$')#kpc
#ax1.set_title('Stars from center of host galaxy')
#plt.xlim(-1000,1000)
#plt.ylim(-1000,1500)
#plt.show()

#fig=plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.scatter(Xstar, Ystar, s=2, c='r', marker="o", label='second')
#ax1.set_xlabel('$\Delta X$')#kpc
#ax1.set_ylabel('$\Delta Y$')#kpc
#ax1.set_title('Stars from center of host galaxy')
#plt.xlim(-500,500)
#plt.ylim(-1000,1000)
#plt.show()

#fig=plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.scatter(Xstar, Ystar, s=2, c='r', marker="o", label='second')
#ax1.set_xlabel('$\Delta X$')#kpc
#ax1.set_ylabel('$\Delta Y$')#kpc
#ax1.set_title('Stars from center of host galaxy')
#plt.xlim(-200,200)
#plt.ylim(-500,500)
#plt.show()

#fig=plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.scatter(Xstar, Ystar, s=1, c='r', marker="o", label='second')
#ax1.set_xlabel('$\Delta X$')#kpc
#ax1.set_ylabel('$\Delta Y$')#kpc
#ax1.set_title('Stars from center of host galaxy')
#plt.xlim(-100,100)
#plt.ylim(-500,500)
#plt.show()

####################plot XvsZ####################
#fig=plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.scatter(X, Z, s=10, c='r', marker="o", label='second')
#ax1.set_xlabel('$\Delta X$')#kpc
#ax1.set_ylabel('$\Delta Z$')#kpc
#ax1.set_title('Eu particles from center of host galaxy')
#plt.xlim(-2500,2500)
#plt.ylim(-1500,3000)
#plt.show()

#fig=plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.scatter(X, Y, s=10, c='r', marker="o", label='second')
#ax1.set_xlabel('$\Delta X$')#kpc
#ax1.set_ylabel('$\Delta Z$')#kpc
#ax1.set_title('Stars from center of host galaxy')
#plt.xlim(-2000,2000)
#plt.ylim(-1500,2000)
#plt.show()

#fig=plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.scatter(X, Y, s=10, c='r', marker="o", label='second')
#ax1.set_xlabel('$\Delta X$')#kpc
#ax1.set_ylabel('$\Delta Z$')#kpc
#ax1.set_title('Stars from center of host galaxy')
#plt.xlim(-1000,1000)
#plt.ylim(-1000,1500)
#plt.show()

#fig=plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.scatter(X, Y, s=10, c='r', marker="o", label='second')
#ax1.set_xlabel('$\Delta X$')#kpc
#ax1.set_ylabel('$\Delta Z$')#kpc
#ax1.set_title('Stars from center of host galaxy')
#plt.xlim(-500,500)
#plt.ylim(-1000,1000)
#plt.show()

#fig=plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.scatter(X, Y, s=10, c='r', marker="o", label='second')
#ax1.set_xlabel('$\Delta X$')#kpc
#ax1.set_ylabel('$\Delta Z$')#kpc
#ax1.set_title('Stars from center of host galaxy')
#plt.xlim(-200,200)
#plt.ylim(-500,500)
#plt.show()

#fig=plt.figure()
#ax1 = fig.add_subplot(111)
#ax1.scatter(X, Y, s=10, c='r', marker="o", label='second')
#ax1.set_xlabel('$\Delta X$')#kpc
#ax1.set_ylabel('$\Delta Z$')#kpc
#ax1.set_title('Stars from center of host galaxy')
#plt.xlim(-100,100)
#plt.ylim(-500,500)
#plt.show()
