#trying stellar to halo mass ratio (z=0) for MR + MRII

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

# First need to define hubble
# Use same value as script_lf # 'where the Hubble parameter is 100h km/s/Mpc - h is called "hubble" in the script' #
h=0.673

#MR & MRII datasets old model (Hen14) for z=0.00,1.04,2.07,3.11
boxside_MR=480.8 #In Mpc/h
Hen0_MR='data/Hen14_sfh2/snap_58.pkl'
Hen1_MR='data/Hen14_sfh2/snap_38.pkl'
Hen2_MR='data/Hen14_sfh2/snap_30.pkl'
Hen3_MR='data/Hen14_sfh2/snap_25.pkl'

boxside_MRII=96.0558 
Hen0_MRII='data/Hen14_MRII/snap_62.pkl'
Hen1_MRII='data/Hen14_MRII/snap_42.pkl'
Hen2_MRII='data/Hen14_MRII/snap_34.pkl'
Hen3_MRII='data/Hen14_MRII/snap_29.pkl'

#MR & MRII data for new model -- Dev17
Dev0_MR='data/Dev17/snap_58.pkl'
Dev1_MR='data/Dev17/snap_38.pkl'
Dev2_MR='data/Dev17/snap_30.pkl'
Dev3_MR='data/Dev17/snap_25.pkl'

Dev0_MRII='data/Dev17_MRII/snap_62.pkl'
Dev1_MRII='data/Dev17_MRII/snap_42.pkl'
Dev2_MRII='data/Dev17_MRII/snap_34.pkl'
Dev3_MRII='data/Dev17_MRII/snap_29.pkl'

#redshifts=[0,1,2,3]
redshifts=1,

# Variables used for the graph
Xlabel=r'$M_h$/$h^{-1}10^{10}M_\odot$'
Ylabel=r'$\log_{10}($M$_{*}$/$M_h$)'
ylimits=[1e-4,1]

# Define mass bins in log space
bin_edges=np.logspace(-2,5,50)
nbins=len(bin_edges)-1
    
#Loop over all redshifts
for redshift in redshifts:

#Loading in galaxy data -- Load in Stellar Mass for MR & MRII
    pickleFile=eval('Hen'+str(redshift)+'_MR')
    exec(open('script_unpickle.py').read())
    index=np.where((gals['Type']==0) & (gals['Mvir']>20))
    gals0=gals[index]
    pickleFile=eval('Hen'+str(redshift)+'_MRII')
    exec(open('script_unpickle.py').read())
    index=np.where((gals['Type']==0) & (gals['Mvir']<=20))
    gals0=np.append(gals0,gals[index])
    mass_stellar=gals0['StellarMass']
    mass_virial=gals0['Mvir']

#Defining axis of plot
#    y1_stellar=mass_stellar
#    y1_virial=mass_virial

#    x1=(y1_virial)
#    y1=np.log10((y1_stellar/y1_virial))

    plt.close()
    plt.figure(redshift)
    plt.fig1=plt.semilogx(x1,y1,'r,',zorder=-1)
    plt.grid(True)

# Divide galaxies into mass bins using digitize function & average properties
    bins=np.digitize(y1_virial,bin_edges)

#Define percentile arrays -- create empty array for 
    low=np.empty(nbins)
    median=np.empty(nbins)
    high=np.empty(nbins)

# Loop over mass bins calculating properties
    for ibin in np.arange(1,len(bin_edges)):
        index=np.where(bins==ibin)[0]
        mvir_mean=np.mean(y1_virial[index]) #mean of viral mass
        mratio_mean=np.mean(y1[index]) #mean of stellar mass/viral mass
        std1=np.std(y1[index]) #standard deviation of y axis (ratio of stellar/viral)
        plt.plot(mvir_mean,mratio_mean,marker='o',color='c',linestyle='-')
#        plt.errorbar(0.5*(bin_edges[ibin-1]+bin_edges[ibin]),mratio_mean,yerr=std1,fmt='o-',capthick=2,color='c')#,label='Hen14 - Mean & Std')
#        if len(index) > 0:
#            low[ibin-1],median[ibin-1],high[ibin-1]=np.percentile(y1[index],[16,50,84])
#    plt.errorbar(0.5*(bin_edges[ibin-1]+bin_edges[ibin]),mratio_mean,yerr=std1,fmt='o',capthick=2,color='c',label='Mean+stdDev')

#Median & percentiles
#    plt.plot(0.5*(bin_edges[0:nbins]+bin_edges[1:nbins+1]),low,':c')
#    plt.plot(0.5*(bin_edges[0:nbins]+bin_edges[1:nbins+1]),median,'-c',label='median')
#    plt.plot(0.5*(bin_edges[0:nbins]+bin_edges[1:nbins+1]),high,':c',label='16 & 84 percentiles')

#--------------------------------------------------------------------#
# ADDING OBS DATA
#SK, More 2011 for z=0:
    SKMore2011x=[(h*10**1.930),(h*10**2.035),(h*10**2.151),(h*10**2.314),(h*10**2.500),(h*10**2.895),(h*10**2.651),(h*10**3.174),(h*10**3.384),(h*10**3.698)]
    SKMore2011y=[-1.620,-1.552,-1.484,-1.450,-1.450,-1.433,-1.365,-1.522,-1.531,-1.701]

    if redshift == 0:
        plt.plot(SKMore2011x,SKMore2011y,'k^', label='SK, More (2011) - z=0')

#Foucaud (2010), z=1:
    Foucaudx=[(h*10**2.698),(h*10**2.791)]
    Foucaudy=[-1.204,-1.280]

    if redshift == 1:
        plt.plot(Foucaudx,Foucaudy,'orchid^')

#Behroozi (2013), z=1,2,3,4,5,6,7
    Behroozix1=[h*10**0.873,h*10**1.248,h*10**1.622,h*10**1.996,h*10**2.371,h*10**2.745,h*10**3.119,h*10**3.494,h*10**3.868,h*10**4.242,h*10**4.626]
    Behrooziy1=[-2.588,-2.334,-1.881,-1.579,-1.606,-1.811,-2.098,-2.420,-2.720,-3.049,-3.404]

    if redshift == 1:
        plt.plot(Behroozix1,Behrooziy1,'navyv')

    Behroozix2=[h*10**1.248,h*10**1.622,h*10**1.996,h*10**1.371,h*10**2.745,h*10**3.119,h*10**3.494,h*10**3.868]
    Behrooziy2=[-2.484,-2.093,-1.695,-1.606,-1.777,2.057,-2.378,-2.727]

    if redshift == 2:
        plt.plot(Behroozix2,Behrooziy2,'navyv')

    Behroozix3=[h*10**1.497,h*10**1.872,h*10**2.246,h*10**2.620,h*10**2.995,h*10**3.369]
    Behrooziy3=[-2.169,-1.778,-1.627,-1.763,-2.030,-2.337]

    if redshift == 3:
        plt.plot(Behroozix3,Behrooziy3,'navyv')

    Behroozix4=[h*10**1.373,h*10**1.747,h*10**2.121,h*10**2.496,h*10**2.870,h*10**3.244]
    Behrooziy4=[-2.128,-1.778,-1.654,-1.784,-2.057,-2.365]

    if redshift == 4:
        plt.plot(Behroozix4,Behrooziy4,'navyv')

    Behroozix5=[h*10**1.248,h*10**1.622,h*10**1.988,h*10**2.371,h*10**2.745]
    Behrooziy5=[-2.039,-1.737,-1.682,-1.852,-2.126]

    if redshift == 5:
        plt.plot(Behroozix5,Behrooziy5,'navyv')

    Behroozix6=[h*10**1.123,h*10**1.488,h*10**1.872,h*10**2.246]
    Behrooziy6=[-1.943,-1.703,-1.716,-1.935]

    if redshift == 6:
        plt.plot(Behroozix6,Behrooziy6,'navyv')

    Behroozix7=[h*10**0.374,h*10**0.749,h*10**1.123,h*10**1.497]
    Behrooziy7=[-2.458,-2.088,-1.779,-1.689]

    if redshift == 7:
        plt.plot(Behroozix7,Behrooziy7,'navyv')

#---------------------------------------------------------------------#
#Dev17 loop 
for redshift in redshifts:

#Loading in galaxy data -- Load in Stellar Mass for MR
    pickleFile=eval('Dev'+str(redshift)+'_MR')
    exec(open('script_unpickle.py').read())
    index=np.where((gals['Type']==0) & (gals['Mvir']>20))
    gals0=gals[index]
    pickleFile=eval('Dev'+str(redshift)+'_MRII')
    exec(open('script_unpickle.py').read())
    index=np.where((gals['Type']==0) & (gals['Mvir']<=20))
    gals0=np.append(gals0,gals[index])
    mass_stellar=gals0['StellarMass']
    mass_virial=gals0['Mvir']

    #Defining axis of plot
    y1_stellar=mass_stellar
    y1_virial=mass_virial

    x1=(y1_virial)
    y1=np.log10((y1_stellar/y1_virial))

    Title='SHMR for z='+str(redshift)+' for Hen14 & Dev17 (MR & MRII)'
    plt.figure(redshift)
    plt.fig1=plt.semilogx(x1,y1,'r,',zorder=-1)
    plt.grid(True)
    plt.xlabel(Xlabel)
    plt.ylabel(Ylabel)
    plt.title(Title)

# Divide galaxies into mass bins using digitize function & average properties
    bins=np.digitize(y1_virial,bin_edges)

#Define percentile arrays -- create empty array for 
    low=np.empty(nbins)
    median=np.empty(nbins)
    high=np.empty(nbins)

# Loop over mass bins calculating properties
    for ibin in np.arange(1,len(bin_edges)):
        index=np.where(bins==ibin)[0]
        mvir_mean=np.mean(y1_virial[index]) #mean of viral mass
        mratio_mean=np.mean(y1[index]) #mean of stellar mass/viral mass
        std1=np.std(y1[index]) #standard deviation of y axis (ratio of stellar/viral)
        plt.plot(mvir_mean,mratio_mean,marker='o',color='g',linestyle='-')
#        plt.errorbar(0.5*(bin_edges[ibin-1]+bin_edges[ibin]),mratio_mean,yerr=std1,fmt='o-',capthick=2,color='lime')#,label='Dev17 - Mean & Std')
#        if len(index) > 0:
#            low[ibin-1],median[ibin-1],high[ibin-1]=np.percentile(y1[index],[16,50,84])
#    plt.errorbar(0.5*(bin_edges[ibin-1]+bin_edges[ibin]),mratio_mean,yerr=std1,fmt='o',capthick=2,color='g',label='Mean+stdDev')

#Median & percentiles
#    plt.plot(0.5*(bin_edges[0:nbins]+bin_edges[1:nbins+1]),low,':c')
#    plt.plot(0.5*(bin_edges[0:nbins]+bin_edges[1:nbins+1]),median,'-c',label='median')
#    plt.plot(0.5*(bin_edges[0:nbins]+bin_edges[1:nbins+1]),high,':c',label='16 & 84 percentiles')

#Plot
    sns.set_context('notebook')    # notebook, talk, paper
    sns.set_style('darkgrid')

    plt.legend(loc=1,numpoints=1)
    plt.ylim(np.log10(ylimits))

    plt.show()
    plt.savefig('SHMR'+str(redshift)+'.png')