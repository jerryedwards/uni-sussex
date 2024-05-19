#Plotting specific star formation rate (sSFR) for different red shifts

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import zoom
from astropy.table import Table
from plots_input import *

#MR & MRII datasets for z=0.00,1.04,2.07,3.11
boxside_MR=480.8 #In Mpc/h
datafile0_MR='data/Hen14_sfh2/snap_58.pkl'
datafile1_MR='data/Hen14_sfh2/snap_38.pkl'
datafile2_MR='data/Hen14_sfh2/snap_30.pkl'
datafile3_MR='data/Hen14_sfh2/snap_25.pkl'

datafile0_MRII='data/Hen14_MRII/snap_62.pkl'
datafile1_MRII='data/Hen14_MRII/snap_42.pkl'
datafile2_MRII='data/Hen14_MRII/snap_34.pkl'
datafile3_MRII='data/Hen14_MRII/snap_29.pkl'

#Loop to run for all redshifts
#redshifts=[0,1,2,3] 
redshifts=2,

# Active sSFR cuts (log_10(SFR/Mstellar)>sSFR_?)
sSFR_0=-11
sSFR_1=-10
sSFR_2=-9.7
sSFR_3=-9.5

#Variables used for the graph
Xlabel=r'$\log_{10}($M$_{*})$/$h^{-2}$M$_\odot$'
Ylabel='$\log_{10}$(SFR)/$h^{-2}$M$_\odot$ yr$^{-1}$'
pngfile='a_ratio_graph.png'
ylimits=[-2,4]
xlimits=[7.5,11.5]

#Plot Options
Datadir = '/home/j/je/je235/project/Obsdata/'
Hubble_h_WMAP7  = 0.704
Hubble_h_Planck = 0.673  

#Variables for Median & Error Bars
binedge  = np.linspace(8.8,11.8,10) 
NbinsErr = len(binedge)-1

#------------------------------------------------------------------------

#DEFINE LABELS

def plot_label (subplot,label_type, xlim, ylim, x_percentage, y_percentage, color, 
                x2_percentage=0., xlog=0, ylog=0, label='', linestyle='-', linewidth=2, 
                fontsize=16, fontweight='normal', sym='o', sym_size=5, err_size=0.1):    
    
     if xlog==0 & ylog==0:
      
         if label_type=='label':
             x=xlim[0]+(xlim[1]-xlim[0])*x_percentage
             y=ylim[0]+(ylim[1]-ylim[0])*y_percentage             
             plt.text(x,y,label, fontsize=fontsize, fontweight=fontweight)
         else:
             if label_type =='line':
                 x1=xlim[0]+(xlim[1]-xlim[0])*x_percentage
                 x2=xlim[0]+(xlim[1]-xlim[0])*x2_percentage
                 y=ylim[0]+(ylim[1]-ylim[0])*y_percentage
                 plt.plot([x1,x2],[y,y],color=color,linestyle=linestyle, linewidth=2)
                
             else:
                 if label_type=='symbol':
                     x=xlim[0]+(xlim[1]-xlim[0])*x_percentage
                     y=ylim[0]+(ylim[1]-ylim[0])*y_percentage                     
                     plt.errorbar(x,y, yerr=err_size, fmt=sym, markersize=sym_size, color=color)

#INPUT OBSERVATION DATA
def plot_obs(redshift):

    xlim=[8,12]
    ylim=[-2, 4]   

# #OBSERVATIONS
#values at all_z

#ELBAZ2007
    obs_slope_elbaz2007 =[0.77, -99.0, 0.9, -99.0, -99.0]
    obs_offset_elbaz2007=[np.log10(8.7)-(0.77*11.), -99.0, np.log10(7.2)-9, -99.0, -99.0]
    obs_offset_low_elbaz2007 =[np.log10(5.0)-(0.77*11.), -99.0, np.log10(3.6)-9, -99.0, -99.0]
    obs_offset_high_elbaz2007=[np.log10(16.1)-(0.77*11.), -99.0, np.log10(14.4)-9, -99.0, -99.0]

#KARIM2011
    file = Datadir + 'karim2011_sfr_mass_sf.txt'       
    karim2011 = Table.read(file, format='ascii') 
    karim_low_z_limit        = karim2011['col4']
    karim_medium_mass        = karim2011['col3']
    karim_sfr                = karim2011['col19']
    karim_sfr_error_up   = karim2011['col20']
    karim_sfr_error_down = karim2011['col21']
    log_karim_sfr_error_up=np.log10((karim_sfr+karim_sfr_error_up)/karim_sfr)
    log_karim_sfr_error_down=np.log10(karim_sfr/(karim_sfr-karim_sfr_error_down))

    obs_x=np.arange(xlim[0], xlim[1], 0.01)           

    if redshift==0.0:  
#ELBAZ2007
        obs_y=obs_x*obs_slope_elbaz2007[0] + obs_offset_elbaz2007[0] + 2.*np.log10(Hubble_h_WMAP7)
        plt.plot(obs_x+2.*np.log10(Hubble_h_WMAP7), obs_y, color='firebrick', linewidth=2,label='Elbaz et al. (2007)')            
        obs_y=obs_x*obs_slope_elbaz2007[0] + obs_offset_low_elbaz2007[0] + 2.*np.log10(Hubble_h_WMAP7)
        plt.plot(obs_x+2.*np.log10(Hubble_h_WMAP7), obs_y, color='firebrick', linewidth=2, linestyle='--')
        obs_y=obs_x*obs_slope_elbaz2007[0] + obs_offset_high_elbaz2007[0] + 2.*np.log10(Hubble_h_WMAP7)
        plt.plot(obs_x+2.*np.log10(Hubble_h_WMAP7), obs_y, color='firebrick', linewidth=2, linestyle='--')

    if redshift==0.4: 
#KARIM2011
        sel=(karim_low_z_limit==0.2) & (karim_medium_mass>8.8)                
        plt.plot.errorbar(karim_medium_mass[sel]+2.*np.log10(Hubble_h_WMAP7), 
                         np.log10(karim_sfr[sel])+2.*np.log10(Hubble_h_WMAP7), 
                         [log_karim_sfr_error_down[sel], log_karim_sfr_error_up[sel]], 
                         mfc='white', markeredgecolor='limegreen', color='limegreen', fmt='o', markersize=5)
        sel=(karim_low_z_limit==0.4) & (karim_medium_mass>8.9)                
        plt.plot.errorbar(karim_medium_mass[sel]+2.*np.log10(Hubble_h_WMAP7), 
                         np.log10(karim_sfr[sel])+2.*np.log10(Hubble_h_WMAP7), 
                         [log_karim_sfr_error_down[sel], log_karim_sfr_error_up[sel]], 
                         color='limegreen', fmt='o', markersize=5,label='Karim (2011)')

    if redshift==1.0:  
#ELBAZ2007
        obs_y=obs_x*obs_slope_elbaz2007[2] + obs_offset_elbaz2007[2] + 2.*np.log10(Hubble_h_WMAP7)
        plt.plot(obs_x+2.*np.log10(Hubble_h_WMAP7), obs_y, color='firebrick', linewidth=2)            
        obs_y=obs_x*obs_slope_elbaz2007[2] + obs_offset_low_elbaz2007[2] + 2.*np.log10(Hubble_h_WMAP7)
        plt.plot(obs_x+2.*np.log10(Hubble_h_WMAP7), obs_y, color='firebrick', linewidth=2, linestyle='--')            
        obs_y=obs_x*obs_slope_elbaz2007[2] + obs_offset_high_elbaz2007[2] + 2.*np.log10(Hubble_h_WMAP7)
        plt.plot(obs_x+2.*np.log10(Hubble_h_WMAP7), obs_y, color='firebrick', linewidth=2, linestyle='--')

#KARIM2011
        sel=(karim_low_z_limit==0.8) & (karim_medium_mass>9.1)                
        plt.errorbar(karim_medium_mass[sel]+2.*np.log10(Hubble_h_WMAP7), 
                         np.log10(karim_sfr[sel])+2.*np.log10(Hubble_h_WMAP7), 
                         [log_karim_sfr_error_down[sel], log_karim_sfr_error_up[sel]],
                         mfc='white', markeredgecolor='limegreen', color='limegreen', fmt='o', markersize=5)
        sel=(karim_low_z_limit==1.0) & (karim_medium_mass>9.3)                
        plt.errorbar(karim_medium_mass[sel]+2.*np.log10(Hubble_h_WMAP7), 
                         np.log10(karim_sfr[sel])+2.*np.log10(Hubble_h_WMAP7), 
                         [log_karim_sfr_error_down[sel], log_karim_sfr_error_up[sel]],  
                         color='limegreen', fmt='o', markersize=5,label='Karim (2011)')

#Whitaker2013
        file = Datadir + 'whitaker2013_mass_vs_sfr_z0.5_1.0.txt'       
        obs = Table.read(file, format='ascii') 
        plt.errorbar(obs['col1']+2.*np.log10(Hubble_h_WMAP7), obs['col2']+2.*np.log10(Hubble_h_WMAP7),
                         obs['col3'], mfc='white', markeredgecolor='blue', color='blue', fmt='o', markersize=5)
        file = Datadir + 'whitaker2013_mass_vs_sfr_z0.5_1.0.txt'       
        obs = Table.read(file, format='ascii') 
        plt.errorbar(obs['col1']+2.*np.log10(Hubble_h_WMAP7), obs['col2']+2.*np.log10(Hubble_h_WMAP7),
                         obs['col3'], color='blue', fmt='o', markersize=5,label='Whitaker (2013)')


    if redshift==2.0:  
#KARIM2011
        sel=(karim_low_z_limit==1.6) & (karim_medium_mass>9.6)                
        plt.errorbar(karim_medium_mass[sel]+2.*np.log10(Hubble_h_WMAP7), 
                         np.log10(karim_sfr[sel])+2.*np.log10(Hubble_h_WMAP7), 
                         [log_karim_sfr_error_down[sel], log_karim_sfr_error_up[sel]],
                         mfc='white', markeredgecolor='limegreen', color='limegreen', fmt='o', markersize=5,label='Karim (2011) - 1.6<z<2.0')
        sel=(karim_low_z_limit==2.0) & (karim_medium_mass>9.8)                
        plt.errorbar(karim_medium_mass[sel]+2.*np.log10(Hubble_h_WMAP7), 
                         np.log10(karim_sfr[sel])+2.*np.log10(Hubble_h_WMAP7), 
                         [log_karim_sfr_error_down[sel], log_karim_sfr_error_up[sel]],  
                         color='limegreen', fmt='o', markersize=5,label='Karim (2011) - 2.0<z<2.5')

#Whitaker2013
        file = Datadir + 'whitaker2013_mass_vs_sfr_z1.5_2.0.txt'       
        obs = Table.read(file, format='ascii') 
        plt.errorbar(obs['col1']+2.*np.log10(Hubble_h_WMAP7), obs['col2']+2.*np.log10(Hubble_h_WMAP7),
                         obs['col3'], mfc='white', markeredgecolor='blue', color='blue', fmt='o', markersize=5,label='Whitaker (2013) - 1.5<z<2.0')
        file = Datadir + 'whitaker2013_mass_vs_sfr_z2.0_2.5.txt'       
        obs = Table.read(file, format='ascii') 
        plt.errorbar(obs['col1']+2.*np.log10(Hubble_h_WMAP7), obs['col2']+2.*np.log10(Hubble_h_WMAP7),
                         obs['col3'], color='blue', fmt='o', markersize=5,label='Whitaker (2013) - 20.<z<2.5')

    if redshift==3.0:
#KARIM2011
        sel=(karim_low_z_limit==2.5) & (karim_medium_mass>10.0)                
        plt.errorbar(karim_medium_mass[sel]+2.*np.log10(Hubble_h_WMAP7), 
                         np.log10(karim_sfr[sel])+2.*np.log10(Hubble_h_WMAP7), 
                         [log_karim_sfr_error_down[sel], log_karim_sfr_error_up[sel]],
                         color='limegreen', fmt='o', markersize=5,label='Karim (2011) - 2.5<z<3.0')

#labels

    if redshift==0:
        plot_label (plt.plot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.82, 
                    color='black', xlog=0, ylog=0, label='Karim 2011', 
                    fontsize=13, fontweight='normal') 
        plot_label (plt.plot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.845, 
                    color='limegreen', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15) 
        
        plot_label (plt.plot, 'label', xlim, ylim, x_percentage=0.075, y_percentage=0.9, 
                    color='black', xlog=0, ylog=0, label='Whitaker 2013', 
                    fontsize=13, fontweight='normal') 
        plot_label (plt.plot, 'symbol', xlim, ylim, x_percentage=0.05, y_percentage=0.925, 
                    color='blue', xlog=0, ylog=0, sym='o', sym_size=5, err_size=0.15) 
        
        plot_label (plt.plot, 'label', xlim, ylim, 
                    x_percentage=0.15, y_percentage=0.74, color='black', xlog=0, ylog=0, 
                    label='Elbaz 2007', fontsize=13, fontweight='normal') 
        plot_label (plt.plot, 'line', xlim, ylim,
                    x_percentage=0.04, y_percentage=0.76, color='firebrick', x2_percentage=0.13, 
                    xlog=0, ylog=0, linestyle='-', linewidth=2)

#------------------------------------------------------------------------

# Loop over redshifts
for redshift in redshifts:

#Load in galaxy data -- SFR & Stellar Mass for MR -- pickling for snap_?
    pickleFile=eval('datafile'+str(redshift)+'_MR')
    pickleFile=eval('datafile'+str(redshift)+'_MRII')
    exec(open('script_unpickle.py').read())
    index=np.where(gals['Type']==0)
    gals=gals[index]
    mass_stellar=gals['StellarMass']
    Sfr=gals['Sfr']

#Define axis of plot
    x=np.log10(mass_stellar*Hubble_h_Planck)+10
    y=np.log10(Sfr*Hubble_h_Planck**2)

#Plot 1 for z=0.00,1.04,2.07,3.11 (MR data)
#plt.interactive(True) #For contour
    plt.figure(redshift)
    plt.close(redshift)
    plt.clf()
    plt.fig1=plt.plot(x,y,'r,',zorder=-1)
# Contouring code -- split data using histogram
    xlim=[7.5,11.5]
    ylim=[-2, 4]   
    bin=[0.1,0.05]
    Nbins=[int((xlim[1]-xlim[0])/bin[0]),int((ylim[1]-ylim[0])/bin[1])]
    Nlevel=20
    #Construct 2D histogram with variable bin width -- define bin edges:
    Ngals=len(x)
    H, xedges, yedges = np.histogram2d(x, y, bins=Nbins, range=[xlim,ylim])         
    extent = [xedges[0], xedges[-1],yedges[0], yedges[-1]]       
    mylevels = np.linspace(1., Nlevel, Nlevel)*0.7*Ngals/(Nbins[0]*Nbins[1])        
    #H = zoom(H, 20)        
    cont=plt.contourf(H.transpose()[::], cmap='Greys_r', levels=mylevels, extent=extent)  

    plt.grid(True)
    plt.xlabel(Xlabel)
    plt.ylabel(Ylabel)
    Title='SFR for z='+str(redshift)+'(MR data)'
    plt.title(Title)

#Divide galacies into mass bins -- digitize function
    bins=np.digitize(x,binedge)

#Define percentile arrays -- create empty array for 
    low    = np.empty(NbinsErr)
    median = np.empty(NbinsErr)
    high   = np.empty(NbinsErr)

# Loop over mass bins calculating properties
    for ibin in np.arange(1,len(binedge)):
        index= np.where(((y-x)>eval('sSFR_'+str(redshift)))&
                        (bins==ibin))[0] # For z=0 - needs ot be redshift dependent
        meanmass=np.mean(x[index])
        sfrmean=np.mean(y[index])
        std=np.std(y[index])
        plt.plot(meanmass,sfrmean,'b^')
        plt.errorbar(0.5*(binedge[ibin-1]+binedge[ibin]),sfrmean,yerr=std,fmt='o',capthick=2,color='y')
        if len(index) > 0:
            low[ibin-1],median[ibin-1],high[ibin-1]=np.percentile(y[index],[16,50,84])
    plt.errorbar(0.5*(binedge[ibin-1]+binedge[ibin]),sfrmean,yerr=std,fmt='o',capthick=2,color='y',label='Mean+stdDev')

#Median & percentiles
    plt.plot(0.5*(binedge[0:NbinsErr]+binedge[1:NbinsErr+1]),low,':y')
    plt.plot(0.5*(binedge[0:NbinsErr]+binedge[1:NbinsErr+1]),median,'-y',label='median')
    plt.plot(0.5*(binedge[0:NbinsErr]+binedge[1:NbinsErr+1]),high,':y',label='16 & 84 percentiles')

#Plot active-passive cut-off
    a=np.arange(8,12,0.01)              # Limits of plot - x-values
    b=a+(eval('sSFR_'+str(redshift)))   #y=mx+c where m=1
    plt.plot(a,b,'-c',label='Cut-off')

    plt.legend(frameon=False,fancybox=True)   #??????? Remove boarder & Reduce font size ????????

    plot_obs(redshift)
    plt.tight_layout()                                            
    plt.legend(loc=2)
    plt.ylim(ylimits)
    plt.xlim(xlimits)
    plt.show()
    plt.savefig('ssfrz'+str(redshift)+'.png')