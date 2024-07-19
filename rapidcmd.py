# -*- coding: utf-8 -*-
"""
Created on Fri May 31 15:27:14 2024

@author: hgtra
"""

import seaborn as sns
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.coordinates import SkyCoord
import numpy as np
import os
from scipy.stats import gaussian_kde
from matplotlib import rcParams

rcParams.update({'font.size':10})
pkde = 1
target = 13.73004, -37.60999
gtarget = "NGC 300"
target_name = "MLT 28037995"
dircats = "../../../Downloads/Catalogs/"
dirhsc = "../../../Downloads/Catalogs/HSCv3/"
dirmast = "../../../Downloads/MAST_Catalogs/"
verbose = 0

gals = Table.read(dircats+"galaxy_summary_new2.fits")
gnames = [' '.join(g.replace('NAME','').split()) for g in gals['main_id']]
gals['galaxy'] = gnames
Cands = Table.read(dircats+"YSG_candidates_HSC_280524.fits")
if target:
    c1 = SkyCoord(target[0],target[1],frame='icrs',unit='deg')
    c2 = SkyCoord(Cands['MatchRA'],Cands['MatchDec'],frame='icrs',unit='deg')
    d2d = c2.separation(c1)
    if d2d.arcsec.min()<1:
        target = Cands[d2d.deg.argmin()]
    else:
        target = 0
    
gals2 = Table.read(dircats+"galaxies_with_candidates_280524.csv")

Avfile = open("Av_stds.csv","w")

mapfilt = {}
filts = ['A_F475W','A_F555W','A_F606W','W3_F475W','W3_F555W','W3_F606W','W2_F555W','W2_F606W']

for filt in filts:
    mapfilt[filt]=filt
    mapfilt[filt[:5]]=filt

for j in range(len(gals2)):
    gname = gals2['galaxy'][j]
    if gtarget and gname!=gtarget:
        continue
    
    gal = gals[gals['galaxy']==gname]
    
    cands = Cands[Cands['galaxy']==gname] # let's remove some faint sources for visualization
    size = 4*np.log10(len(cands))-10
    gal2 = gals2[gals2['galaxy']==gname]
    if len(gal2)>=1:
        gal2 = gal2[0]
    else:
        continue
    cands = cands[(cands['y']-cands['y'].max())**2>size*np.random.random(len(cands))]
        
    
    try:
        # GALAXY PARAMETERS
        source = gal2['source']
        if gname=="M 31" and 0:
            source="HSCv3"
        if gname=="NGC 3310":
            source = "MAST"
            
        gcol = mapfilt[gal2['filter']]
        inst = gcol.split("_")[0]
        icol = inst+"_F814W"
        Av = gal2['Av']
        metal = gal2['metal']
        mstarref = gal['mstarref'][0]
        if gname=="WLM Galaxy":
            mstarref = 3
        
        # PARAMETERS OF MIST AND EXTINCTION+APERTURE CORRECTIONS
        metals = np.array([-1.2,-0.9,-0.6,-0.3])
        Avs = np.array([0,0.2,0.4,0.6,0.8,1,2,2.5])
        metal = metals[np.abs(metals-metal).argmin()]
        Avref = Avs[np.abs(Avs-Av).argmin()]
        prefix = {"A":"ACS_WFC","W3":"WFC3_UVIS","W2":"WFPC2"}
        cam = prefix[inst]
        dirmist = "../../../Downloads/MIST/MIST_%.1f_%.1f_%s/"%(metal,Avref,cam.split('_')[0])
        
        Alambda = (0.72e3/float(gcol.split('_')[-1][1:-1])-0.31)*Av
        Rlambda = 1/(0.72e3*(1/float(gcol.split('_')[-1][1:-1])-1/814))
        
        EBV = Av/Rlambda
        apcorr= {'A_F435W':-0.28,'A_F475W':-0.24,'A_F555W':-0.25,'A_F606W':-0.25,'A_F814W':-0.29,
               'W2_F439W':-0.21,'W2_F555W':-0.21,'W2_F606W':-0.23,'W2_F814W':-0.24,
               'W3_F336W':-0.24,'W3_F438W':-0.21,'W3_F475W':-0.2,'W3_F555W':-0.2,'W3_F606W':-0.2,'W3_F814W':-0.24}
        vegaconv = {'A_F435W':25.779-25.673,'A_F475W':26.168-26.068,'A_F555W':25.724-25.718,'A_F606W':26.398-26.486,'A_F814W':25.501-25.937,
                    'W3_F336W':23.46-24.64,'W3_F438W':24.98-24.83,'W3_F475W':25.79-25.69,'W3_F555W':25.81-25.78,'W3_F606W':25.99-26.08,'W3_F814W':24.67-25.09,
                    'W2_F439W':24.98-24.83,'W2_F555W':25.81-25.78,'W2_F606W':25.99-26.08,'W2_F814W':24.67-25.09}
        
        
        # LOAD CATALOG
        if source=="HSCv3":
            catname = dirhsc+"%s_data.fits"%gname
            if os.path.isfile(catname[:-1]):
                flatten = 1
                data = Table.read(catname[:-1])
            elif os.path.isfile(catname):
                data = Table.read(catname)
            else:
                continue
                
        elif source=="MAST":
            catname = dirmast+"%s_data.fits"%gname
            if os.path.isfile(catname):
                data = Table.read(catname)
            else:
                continue
    
    
    
        plt.close('all')
        x = np.asarray(data[gcol])-np.asarray(data[icol])-EBV
        x += vegaconv[gcol]+apcorr[gcol]-(vegaconv[icol]+apcorr[icol])
        y = np.asarray(data[gcol])-Alambda-5*np.log10(gal2['distance']*1e5)
        y += vegaconv[icol]+apcorr[icol]
        
        data['x'] = x
        data['y'] = y
        
        ind = ~np.isnan(x)
        data = data[ind.flatten()]
        lendata = len(data)
        if target_name != "2018mmb":
            data = data[np.argsort(np.random.random(lendata))[:50000]]
        iplot = np.argsort(np.random.random(len(data)))[:5000]
        for col in data.colnames:
            data[col] = data[col].flatten()
        x, y = data['x'], data['y']
        
        # PLOT CMD
        if verbose:print(gname,"plotting CMD of %d/%d sources"%(len(data[iplot]),lendata))
        plt.figure(1)
        sns.kdeplot(data[iplot].to_pandas(),x='x',y='y',fill=True, cmap='viridis_r')
        plt.scatter(cands['x'], cands['y'], s=max(1,int(10-2*size)), color='red', marker='.', label='candidates')
        if target:
            plt.scatter(target['x'], target['y'], s=50, color='k', marker='*', label=target_name)
        if target_name == "2018mmb":
            # plot the better, brighter match having MatchID=33625884
            print(data[data['MatchID']==33625884]['x'],data[data['MatchID']==33625884]['y'])
            plt.scatter(data[data['MatchID']==33625884]['x'],data[data['MatchID']==33625884]['y'], s=50, color='k', marker='*', label="2018mmb (alt)")
        plt.xlabel("%s - %s"%(gcol,icol))
        plt.ylabel("%s_abs"%(gcol))
        
        
        # PLOT MIST CURVES
        colors = {3:'gray', 4:'m',6:'k',8:'cyan',10:'C2',12:'C1',16:'y',20:'lime'}
        colors[mstarref] = 'blue'
        ls = ['--','-.','-',':']
        gcol = prefix[inst]+gcol[-6:]
        icol = prefix[inst]+icol[-6:]
        
        deltax = np.array([0,0.4/Rlambda])
        deltay = np.array([0,Alambda/Av*0.4])
        for mstar in np.unique([mstarref,4,6,8,10,12,16,20]):
            mist = np.genfromtxt(dirmist+"/0%03d000M.track.eep.cmd"%(mstar*10),names=True,skip_header=14)
            mist = mist[mist['phase']<=4]
            if mstar<10:
                mist = mist[mist['phase']<=2]
            mist = mist[mist['phase']>=0]
            
            xm, ym = mist[gcol]-mist[icol], mist[gcol]
            
            mid = 0.8*(max(xm)+min(xm))
            wid = max(xm)-min(xm)
            # stop gap at faintest magnitude of the red part
            gapstop = np.argmax(xm>mid)+np.argmax(ym[xm>mid])
            # start gap at 0.05 mag before brightest magnitude
            gapstart = mist[gcol][:gapstop].argmin()
            gapstart = np.argmin(np.abs(mist[gcol][:gapstart]-mist[gcol][gapstart]-0.05))
            if mstar==mstarref+0:
                xref, yref = xm[gapstop],np.mean(ym[(gapstart,gapstop),])
            plt.scatter(xm[(gapstart,gapstop),]+deltax,ym[(gapstart,gapstop),]+deltay,color=colors[mstar],s=5)
            for phase in range(0,3):
                ind = mist['phase']==phase
                if phase==2:
                    ind = mist['phase']>=phase
                    plt.plot(xm[ind],ym[ind], ls=ls[phase],lw=3,color=colors[mstar],label=r"%d M$_\odot$"%mstar)
                else:
                    plt.plot(xm[ind],ym[ind], ls=ls[phase],lw=3,color=colors[mstar])
                    
                    
        rgbquadrant = np.logical_and(x>xref,y<1e-1/(xref-x)+yref)
        
        xrgb,yrgb = x[rgbquadrant], y[rgbquadrant]
        from numpy.polynomial.polynomial import Polynomial
        # Fit two polynomials of degree 2, y=f(x) and x=f(y)
        p1 = Polynomial.fit(xrgb, yrgb, 1)
        a1 = p1.deriv().coef[0]
        p2 = Polynomial.fit(yrgb, xrgb, 1)
        a2 = 1/p2.deriv().coef[0]
        stds = []
        
        x_fit1 = np.linspace(min(xrgb), max(xrgb), 100)
        plt.plot(x_fit1,)
        y_fit2 = np.linspace(max(yrgb), min(yrgb), 100)
        y_fit1 = p1(x_fit1)
        x_fit2 = p2(y_fit2)
        x_fit12, y_fit12 = (x_fit1+x_fit2)/2, (y_fit1+y_fit2)/2
        p12 = Polynomial.fit(x_fit12, y_fit12, 1)
        
        plt.plot(x_fit12, p12(x_fit12),color='C3',ls=':', label='RGB')
        
        a = p12.deriv().coef[0]
        
        plt.xlim((np.percentile(x,0.1)-1.5,
                      np.percentile(x,99.9)-0.3))
        plt.ylim((np.percentile(y,99.9)+0.5,
                      np.percentile(y,0.1)-1.5))
        plt.legend(frameon=False,loc=2)
        
        
        
        
        if target:
            suffix = "_%s"%target_name
        else:
            suffix = ""
        plt.savefig("../figures/%s_cmd_light%s.png"%(gname,suffix))
        
        # PLOT Av SCATTER
        data = data[rgbquadrant]
        x, y = xrgb, yrgb
        
        # Rotate the data
        x_translated = x - x.mean()
        y_translated = y - y.mean()
        theta = -np.arctan(Rlambda)
        
        rotation_matrix = np.array([
            [np.cos(theta), -np.sin(theta)],
            [np.sin(theta), np.cos(theta)]])
        data_translated = np.vstack((x_translated, y_translated))
        data_rotated = rotation_matrix @ data_translated
        x_new, y_new = data_rotated[0, :], data_rotated[1, :]
        
        # correct slope to de-bias x
        atest = np.linspace(-5,5,100)
        for asl in atest:
            X = (x_new-y_new/asl)
            X = X[np.percentile(X,2)<X]
            X = X[X<np.percentile(X,98)]
            stds.append(X.std())
            
        imin = np.argmin(stds)
        plt.figure()
        plt.plot(atest,stds-stds[imin]+1e-3,color='k')
        plt.axvline(atest[imin],color='r')
        plt.xlabel('slope')
        plt.ylabel('std')
        if verbose:print("Minimum std: %.3f (%d) a=%.2f"%(stds[imin],imin,atest[imin]))
        plt.gca().set_yscale('log')
        x_new -= y_new/atest[imin]
        ind = np.logical_and(np.percentile(x_new,2)<x_new,x_new<np.percentile(x_new,98))
        data, x_new, y_new = data[ind], x_new[ind], y_new[ind]
        
        data['a'] = x_new
        data['b'] = y_new
        if verbose:print("a=%.2f, std=%.3f"%(a,x_new.std()))
        x_new*= abs(np.sin(theta)) # Av '+' E(B-V) to Av
        Avsig = x_new.std()
        
        if pkde:
            sns.jointplot(data.to_pandas(),x='a',y='b', cmap='viridis_r',kind='kde')
        else:
            sns.jointplot(data.to_pandas(),x='a',y='b', cmap='viridis_r')
        
        plt.axvline(np.median(x_new)+Avsig,color='k',ls=':')
        plt.axvline(np.median(x_new)-Avsig,color='k',ls=':')
        plt.gca().invert_yaxis()
        plt.xlabel("A$_V$ scatter in the RGB (%.2f)"%Avsig,fontsize=15)
        plt.ylabel('Magnitude',fontsize=15)
        plt.tight_layout()
        if pkde:
            plt.savefig("../figures/%s_ebv_light.png"%(gname))
        else:
            plt.savefig("../figures/%s_Av_light.png"%(gname))
        print(gname,Avsig,file=Avfile)
        
        #aleft = slope gapstart points
        #cleft = offset gapstart points
        #cAv = cands['y']-Rlambda*(Alambda/Av)*cands['x']
        #xinter,yinter = (cAv-cleft)/(aleft-Rlambda*(Alambda/Av)), (Rlambda*(Alambda/Av)*cleft-aleft*cAv)/(Rlambda*(Alambda/Av)-aleft)
        #Avsig_corr = np.sqrt(max(Avsig**2-0.1**2,1e-4))
        #cands['Pextin'] = 2*(1-norm.cdf(x[cand]-xinter,0,Avsig_corr/Rlambda))
    except Exception as e:
        print(gname,e)

Avfile.close()