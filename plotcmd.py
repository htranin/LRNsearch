# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 23:30:21 2023

@author: hgtra
"""

from astropy.table import Table
import numpy as np
import matplotlib.pyplot as plt
from fastkde import fastKDE
from scipy.interpolate import RectBivariateSpline, LinearNDInterpolator
from scipy.stats import norm
import os
from astropy.coordinates import SkyCoord
#import sys
from sklearn.mixture import GaussianMixture
from warnings import filterwarnings
filterwarnings('ignore')
import json


# user defined parameters
target = "M 31"
itarget = 0
forcefilt = ""#"ACS_F606W"#"WFC3_F438W"
forcecoord = 0 #190.582, 32.544
beam = 0  # pencil beam (beam)*arcmin at random
dmetal = -0.3 # manual metallicity offset
dAv = 0.4# manual Av offset
icoly = 0 # 0 = green 1 = red
pct = 90  # pct of likelihood distribution to keep, along the gap of the reference MIST track
          # warning: if likelihood<2th percentile of gaussian mixture values, then contour is taken at the 2th percentile
factor = 6.5 # compression factor of y-axis for Gaussian Mixture
mstar0 = 4 # stellar mass of reference MIST track
mstar1 = 6 # stellar mass of auxiliary MIST track (slope computation)
reduce = 1 # maximum factor by which sample size is reduced, for GM and computing time purpose
rowlim = 100 # number of rows below which select MAST data rather than HSCv3
forcedist = 0 # force value of distance
fullplot = 1. # force xlim and ylim to show all data

apcorr= {'A_F435W':-0.28,'A_F475W':-0.24,'A_F555W':-0.25,'A_F606W':-0.25,'A_F814W':-0.29,
       'W2_F439W':-0.21,'W2_F555W':-0.21,'W2_F606W':-0.23,'W2_F814W':-0.24,
       'W3_F336W':-0.24,'W3_F438W':-0.21,'W3_F475W':-0.2,'W3_F555W':-0.2,'W3_F606W':-0.2,'W3_F814W':-0.24}
#https://archive.stsci.edu/hst/hsc/help/HSC_faq/ci_ap_cor_table_2016.txt
vegaconv = {'A_F435W':25.779-25.673,'A_F475W':26.168-26.068,'A_F555W':25.724-25.718,'A_F606W':26.398-26.486,'A_F814W':25.501-25.937,
            'W3_F336W':23.46-24.64,'W3_F438W':24.98-24.83,'W3_F475W':25.79-25.69,'W3_F555W':25.81-25.78,'W3_F606W':25.99-26.08,'W3_F814W':24.67-25.09,
            'W2_F439W':24.98-24.83,'W2_F555W':25.81-25.78,'W2_F606W':25.99-26.08,'W2_F814W':24.67-25.09}
#https://iopscience.iop.org/article/10.1086/444553/pdf
#https://www.stsci.edu/files/live/sites/www/files/home/hst/instrumentation/wfc3/documentation/instrument-science-reports-isrs/_documents/2009/WFC3-2009-31.pdf
#For WFPC2, the best documentation https://www.stsci.edu/files/live/sites/www/files/home/hst/documentation/_documents/wfpc2/wfpc2_dhb_v5.pdf
# does not give the AB zeropoints > I use those of WFC3... Should be ok since WFC3 ~ ACS conversions.

prefix = {"ACS":"A","WFC3":"W3","WFPC2":"W2"}
prefix2 = {"ACS":"ACS_WFC","WFC3":"WFC3_UVIS","WFPC2":"WFPC2"}
oprefix = {v: k for k, v in prefix.items()}
metals = np.array([-1.2,-0.9,-0.6,-0.3])
Avs = np.array([0,0.2,0.4,0.6,0.8,1,2,2.5])
#source = "HSCv3"

gals = Table.read("../../../Downloads/MAST_catalogs/galaxies_selected_filter_0923_5.fits")
gals = Table.read("../../../Downloads/Catalogs/galaxy_mast_stats_addendum_1123.fits")
if target==" ":
    target = " ".join(gals['main_id'][itarget].split())
gals2 = Table.read("../../../Downloads/Catalogs/galaxy_compilation_Simbad_unique_1023.fits")
gals2 = Table.read("../../../Downloads/Catalogs/galaxies_hecate_20Mpc_new.fits")
gals2['selected_filter']="           " # Let's be agnostic about non-MAST galaxies
gals2['selected_filter'][1:]=" "

gals = Table.read("../../../Downloads/Catalogs/galaxy_summary_new.fits")

def rotate(x,y,th):
    return np.sin(th*np.pi/180)*x+np.cos(th*np.pi/180)*y,np.sin(th*np.pi/180)*y-np.cos(th*np.pi/180)*x

EBV_std_dic = {}

Igal = 0
for gal in gals[:]:
    source = "MAST"
    tget = " ".join(gal['main_id'].replace("NAME","").split())
    if target!="" and tget!=target:
        continue
    print(tget)
    if target=="" and np.isin(['Av','dAv','metal','dmetal','pct','forcedist','mstarref'],gals.colnames).all():
        Av, metal, dAv, dmetal, pct, forcedist, mstar0 = gal['Av','metal','dAv','dmetal','pct','forcedist','mstarref']
        Av, metal = Av+dAv, metal+dmetal
        mstar1 = 6 if mstar0!=6 else 8
        if gal['R1']>1800:
            print('galaxy too large, skipping')
            continue
        if np.ma.is_masked(mstar0) and target=="":
            print('skipping galaxy: "%s"'%gal['commentaire'])
            continue
    else:
        #gal = gals[np.array([" ".join(g.replace("NAME","").split())==tget for g in gals["main_id"]])][0]
        gal2 = gals2[np.array([" ".join(g.replace("NAME","").split())==tget for g in gals2["main_id"]])][0]
        if 0 and gal2['R1']>1800 and not(beam):
            print("galaxy too large, skipping")
            continue
        metal = gal2["logMH_Guseva"]+dmetal # offset to account for stellar evolution
        Av = gal2["A0"]+dAv
        
    
    metal = metals[np.abs(metals-metal).argmin()]
    print(gal['main_id','Dist','selected_filter'], metal)    
    Avref = Avs[np.abs(Avs-Av).argmin()]
    Avmax = Avs[np.abs(Avs-Av-0.4).argmin()]
    
    if not(forcedist):
        mod = 5*np.log10(gal['Dist']*1e6)-5
    else:
        mod = 5*np.log10(forcedist*1e6)-5
    if source=="HSCv3":
        filename = "../../../Downloads/Catalogs/HSCv3/%s_data.fit"%tget
        if not(os.path.isfile(filename)):
            filename = filename.replace(".fit",".fits")
    else:
        filename = "../../../Downloads/MAST_catalogs/%s_data2.fits"%tget
        
    filtcols = ["A_F435W","W3_F438W","W2_F439W","A_F475W","W3_F475W","A_F555W","W3_F555W","W2_F555W","A_F606W","W3_F606W","W2_F606W"]
    if os.path.isfile(filename):
        data = Table.read(filename)
        col_v = [c for c in data.colnames if c[-1]=='W' and not('814' in c)][np.argmax([len(np.where(~np.isnan(data[c]))[0]) for c in data.colnames if c[-1]=='W' and not('814' in c)])]
        if len(np.where(~np.isnan(np.maximum(data[col_v],data[col_v.split('_')[0]+'_F814W'])))[0])<rowlim:
            try:
                data = Table.read("../../../Downloads/MAST_catalogs/%s_data.fits"%(tget))
                col_v = [c for c in data.colnames if c[-1]=='W' and not('814' in c)][np.argmax([len(np.where(~np.isnan(data[c]))[0]) for c in data.colnames if c[-1]=='W' and not('814' in c)])]
                if target!='':
                    forcefilt = "%s_%s"%(oprefix[col_v.split('_')[0]],col_v.split('_')[1])
                source = "MAST"
                print('HSCv3 catalog has less than %d sources, trying with MAST'%rowlim)
            except:
                print('Catalog has less than %d sources, skipping'%rowlim)
                continue
    
    
    
    
    if not(os.path.isfile(filename)) or not(np.in1d(filtcols,data.colnames).all()) and len(data)<1e6:
        # no file or deprecated file, need to download again
        import casjobs_query
        print("running CasJobs query")
        casjobs_query.run_query(tget.replace(" ","").replace("NAME",""))
        galsfailed = np.genfromtxt("GALAXIES_TO_DOWNLOAD_ON_CASJOBS.txt",delimiter=",",encoding="utf-8",dtype=None)
        if len(galsfailed)>0 and galsfailed[-1]==tget:
            continue
    if not(os.path.isfile(filename)):
        # still no data after CasJobs query: not in HSCv3
        try:
            data = Table.read("../../../Downloads/MAST_catalogs/%s_data.fits"%(tget))
            col_v = [c for c in data.colnames if c[-1]=='W' and not('814' in c)][np.argmax([len(np.where(~np.isnan(data[c]))[0]) for c in data.colnames if c[-1]=='W' and not('814' in c)])]
            if target!='':
                forcefilt = "%s_%s"%(oprefix[col_v.split('_')[0]],col_v.split('_')[1])
            source = "MAST"
        except:
            print('Catalog neither in HSCv3 nor in MAST local data, skipping')
            continue
    elif source!="MAST":
        data = Table.read(filename)
        col_v = [c for c in data.colnames if c[-1]=='W' and not('814' in c)][np.argmax([len(np.where(~np.isnan(data[c]))[0]) for c in data.colnames if c[-1]=='W' and not('814' in c)])]
        
    Alambda = (0.72e3/float(col_v.split('_')[-1][1:-1])-0.31)*Av
    Rlambda = 1/(0.72e3*(1/float(col_v.split('_')[-1][1:-1])-1/814))
    EBV = Av/Rlambda
    
    if 1 or gal['selected_filter']==" ":
        print("Exposure is too shallow, skipping")
        if not(forcefilt):
            gal['selected_filter'] = oprefix[col_v.split('_')[0]]+'_'+col_v.split('_')[-1]
        #continue
    
    if forcefilt:
        gal['selected_filter'] = forcefilt
        
    inst = gal['selected_filter'].split("_")[0]
    col = prefix[inst]+"_"+gal['selected_filter'].split("_")[-1]
    col2 = prefix[inst]+"_F814W"
    colm = prefix2[inst]+"_"+gal['selected_filter'].split("_")[-1]
    colm2 = prefix2[inst]+"_F814W"
        
    
    
    
        
    if len(np.where(~np.isnan(data[col]))[0])<len(data)/4 and not(forcefilt):
        nsrc = np.zeros(len(filtcols))
        magl = 99*np.ones(len(filtcols))
        for i,filt in zip(range(len(filtcols)),filtcols):
            if int(filt[-4:-1])<int(col[-4:-1]):
                # exposure too shallow in bluer filters
                continue
            nsrc[i] = len(np.where(~np.isnan(data[filt]))[0])
            if nsrc[i]>0:
                #limiting magnitude in this filter
                magl[i] = np.nanpercentile(np.asarray(data[filt]),99)
        if not(nsrc.any()):
            # relax condition on exposure depth
            for i,filt in zip(range(len(filtcols)),filtcols):
                nsrc[i] = len(np.where(~np.isnan(data[filt]))[0])
                if nsrc[i]>0:
                    #limiting magnitude in this filter
                    magl[i] = np.nanpercentile(np.asarray(data[filt]),99)*606/int(filt[-4:-1])
            
        
        magl = magl[nsrc>max(nsrc)/4]
        filtcols2 = np.array(filtcols)[nsrc>max(nsrc)/4]
        col = filtcols2[magl.argmax()]
        #col =  filtcols[nsrc.argmax()]
        print("\t".join(["%s: %d"%(filt,i) for i,filt in zip(nsrc,filtcols)]))
        newfilt = "%s_%s"%({j:i for (i,j) in prefix.items()}[col[-8:-6]],col[-5:])
        print("old / new selected filter:",gal['selected_filter'],newfilt)
        
        gal['selected_filter'] = newfilt
        inst = gal['selected_filter'].split("_")[0]
        col = prefix[inst]+"_"+gal['selected_filter'].split("_")[-1]
        col2 = prefix[inst]+"_F814W"
        colm = prefix2[inst]+"_"+gal['selected_filter'].split("_")[-1]
        colm2 = prefix2[inst]+"_F814W"
    
    colmy = [colm,colm2][icoly]
    coly = [col,col2][icoly]
    #gal['selected_filter'] = "ACS_F555W"
    
    inst = gal['selected_filter'].split("_")[0]
    
    dirmist = "../../../Downloads/MIST/MIST_%.1f_%.1f_%s/"%(metal,Avref,inst)
    dirmistAv = "../../../Downloads/MIST/MIST_%.1f_%.1f_%s/"%(metal,Avmax,inst)
    plt.rcParams.update({'font.size': 20, "text.usetex": False})
    plt.figure(figsize=(12,8))
    

    
    
    
    print("Limiting mag in aperture",np.max(data[col]-mod))
    
    print("Unabsorbed limiting mag",np.max(data[col]+apcorr[col]-mod-Alambda))
    
    data[col] += vegaconv[col]+apcorr[col]
    
    print("Unabsorbed limiting Vega mag",np.max(data[col]-mod-Alambda))
    
    data[col2] += vegaconv[col2]+apcorr[col2]
    
    if len(data[col].shape)==2:
        data = data[~np.isnan(np.asarray(data[col][:,0]))]
        data = data[~np.isnan(np.asarray(data[col2][:,0]))]
        
        x = np.asarray(data[col]-data[col2]-EBV)[:,0]
        y = np.asarray(data[coly]-mod-Alambda)[:,0]
        # 0.6*Av = A_814 approximately
    else:
        data = data[~np.isnan(np.asarray(data[col]))]
        data = data[~np.isnan(np.asarray(data[col2]))]
        
        x = np.asarray(data[col]-data[col2]-EBV)
        y = np.asarray(data[coly]-mod-Alambda)
        # 0.6*Av = A_814 approximately
    if len(data)==0:
        print("empty red+green selection")
        continue
        
    ra, dec = data['MatchRA','MatchDec'][np.random.randint(len(data))]
    if forcecoord:
        ra, dec = forcecoord[0],forcecoord[1]
    if beam:
        ind = SkyCoord(np.reshape(data['MatchRA'],-1),np.reshape(data['MatchDec'],-1),frame='icrs',unit="deg").separation(SkyCoord(ra,dec,frame='icrs',unit="deg")).arcmin<beam
        data, x, y = data[ind], x[ind], y[ind]
    
    print("compute PDF and Gaussian Mixture")
    myPDF,axes = fastKDE.pdf(x,y,numPoints=2**9+1)
    finterpKDE = RectBivariateSpline(axes[1],axes[0],myPDF,s=3)
    zKDE = finterpKDE(y,x,grid=False)
    # alleviate overdense areas to better fit outskirts and speed computation
    
    zKDEcut = np.percentile(zKDE, max(100/reduce,300-55*np.log10(max(len(x),4.4e3))))
    
    rem = np.random.random(len(x))*zKDE>zKDEcut
    x,y,data, zKDE = x[~rem], y[~rem],data[~rem], zKDE[~rem]
    
    X = np.vstack([x, y/factor]).T
    # truncate for speed
    X = X[::5]
    
    ncomp = int(4*np.log10(len(x))*(1-70/len(x)))
    print("number of components",ncomp)
    if ncomp<=1:
        print("too few sources (%d)"%len(data))
        continue
    models = [GaussianMixture(n_components=ncomp, max_iter=1000,
                              covariance_type="full") for k in range(15)]
    BIC = []
    for k in range(len(models)):
        models[k].fit(X)
        for p, xcyc in enumerate(models[k].means_):
            if max(finterpKDE(xcyc[1]*factor,xcyc[0],grid=False),models[k].weights_[p])<1e-2 or abs(models[k].covariances_[p][0,1]/np.sqrt(models[k].covariances_[p].diagonal().prod()))>0.99:
                models[k].weights_[p] /= 100
                print("removing component in %d"%p)
            
        BIC.append(models[k].bic(X))
        
    model = models[np.argmin(BIC)]
    print("models are computed (best model %d)"%np.argmin(BIC))

    
    Xgrid = np.array(np.meshgrid(axes[0],axes[1]/factor)).reshape(2,-1).T
    myGM = 2**model.score_samples(Xgrid).T.reshape((len(axes[0]),len(axes[1])))
    finterp = RectBivariateSpline(axes[1]/factor,axes[0],myGM)
    
    
    # interpolate to get z values at points
    z = finterp(y/factor,x,grid=False)
    
    print("interpolation done")
    
    # Sort the points by density, so that the densest points are plotted last
    idx = z.argsort()
    xs, ys, zs = x[idx], y[idx], zKDE[idx]
    
    zmed = np.median(zs)
    plt.scatter(xs[zs>zmed], ys[zs>zmed], c=zs[zs>zmed], s=5+10*(int(len(x)<1e3)+int(len(x)<1e4)), cmap='viridis_r', vmin=zs.min(), vmax=zs.max())
    plt.scatter(xs[zs<zmed], ys[zs<zmed], c=zs[zs<zmed], s=3+10*(int(len(x)<1e3)+int(len(x)<1e4)), cmap='viridis_r', vmin=zs.min(), vmax=zs.max(), label="%s data"%source)
    
    #plt.contour(axes[0], axes[1], myPDF, levels=[np.percentile(myPDF[myPDF>0],70),np.percentile(myPDF[myPDF>0],75),np.percentile(myPDF[myPDF>0],80)])
    
    #plt.scatter(data[col]-data[col2]-EBV,data[col]-mod-Av,label="HSCv3 data")
    
    plt.title("%s, Z= %.1f Av= %.1f, %s"%(tget, metal,Avref,inst)+int(beam!=0)*(", around ra=%.3f dec=%.3f"%(ra,dec)))
    
    plt.xlabel("%s-%s-E(%s-F814W)"%(col,col2,col.split("_")[-1]))
    plt.ylabel("%s-A(%s)"%(coly,coly.split("_")[-1]))
    
    colors = {3:'gray', 4:'m',6:'k',8:'cyan',10:'C2',12:'C1',16:'y',20:'lime'}
    colors[mstar0] = 'blue'
    ls=['--','-.','-',':']
    
    mistx,misty,mistm,mistl,mistt=[],[],[],[],[]
    shallow = 0
    mstarlist = [4,6,8,10,12,16,20][::int((mstar1-mstar0)/abs(mstar1-mstar0))]
    if mstar0==3:
        mstarlist = [3]+mstarlist
    for mstar in mstarlist:
        print(mstar)
        mist = np.genfromtxt(dirmist+"/0%03d000M.track.eep.cmd"%(mstar*10),names=True,skip_header=14)
        mist = mist[mist['phase']<=4]
        mist = mist[mist['phase']>=0]
        
        xm, ym = mist[colm]-mist[colm2], mist[colmy]
        
        mid = 0.5*(max(xm)+min(xm))
        wid = max(xm)-min(xm)
        gapstop = min(np.where(ym==max(ym[xm>mid]))[0][0],len(ym)-1)
        gapstart = np.argmin(mist[colm][:gapstop]-0.5*mist[colm2][:gapstop])
        mistx = mistx+list(xm[gapstart:gapstop])
        misty = misty+list(ym[gapstart:gapstop])
        mistm = mistm+(gapstop-gapstart+1)*[mstar]
        mistl = mistl+list(mist[gapstart:gapstop+1]['log_L'])
        mistt = mistt+list(mist[gapstart:gapstop+1]['log_Teff'])
        
        if 1:
            mistAv = np.genfromtxt(dirmistAv+"/0%03d000M.track.eep.cmd"%(mstar*10),names=True,skip_header=14)
            mistAv = mistAv[mistAv['phase']<4]
            mistAv = mistAv[mistAv['phase']>=0]
            
            xmAv, ymAv = mistAv[colm]-mistAv[colm2], mistAv[colmy]
            
            midAv = 0.5*(max(xmAv)+min(xmAv))
            widAv = max(xmAv)-min(xmAv)
            gapstopAv = np.where(ymAv==max(ymAv[xmAv>midAv]))[0][0]
            gapstartAv = np.argmin(mistAv[colm][:gapstopAv]-0.5*mistAv[colm2][:gapstopAv])
            if mstar==mstar0:
                # assign last gap point
                xAv0,yAv0 = xmAv[gapstopAv],ymAv[gapstopAv]
            elif mstar==mstar1:
                # assign last gap point
                xAv2,yAv2 = xmAv[gapstopAv],ymAv[gapstopAv]
                
        # gross extrapolation to include extincted candidates
        mistx = mistx + [mistx[-1]+0.5]
        misty = misty + [misty[-1]+1.5]
        
        
        if mstar==mstar0:
            xmref, ymref = xm[gapstart:gapstop], ym[gapstart:gapstop]
            Mlim = max(ym[xm>mid])
            if Mlim>np.nanpercentile(y,95):
                print("exposure too shallow, sample not complete at %.1f solar masses"%mstar0)
                shallow = 1
            # compute likelihood along MIST track in the gap
            interp_mist = finterp(ym[gapstart:gapstop]/factor,xm[gapstart:gapstop],grid=False)
            minlevel = min(interp_mist)
            maxlevel = max(interp_mist)
            if maxlevel<np.percentile(z,2):
                medlevel = np.percentile(z,2)
            else:
                medlevel = np.percentile(interp_mist,pct)
            minlevel,maxlevel = medlevel-abs(maxlevel-medlevel)/2,medlevel+abs(maxlevel-medlevel)/2
            if maxlevel==medlevel:
                maxlevel=maxlevel+abs(maxlevel)/2.
            # plot likelihood contour crossing the MIST track in the gap
            
            


            # assign first gap point
            xs1,ys1 = xm[gapstart],ym[gapstart]
            # assign last point of MS
            xms1, yms1 = xm[mist['phase']==0][-1],ym[mist['phase']==0][-1]
            # assign last gap point
            x0,y0 = xm[gapstop],ym[gapstop]
            # assign last gap point in low-density area i.e. at medlevel
            #i0 = np.where(finterp_mist<medlevel)[0][-1] 
            #x0,y0 = xm[gapstart+i0],ym[gapstart+i0]
            #x0,y0 = x1,y1
           
        elif mstar==mstar1:
            # assign first gap point
            xs2,ys2 = xm[gapstart],ym[gapstart]
            # assign last point of MS
            xms2, yms2 = xm[mist['phase']==0][-1],ym[mist['phase']==0][-1]
            
            aleft=(ys2-ys1)/(xs2-xs1)
            cleft = ys2-aleft*xs2
            aright=(yAv2-yAv0)/(xAv2-xAv0)
            
            # compute x of the red selection line at magnitude -9
            xf = xAv0+(-9-yAv0)*(xAv2-xAv0)/(yAv2-yAv0)
            
            # compute x of the blue selection line at magnitude -9
            xsf = xs1+(-9-ys1)*(xs2-xs1)/(ys2-ys1) # parallel to gap starts
            #xsf = xs1+(-9-yms1)*(xms2-xms1)/(yms2-yms1) # parallel to MS
            
            
        #elif mstar==12:
            #x3,y3 = xm[gapstop],ym[gapstop]
            #xc = ((x3**2-x2**2+y3**2-y2**2)/(2*(y3-y2))-(x2**2-x1**2+y2**2-y1**2)/(2*(y2-y1)))/((x3-x2)/(y3-y2)-(x2-x1)/(y2-y1))
            #yc = -(x2-x1)/(y2-y1)*xc+(x2**2-x1**2+y2**2-y1**2)/(2*(y2-y1))
            #Rc = np.sqrt((x1-xc)**2+(y1-yc)**2)
            #xc += x0-x1 # translation to low-density region
            #yc += y0-y1 # translation to low-density region
            
            #plt.gca().add_patch(plt.Circle((xc,yc),Rc,color='cyan',fill=False))
            
        #print(mid,wid)
        
        
        
        #plot gap
        #plt.plot(xm[gapstart:gapstop],ym[gapstart:gapstop],color="lightgrey",alpha=0.5,lw=6)
        for phase in range(0,3):
            ind = mist['phase']==phase
            if phase==2 and mstar>=10:
                ind = np.logical_and(mist['phase']>=phase, np.arange(len(mist))<gapstop)
            
            
            if phase==0:
                plt.plot(xm[ind],ym[ind], ls=ls[phase],lw=3,color=colors[mstar],label=r"%d M$_\odot$"%mstar)
            else:
                plt.plot(xm[ind],ym[ind], ls=ls[phase],lw=3,color=colors[mstar])
    
    interpm = LinearNDInterpolator(list(zip(mistx, misty)), mistm)
    interpt = LinearNDInterpolator(list(zip(mistx, misty)), mistt)
    interpl = LinearNDInterpolator(list(zip(mistx, misty)), mistl)
    
    
    
    if shallow:
        plt.xlim((min(x)-0.5,max(x)+1))
        plt.ylim((max(y)+0.5,min(y)-1))
        plt.quiver(x0,y0,xAv0-x0,yAv0-y0,color='blue',scale_units='xy', angles='xy', scale=1)
        plt.plot([xAv0,xf],[yAv0,-9],color='blue', lw=3)
        plt.plot([xs1,xsf],[ys1,-9],color='blue', lw=3, label="selection")
        plt.legend(loc=1, fontsize=17, framealpha=1, edgecolor='white')
        plt.savefig("../figures/cmd_%s_%s_%s.pdf"%(tget,source,col))
        # exposure too shallow
        continue
    
    CS = plt.contour(axes[0], axes[1], myGM, levels=[medlevel,maxlevel], cmap="coolwarm", vmax=maxlevel+(maxlevel-medlevel)/4, linewidths=3, linestyles="dotted")
    pcs = CS.collections[0].get_paths() # paths at medlevel
    pcs = pcs[np.argmax([len(pc.vertices) for pc in pcs])] # longest path
    xcs,ycs = pcs.vertices[:,0],pcs.vertices[:,1]
    
    # First find the angle of the RGB
    th = -np.arctan(aright)/np.pi*180-5 # degrees, orientation of the RGB
    
    # rotate the contour
    Xcs, Ycs = rotate(xcs, ycs, th)
    # now rotate the CMD
    Xrot,Yrot = rotate(x, y, th)
    # tip of the rotated right selection line
    XAvrot = rotate(xAv0, yAv0, th)[0]
    ired = abs(Xcs-XAvrot)< 0.25# RGB part of the contour
    iredspan = 0.5
    while len(Ycs[ired])==0:
        iredspan*=1.5
        ired = abs(Xcs-XAvrot)< iredspan
    irgb = np.argmin(Ycs[ired])
    # fine search boxes
    xstep = 0.8
    irgb_box_xy = np.logical_and(abs(x-xcs[ired][irgb])<xstep, abs(y-ycs[ired][irgb])<1)
    
    # find the optimal rotation angle
    srgb = []
    # try: 
    #     srgb = np.array([rotate(x[irgb_box_xy],y[irgb_box_xy],th)[0].std() for th in np.linspace(80,90,21)])
    #     th = np.linspace(60,90,61)[np.nanargmin(srgb)]
    #     print('Angle of the RGB:',th)
    #     #rotate the contour
    #     Xcs, Ycs = rotate(xcs, ycs, th)
    #     #now rotate the CMD
    #     Xrot,Yrot = rotate(x, y, th)
    #     #find the tip of the rotated right selection line
    #     XAvrot = rotate(xAv0, yAv0, th)[0]
    #     ired = abs(Xcs-XAvrot)< 0.5# RGB part of the contour
    #     irgb = np.argmin(Ycs[ired]) # index of the center of the RGB from the GM contour
    #     irgb_box_xy = np.logical_and(abs(Xrot-Xcs[ired][irgb])<0.5, abs(Yrot-Ycs[ired][irgb])<0.5)
    # except:
    #     print('warning: failed to fit the angle of the rgb')
        
    # if len(srgb)>0 and np.nanmax(srgb)/np.nanmin(srgb)<1.2:    
    #     print('warning: failed to fit the angle of the rgb')
    #     th = -np.arctan(aright)/np.pi*180 # degrees, orientation of the RGB
    #     print('Angle of the RGB:',th)
        
        
    # fine search boxes
    irgb_box = np.logical_and(abs(xcs-xcs[ired][irgb])<0.5, abs(ycs-ycs[ired][irgb])<0.5)
   
    h = np.histogram(Xrot[irgb_box_xy],bins=50) # histogram of rotated V-I
    Xrgb = h[1][np.argmax(h[0])]
    if XAvrot>Xrgb+0.5:
        Xrgb = XAvrot
    print('Xrgb=',Xrgb)
    irgb = np.argmin(abs(Xcs[irgb_box]-Xrgb))
    # index of the center of the RGB from the GM density
    
    xrgb, yrgb = xcs[irgb_box][irgb], ycs[irgb_box][irgb] # true position
    Xrgb, Yrgb = Xcs[irgb_box][irgb], Ycs[irgb_box][irgb] # rotated position
    
    ycenrgb = np.linspace(Yrgb-2,Yrgb+1.5,8)
    rgb = np.logical_and(abs(Xrot-Xrgb-0.1)<xstep,abs(Yrot-Yrgb)<5)
    lrgb = np.array([len(Xrot[rgb][abs(Yrot[rgb]-Ycenrgb)<0.5]) for Ycenrgb in ycenrgb])
    if lrgb.mean()>200:
        ystep = 0.5
    else:
        ystep = 1
    
    lrgb = np.array([len(Xrot[rgb][abs(Yrot[rgb]-Ycenrgb)<ystep]) for Ycenrgb in ycenrgb])
    scenrgb = np.array([Xrot[rgb][abs(Yrot[rgb]-Ycenrgb)<ystep].std() for Ycenrgb in ycenrgb])
    xcenrgb = np.array([Xrot[rgb][abs(Yrot[rgb]-Ycenrgb)<ystep].mean() for Ycenrgb in ycenrgb])
    
    
    # remove Xrot-incomplete windows (at 1.5 sigma) and windows with few points
    scenrgb[xcenrgb-1.5*scenrgb<Xrgb-xstep-0.1] = 99.
    scenrgb[xcenrgb+1.5*scenrgb>Xrgb+xstep-0.1] = 99.
    scenrgb[ycenrgb>-0.8] = 99.

    if np.nanmin(scenrgb)<99:
        scenrgb[np.where(scenrgb<99)[0][lrgb[scenrgb<99]<min(100,max(lrgb[scenrgb<99]))]] = 99.
    
    Yrgb2 = ycenrgb[np.nanargmin(scenrgb)]
    print('Yrgb=',Yrgb2)
    rgb = np.logical_and(abs(Xrot-Xrgb-0.1)<xstep,abs(Yrot-Yrgb2)<ystep)
    
    EBV_std = Xrot[rgb].std()
    if EBV_std==0 or np.isnan(EBV_std) or np.nanmin(scenrgb)==99:
        # something wrong happened, back to default value
        EBV_std = 0.6
    else:
        box_rgb = (np.argmin(Xrot[rgb]-Yrot[rgb]),
               np.argmin(Xrot[rgb]+Yrot[rgb]),
               np.argmin(-Xrot[rgb]+Yrot[rgb]),
               np.argmin(-Xrot[rgb]-Yrot[rgb]),
               np.argmin(Xrot[rgb]-Yrot[rgb]))
        plt.plot(x[rgb][[box_rgb]][0],y[rgb][[box_rgb]][0],color='red',label='RGB')
    #plt.scatter(x[rgb],y[rgb],color='red',label='RGB ref',s=4)
    plt.plot(x[:1],y[:1],ls=':',lw=3,color='C0',label='GM contour')
    plt.quiver(x0,y0,xAv0-x0,yAv0-y0,color='blue',scale_units='xy', angles='xy', scale=1)
    plt.plot([xAv0,xf],[yAv0,-9],color='blue', lw=3)
    plt.plot([xs1,xsf],[ys1,-9],color='blue', lw=3, label="selection")
    
    
        
    # HR gap color cut
    cut_x = np.logical_and((aleft*(x-xs1)+ys1-y)*aleft/abs(aleft)>=0,(aright*(x-xAv0)+yAv0-y)*aright/abs(aright)<=0)
    # density and stellar mass cuts            
    cut_y = np.logical_and(y<np.interp(x,xmref,ymref),finterp(y/factor,x,grid=False)<medlevel)
    cand = np.logical_and(cut_x,cut_y)
    cAv = y[cand]-Rlambda*(Alambda/Av)*x[cand]
    xinter,yinter = (cAv-cleft)/(aleft-Rlambda*(Alambda/Av)), (Rlambda*(Alambda/Av)*cleft-aleft*cAv)/(Rlambda*(Alambda/Av)-aleft)
    d_cand = data[cand]
    d_cand['x'] = x[cand]
    d_cand['y'] = y[cand]
    EBV_std_corr = np.sqrt(max(EBV_std**2-0.12**2,1e-4))
    d_cand['Pextin'] = 2*(1-norm.cdf(x[cand]-xinter,0,EBV_std_corr))
    #Mstar = finterpMstar(y[cand],x[cand],grid=False)
    d_cand['Mstar'] = interpm(x[cand],y[cand])
    d_cand['logTeff'] = interpt(x[cand],y[cand])
    d_cand['logL'] = interpl(x[cand],y[cand])
    
    ncand = len(d_cand)
    print("number of candidates:",ncand)
    d_cand.write("../../../Downloads/Catalogs/HSCv3/%s_%s_candidates_%s_mist.fits"%(tget,source,col), overwrite=1)
    #plt.scatter(x[cand][::10],y[cand][::10],s=5,color='k',marker='x')
    plt.legend(loc=1, fontsize=17, framealpha=1, edgecolor='white')
    print("Required (MIST) magnitude", Mlim)
    print("Number of sources", len(x))
    #plt.gca().invert_yaxis()
    
    if len(np.where(cand)[0])>0 and zKDE[cand].max()>zKDEcut:
        print("Warning, there may be lost candidates")
        
    if fullplot:
        plt.xlim((np.percentile(x,1)-0.5,np.percentile(x,99)+1))
        plt.ylim((np.percentile(y,99)+0.5,min(y)-2))
    else:
        plt.xlim(np.percentile(x,2)-0.5,x[np.argsort(x)][-min(200,len(x)//20)]+1.5)
        plt.ylim(np.percentile(y,99.9)+0.5,y[np.argsort(y)][min(200,len(x)//20)]-2)
    #plt.savefig("../figures/cmd_%s.png"%(tget))
    plt.savefig("../figures/cmd_%s_%s_%s.png"%(tget,source,col))
    
    # Plot the RGB and the V-I scatter
    fig, axs = plt.subplots(2)
    axs[0].scatter(Xrot[rgb],Yrot[rgb],s=10,alpha=0.5)
    axs[0].scatter(Xrgb,Yrgb2,s=40,color='red')
    axs[0].invert_yaxis()
    #axs[0].xlabel('Rotated V-I')
    axs[0].set_ylabel('Rotated y')
    axs[1].hist(Xrot[rgb],bins='auto')
    axs[1].axvline(Xrot[rgb].mean(),color='k')
    axs[1].axvline(Xrot[rgb].mean()-EBV_std,ls='--',color='0.5')
    axs[1].axvline(Xrot[rgb].mean()+EBV_std,ls='--',color='0.5')
    axs[1].set_ylabel('N')
    axs[1].set_xlabel('Rotated x')
    for ax in fig.get_axes():
        ax.label_outer()
    fig.tight_layout()
    fig.savefig('../figures/%s_%s_rgb_scatter.png'%(tget,source))
    print("Scatter of V-I in the RGB:",EBV_std)
    print("Corrected scatter:",EBV_std_corr)
    
    # correct for intrinsic scatter ~ 0.12
    EBV_std_dic[gal['main_id']] = EBV_std_corr
    print("Contamination level: %.2f percent"%(np.nansum(d_cand['Pextin'])/ncand*100))
    
    if 1:
        plt.figure()
        plt.scatter(data['MatchRA'],data['MatchDec'])
        plt.gca().set_aspect(1/np.cos(np.mean(data['MatchDec'])*np.pi/180))
        plt.gca().invert_xaxis()
        plt.xlabel('RA (deg)')
        plt.ylabel('Dec (deg)')
        plt.tight_layout()
        plt.savefig("../figures/%s_%s_skyplot.png"%(tget,source))
    Igal+=1
    if Igal>=5:
        plt.close('all')
    if Igal%10==9:
        with open('EBV_std2.txt','w') as EBVfile:
            EBVfile.write(json.dumps(EBV_std_dic))
    
    Dcoo = SkyCoord(data['MatchRA'],data['MatchDec'],frame='icrs')
    dcoo = SkyCoord(d_cand['MatchRA'],d_cand['MatchDec'],frame='icrs')
    print(Dcoo.separation(SkyCoord(10.5333, 40.9169,unit='deg',frame='icrs')).min())
    print(dcoo.separation(SkyCoord(10.5333,40.9169 ,unit='deg',frame='icrs')).min())
    
    