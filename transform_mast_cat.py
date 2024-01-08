# -*- coding: utf-8 -*-
"""
Created on Mon Oct  2 14:21:58 2023

@author: hgtra
"""

from astropy.table import Table
from astropy.coordinates import SkyCoord
import astropy.units as u
import os 
import sys
import numpy as np
import matplotlib.pyplot as plt
import glob
from tqdm import tqdm

target = "M 31"
max_sep = "3mean" # arcsec
forceinst = ""

gals = Table.read(r'C:\Users\hgtra\Downloads\Catalogs\galaxy_mast_stats_addendum_1123.fits')
gals = Table.read(r'C:\Users\hgtra\Downloads\Catalogs\galaxies_hecate_20Mpc.fits')
for gal in tqdm(gals[:]):    
    tget = "".join(gal['main_id'].replace("NAME","").split())
    tget2 = " ".join(gal['main_id'].replace("NAME","").split())
    if target!='' and tget2!=target:
        continue
    
    flist = glob.glob("../../../Downloads/MAST_catalogs/%s/*"%tget)
    inst = ['ACS','WFC3','WFPC2'][np.argmax([np.sum([int('acs' in f) for f in flist]),
                                            np.sum([int('wfc3' in f) for f in flist]),
                                            np.sum([int('wfpc2' in f) for f in flist])])]
    print("Instrument",inst)
    #sys.exit()
    
    prefix = {"ACS":"A","WFC3":"W3","WFPC2":"W2"}
    
    fname = "../../../Downloads/MAST_catalogs/%s/catalog_%s.fits"%(tget,tget)
    if not(os.path.isfile(fname)):
        print("File not found")
        continue
        
    
    cat = Table.read(fname)
    
    filters = np.unique(cat['filter'])
    print("Available filters:",list(filters))
    
    if filters[0]=="f814w":
        col_v = filters[1].upper()
    else:
        col_v = filters[0].upper()
    
    catv = cat[cat['filter']==col_v.lower()]
    catr = cat[cat['filter']=='f814w']
    
    if len(catv)*len(catr)==0:
        print("empty red or green cat")
        continue
    
    coordv = SkyCoord(ra=np.asarray(catv['RA'])*u.degree, dec=np.asarray(catv['DEC'])*u.degree)
    coordr = SkyCoord(ra=np.asarray(catr['RA'])*u.degree, dec=np.asarray(catr['DEC'])*u.degree)
    idx, d2d, d3d = coordv.match_to_catalog_sky(coordr)
    if max_sep=="3mean":
        max_sep = min(2,np.mean(d2d.arcsec)*3) # arcsec
        print("Maximum separation:",max_sep)
    
    sep_constraint = d2d.arcsec <= max_sep
    coordv, coordr = coordv[sep_constraint], coordr[idx][sep_constraint]
    matchv, matchr = catv[sep_constraint], catr[idx][sep_constraint]
    if len(matchv)==0:
        print("empty match")
        continue
    
    plt.figure()
    plt.hist(d2d.arcsec+1.2e-3, bins=np.geomspace(1e-3, 300,30))
    plt.xlabel("Separation (arcsec)")
    plt.ylabel("Number of sources")
    plt.axvline(max_sep,color='k',label='threshold')
    plt.loglog()
    plt.legend(loc=2)
    plt.savefig("../figures/%s_mast_cat_separation.png"%tget2)
    
    print("Number of matches / number of V sources:", len(matchv), "/", len(catv))
    
    filtcols = ["A_F435W","W3_F438W","W2_F439W","A_F475W","W3_F475W","A_F555W","W3_F555W","W2_F555W","A_F606W","W3_F606W","W2_F606W"]
    
    
    if not(forceinst):
        ins = prefix[inst]
    else:
        ins = forceinst
    matchv['MagAp2'].name = '%s_%s'%(ins,col_v)
    matchv['CI'].name = '%s_CI'%col_v
    if not('MagErr2') in matchv.colnames:
        matchv['MagErrAp2'].name = '%s_%s_err'%(ins,col_v) 
    elif not('MagErrAp2') in matchv.colnames:
        matchv['MagErr2'].name = '%s_%s_err'%(ins,col_v) 
    else:
        matchv['%s_%s_err'%(ins,col_v)] = np.minimum(np.nan_to_num(np.asarray(matchv['MagErrAp2']),nan=99),np.nan_to_num(np.asarray(matchv['MagErr2']),nan=99))
        
    
    matchv['%s_F814W'%ins] = matchr['MagAp2']
    
    if not('MagErr2') in matchr.colnames:
        matchv['%s_F814W_err'%ins] = matchr['MagErrAp2']
    elif not('MagErrAp2') in matchr.colnames:
        matchv['%s_F814W_err'%ins] = matchr['MagErr2']
    else:
        matchv['%s_F814W_err'%ins] = np.minimum(np.nan_to_num(np.asarray(matchr['MagErrAp2']),nan=99),np.nan_to_num(np.asarray(matchr['MagErr2']),nan=99))
    
    matchv['F814W_CI'] = matchr['CI']
    
    matchv['RA'].name = 'MatchRA'
    matchv['DEC'].name = 'MatchDec'
    
    
    for c in filtcols:
        if not(c in matchv.colnames):
            matchv[c] = np.nan
            
    
    igoodmag = np.maximum(matchv['%s_F814W_err'%ins],matchv['%s_%s_err'%(ins,col_v)])<0.1
    igoodx = np.maximum(np.abs(matchv["%s_%s"%(ins,col_v)]),np.abs(matchv["%s_F814W"%ins]))<99
    igoodci = matchv['%s_CI'%col_v]>0.9-matchv['%s_%s_err'%(ins,col_v)]/0.5
    igoodci2 = matchv['F814W_CI']>0.9-matchv['%s_F814W_err'%ins]/0.5
    
    igood = np.logical_and(np.logical_and(igoodmag,igoodx),
                           np.logical_and(igoodci,igoodci2))
    #igood = np.logical_and(igoodmag,igoodx)
    #igood = igoodx
    
    coordv = coordv[igood]
    matchv = matchv[igood]
    
    #matchv, coordv = matchv[::-1], coordv[::-1]
    
    k=0
    sep = None
    while len(idx[sep_constraint])>0 and k<100:
        
        print(k, len(coordv))
        idx, d2d, d3d = coordv.match_to_catalog_sky(coordv, nthneighbor=2)
        if sep is None:
            sep = max(min(0.3,np.mean(d2d.arcsec)),0.1)
            plt.figure()
            plt.hist(np.log10(d2d.arcsec+1e-3),bins='auto')
            
            plt.axvline(np.log10(sep),label='threshold',color='k')
            plt.legend()
            plt.gca().set_yscale('log')
            plt.xlabel('Log. internal separation (arcsec)')
            plt.ylabel('Number of sources')
            print("Mean of internal separation: %f arcsec"%sep)
            print("Median internal separation: %f arcsec"%np.median(d2d.arcsec))
            plt.figure()
            plt.scatter(np.log10(d2d.arcsec+1e-3),np.log10(np.abs(matchv['%s_F814W'%ins]-matchv[idx]['%s_F814W'%ins])+np.abs(matchv['%s_%s'%(ins,col_v)]-matchv[idx]['%s_%s'%(ins,col_v)])+1e-3),s=1,alpha=0.1)
            plt.xlabel('Log. internal separation (arcsec)')
            plt.ylabel('Log. internal delta-mag (mag)')
            plt.axvline(np.log10(sep),color='k')
            plt.figure()
            ieqmag = np.abs(matchv['%s_F814W'%ins]-matchv[idx]['%s_F814W'%ins])+np.abs(matchv['%s_%s'%(ins,col_v)]-matchv[idx]['%s_%s'%(ins,col_v)])<1e-5
            plt.hist(d2d[ieqmag].arcsec,bins='auto')
            plt.axvline(sep,label='threshold',color='k')
            plt.legend()
            plt.xlabel('Internal separation of mag-matches (arcsec)')
            plt.gca().set_yscale('log')
            plt.ylabel('Number of sources')
        sep_constraint = d2d.arcsec <= sep
        coordv = coordv[~np.isin(np.arange(len(idx)),idx[sep_constraint])]
        matchv = matchv[~np.isin(np.arange(len(idx)),idx[sep_constraint])]
        if len(matchv)==0:
            break
        k+=1
    if target=="":
        plt.close('all')            
    if len(matchv)==0:
        print("empty match")
        continue
    
    matchv.write('../../../Downloads/MAST_catalogs/%s_data.fits'%(tget2), overwrite=True)
    
        
        