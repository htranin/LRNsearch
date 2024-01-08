# -*- coding: utf-8 -*-
"""
Created on Mon Jun 26 16:36:06 2023

@author: Hugo Tranin


Take as input the initial list of galaxies for dolphot, the aim being to determine those for which it is really necessary to use dolphot
Search for MAST data for each galaxy centered in its R1 circle
Select only
 Hubble data with t_exp>300s, both F814W and green filter, calibration level 2 or 3
Try to select the minimum number of observations giving coverage of the galaxy as complete as possible.
  To do this, create and compare MOCs for each observation
Assign each galaxy its number of sexphot/daophot catalogs and its number of FLC images.
Note that some of these galaxies are present in HSCv3, such as M110.

"""

from astroquery.mast import Observations
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from mocpy import MOC
import os
from astropy.table import Table, vstack
import requests
from tqdm import tqdm
import matplotlib.pyplot as plt
from astropy.time import Time
import casjobs_query

fewflc = 0
verbose = 1
skip_selection = 0
tget = "M 31"
nocat_only = 0
getcat = 1 # shortcut for optimal catalog retrieval, from V+I observations with same obs_id2
dl_products = 0 # not only catalogs
obsids_orig = []#'hst_10585_02_acs_wfc_f555w', 'hst_10585_02_acs_wfc_f814w', 'hst_10585_01_acs_wfc_f555w', 'hst_10585_01_acs_wfc_f814w', 'hst_10585_01_acs_wfc_f555w_j9cd01', 'hst_10585_01_acs_wfc_f814w_j9cd01']
simultaneous = 1 # keep only simultaneous observations in V+I
selected = 'F475W' # force filter choice
casjobs = 0 # run CasJobs before trying MAST -> favour HSCv3 data

#gals = Table.read(r'C:\Users\hgtra\Downloads\Catalogs\sample_galaxies_for_dolphot.fits')
gals = Table.read(r'C:\Users\hgtra\Downloads\Catalogs\galaxy_compilation_Simbad_unique_0823.fits')
gals = Table.read(r'C:\Users\hgtra\Downloads\Catalogs\galaxies_hecate_20Mpc_new.fits')
#gals = gals[~gals['in_compilation']]
#gals = Table.read(r'C:\Users\hgtra\Downloads\Catalogs\galaxy_mast_stats_addendum_1123.fits')
#gals = Table.read(r'C:\Users\hgtra\Downloads\Catalogs\galaxy_addendum_to_fix_retrieve_mast_error.fits')
#gals = Table.read(r'C:\Users\hgtra\Downloads\Catalogs\galaxy_hst_0723.fits')
gal = gals[np.array([tget.replace(' ','')==name.replace(' ','') for name in gals['main_id']])]


#gals = gals[np.argsort(gals['R1'])[::-1]][2:]
# NGC4236, NGC3109, M110...

os.chdir(r'C:\Users\hgtra\Downloads\MAST_catalogs')

selected_filters = ['F336W', 'F435W', 'F438W', 'F439W', 'F475W', 'F555W', 'F606W', 'F814W']

if skip_selection:
    filewrite='galaxies_mast_data_all8.csv'
    if not(os.path.isfile(filewrite)):
        f = open(filewrite,'w')
        print('main_id,nobs,cov_v,cov_i,cov_vi,',
              ','.join(['cov_%s'%filt for filt in selected_filters[:-1]]),
              ','+','.join(['tmin_%s'%filt for filt in selected_filters]),
              ','+','.join(['tmax_%s'%filt for filt in selected_filters]),file=f)
    else:
        f = open(filewrite,'a')
else:
    filewrite='galaxies_mast_data_3.csv'
    if not(os.path.isfile(filewrite)):
        f = open(filewrite,'w')
        print('main_id,ncat,nobs,nflc,cov_sel,cov_hst,col_v',file=f)
    else:
        f = open(filewrite,'a')


g = np.genfromtxt('galaxies_mast_data.csv',names=True,delimiter=',',dtype=None,encoding='utf-8')

if nocat_only:
    gals_no_cat = Table.read("../Catalogs/galaxy_mast_stats.fits")
    gals_no_cat = gals_no_cat[gals_no_cat['ncat']==0]
    
    gals['nocat'] = [int(n in gals_no_cat['main_id']) for n in gals['main_id']]
    gals = gals[gals['nocat']==1]





for gal in tqdm(gals[:]):
    obsids = obsids_orig
    if gal['R1']>6000:
        print("galaxy too large")
        continue
    if tget!='' and gal['main_id'].replace(' ','')!=tget.replace(' ',''):
        continue
    
    if tget=="M 31":
        print("Complement M31")
        # only get MAST complement of HSCv3
        gal['ra'] = 10.1
        gal['dec'] = 40.7
        gal['R1'] = 28*60
        gal['R2'] = 14*60
        
    ra, dec = gal['ra'], gal['dec']
    ddec = gal['R1']/3600 # deg
    dra = gal['R1']/np.cos(np.radians(gal['dec']))/3600
    
    ramin, ramax = ra-dra, ra+dra
    decmin, decmax = dec-ddec, dec+ddec
    
    if gal['R2']>=gal['R1']:
        gal['R2'] = 0.999*gal['R1']
    if gal['PA']<0:
        gal['PA']+=180
        
    
    if verbose:
        print(gal['main_id','R1','R2','PA'])
    target = gal['main_id']
    
    
    if casjobs:
        print("running CasJobs query")
        os.chdir("NGC300") # trick to keep relative path working
        casjobs_query.run_query(target.replace(" ","").replace("NAME",""))
        galsfailed = np.genfromtxt("GALAXIES_TO_DOWNLOAD_ON_CASJOBS.txt",delimiter=",",encoding="utf-8",dtype=None)
        os.chdir("..")
        hscfile = "../Catalogs/HSCv3/%s_data.fits"%(" ".join(target.replace("NAME","").split()))
        if os.path.isfile(hscfile):
            data = Table.read(hscfile)
            print("retrieved HSCv3 data of length",len(data))
            if len(data)>400:
                continue
        
        if len(galsfailed)>0 and galsfailed[-1]==target:
            continue
        print("CasJobs query done, now running MAST query...")
        
    dirtarget = target.replace('NAME','').replace(' ','')
    if g[g['main_id']==target]['nflc']>100 and fewflc:
        print('stopped at nflc step')
        continue
    
    f.close()
    f = open(filewrite,'a')
    r1 = float(gal['R1']/3600)+0.03 #degrees
    # 0.03 accounts for the extent of HST observations (<3.6arcmin squares)
    
    if 0:
        print('skipping query because %d observations already loaded'%len(obs_table_init))
    else:
        try:
            obs_table_init = Observations.query_object(target, radius="%f deg"%r1)
        except Exception as e:
            if skip_selection:
                print(dirtarget+(3+3*len(selected_filters))*',',file=f)
            else:
                print(dirtarget+',,,,,,',file=f)
            print(e)
            continue
        
        
    if obsids:
        obs_table = obs_table_init[np.isin(obs_table_init['obs_id'],obsids)]
        print("Retrieved",len(obs_table),"observations from obsid list")
        continue
    else:
        obs_table = obs_table_init[:] 
        print(0,len(obs_table))
        
        
        obs_table = obs_table[np.logical_or(obs_table['obs_collection']=='HST',
                                            obs_table['obs_collection']=='HLA')]
        print(1,len(obs_table))
        obs_table = obs_table[obs_table['t_exptime']>2e2]
        print(2,len(obs_table))
        obs_table_cat = obs_table[obs_table['calib_level']==2]
        
        # Save for later use
        
        obs_table = obs_table[obs_table['calib_level']>=2]
        print(3,len(obs_table))
        
        obs_table = obs_table[np.isin(obs_table['filters'],selected_filters)]
        print(4,len(obs_table))
        obs_table = obs_table[~np.logical_and(obs_table['filters']=="F336W",np.isin(obs_table['instrument_name'],["WFPC2","WFPC2/WFC","WFPC2/PC"]))]
        print(5,len(obs_table))
        
        obs_table = obs_table[np.argsort(obs_table['t_exptime'])[::-1]]
        
        
    if len(obs_table)==0:
        if verbose:
            print('no HST observations, skipping')
        if skip_selection:
            print(target+(3+3*len(selected_filters))*',',file=f)
        else:
            print(target+',,,,,,',file=f)
        continue
    
    if not('F814W' in obs_table['filters'] and np.isin(obs_table['filters'],selected_filters[:-1]).any()):
        if verbose:
            print('only red or only green observations, skipping')
        if skip_selection:
            print(target+(3+3*len(selected_filters))*',',file=f)
        else:
            print(target+',,,,,,',file=f)
        continue
    
    if verbose:
        print("Working on %d obsids"%len(obs_table))
    
    # add information on observations coverage
    
    mocs = []
    covs = []
    moctot = None
    moc_green = None
    moc_filters = {filt:None for filt in selected_filters[:-1]}
    moc_red = None
    moc_gal = MOC.from_elliptical_cone(lon=gal['ra']*u.deg,lat=gal['dec']*u.deg,a=gal['R1']*u.arcsec,b=gal['R2']*u.arcsec,pa=gal['PA']*u.deg, max_depth=16)
    for o in obs_table:
        m = None
        for s_region in o['s_region'].replace('J2000','').split('POLYGON'):
            if len(s_region)==0: #empty string
                continue
            try:
                vertices = np.array(s_region.split()).astype(float).reshape((-1,2))
            except:
                #if skip_selection:
                #    print(target+',,,,,,,,,',file=f)
                #else:
                #    print(target+',,,,,,',file=f)
                continue
            
            if m is None:
                m = MOC.from_polygon(vertices[:,0]*u.deg,vertices[:,1]*u.deg,max_depth=16)
            else:
                m = m.union(MOC.from_polygon(vertices[:,0]*u.deg,vertices[:,1]*u.deg,max_depth=16))
        
        if moctot is None:
            moctot = m
        else:
            
            moctot = moctot.union(m)
        
        if o['filters']=="F814W":
            if moc_red is None:
                moc_red = m
            else:
                moc_red = moc_red.union(m)
        else:
            if moc_green is None:
                moc_green = m
            else:
                moc_green = moc_green.union(m)
            if moc_filters[o['filters']] is None:
                moc_filters[o['filters']] = m
            else:
                moc_filters[o['filters']] = moc_filters[o['filters']].union(m)
            
        mocs.append(m)
        covs.append(float("%.2e"%(m.sky_fraction*41253*3600))) #squared arcmin
    
    if moc_green is not None:
        moc_greeng = moc_gal.intersection(moc_green)
    else:
        print("empty green MOC (1), skipping")
        continue
    if moc_red is not None:
        moc_redg = moc_gal.intersection(moc_red)
    else:
        print("empty red MOC, skipping")
        continue
    
    for filt in selected_filters[:-1]:
        if moc_filters[filt] is not None:
            moc_filters[filt] = moc_redg.intersection(moc_filters[filt])
            
    
    # Selection of the green color
    cov_filters = [0 if moc_filters[filt] is None else moc_filters[filt].sky_fraction for filt in selected_filters[:-1]]
    col_v = selected_filters[np.argmax(np.array(cov_filters)*np.array([float(f[1:-1]) for f in selected_filters[:-1]]))]
    if selected:
        col_v = selected
        
        
    # multiply by lambda to favour redder colours in the end...
    if verbose:
        print("selected green color:", col_v)
    
    tmins = ["%d"%np.min(obs_table['t_exptime'][obs_table['filters']==filt]) if filt in obs_table['filters'] else '' for filt in selected_filters]
    tmaxs = ["%d"%np.max(obs_table['t_exptime'][obs_table['filters']==filt]) if filt in obs_table['filters'] else '' for filt in selected_filters]
    
    
    
    if skip_selection:
        if moc_red is None:
            print(target+',%d,%f,%f,%f,%s,%s,%s'%(len(obs_table),
                                          moc_green.sky_fraction*41253*3600,
                                          0,
                                          0,
                                          ','.join(['%f'%(cov*41253*3600) for cov in cov_filters]),
                                          ','.join(tmins),
                                          ','.join(tmaxs)),file=f)
        elif moc_green is None:
            print(target+',%d,%f,%f,%f,%s,%s,%s'%(len(obs_table),
                                          0,
                                          moc_red.sky_fraction*41253*3600,
                                          0,
                                          ','.join(['%f'%(cov*41253*3600) for cov in cov_filters]),
                                          ','.join(tmins),
                                          ','.join(tmaxs)),file=f)
        else:
            moc_redgreen = moc_red.intersection(moc_green)
            print(target+',%d,%f,%f,%f,%s,%s,%s'%(len(obs_table),
                                          moc_green.sky_fraction*41253*3600,
                                          moc_red.sky_fraction*41253*3600,
                                          moc_redgreen.sky_fraction*41253*3600,
                                          ','.join(['%f'%(cov*41253*3600) for cov in cov_filters]),
                                          ','.join(tmins),
                                          ','.join(tmaxs)),file=f)
        continue
    
    
    ### Catalog enhancement
    moc_green = moc_filters[col_v]
    if moc_green is None:
        print("empty green MOC (2), skipping")
        continue
    
    obs_table['coverage'] = covs
    moc_redgreen = moc_red.intersection(moc_green)
    i_redgreen = np.logical_or(obs_table['filters']=='F814W',
                               obs_table['filters']==col_v)
    mocs = np.array(mocs)[i_redgreen]
    obs_table = obs_table[i_redgreen]
    if verbose:
        print("Total coverage:",moctot.sky_fraction*41253*3600,"arcmin^2")
        print("Green coverage:",moc_green.sky_fraction*41253*3600,"arcmin^2")
        print("Red coverage:",moc_red.sky_fraction*41253*3600,"arcmin^2")
        print("Red+green coverage:",moc_redgreen.sky_fraction*41253*3600,"arcmin^2")
        
    obs_table['FoM'] = obs_table['t_exptime']#obs_table['coverage']*obs_table['t_exptime']**0.5
    obs_table['FoM'] = obs_table['t_exptime']/10**np.array(['skycell' in o['obs_id'] for o in obs_table])
    
    # Proxy of the number of objects
    mocs = np.array(mocs)[np.argsort(obs_table['FoM'])[::-1]]
    obs_table = obs_table[np.argsort(obs_table['FoM'])[::-1]]
    mocs = mocs[[oid[-3:]!='_pc' for oid in obs_table['obs_id']]]
    obs_table = obs_table[[oid[-3:]!='_pc' for oid in obs_table['obs_id']]]
    #obs_table['obs_id2'] = [o[:2*len(o)//3] for o in obs_table['obs_id']]
    obs_table['obs_id2'] = [o.replace('hst_','')[:5] for o in obs_table['obs_id']]
    if simultaneous:
        mocs = mocs[np.logical_or([o['filters']=='F814W' and o['obs_id2'] in obs_table[obs_table['filters']==col_v]['obs_id2'] for o in obs_table],
                                  [o['filters']==col_v and o['obs_id2'] in obs_table[obs_table['filters']=='F814W']['obs_id2'] for o in obs_table])]
        obs_table = obs_table[np.logical_or([o['filters']=='F814W' and o['obs_id2'] in obs_table[obs_table['filters']==col_v]['obs_id2'] for o in obs_table],
                                            [o['filters']==col_v and o['obs_id2'] in obs_table[obs_table['filters']=='F814W']['obs_id2'] for o in obs_table])]
    
    
    if len(obs_table)==0:
        print("empty obs_table, skipping")
        continue
    
    obs_table['decyear'] = Time(obs_table['t_max'],format='mjd').to_value('decimalyear', subfmt='str')
    #obs_table = obs_table[np.logical_or(
    #                      np.logical_or(["f814w" in o['obs_id'] and o['obs_id'].replace('f814w',col_v.lower()) in obs_table['obs_id'] for o in obs_table],
    #                                    [col_v.lower() in o['obs_id'] and o['obs_id'].replace(col_v.lower(),'f814w') in obs_table['obs_id'] for o in obs_table]),
    #                      np.logical_or([not("_" in o['obs_id']) and str(int(o['obsid'])-1) in obs_table['obsid'] for o in obs_table],
    #                                    [not("_" in o['obs_id']) and str(int(o['obsid'])+1) in obs_table['obsid'] for o in obs_table]))]
    
    
    if verbose:
        print("Starting selection")
        obs_table_inter = obs_table[:]
        print("First ten observations:",obs_table['obs_id','FoM','filters'][:10])
    
    
    
    ### Select observations
    
    
    
    green_sel = np.zeros(len(obs_table))
    green_sel_moc = None
    red_sel = np.zeros(len(obs_table))
    red_sel_moc = None
    keep_1duplicate = 1
            
    # First the red ones
    
    for i, o in zip(range(len(obs_table)),obs_table):
        if o['filters']!='F814W':
            continue
        
        if red_sel_moc is None:
            red_sel_moc = mocs[i]
            red_sel[i] = 1
            print("selected", obs_table[i]['obs_id'])
            print("updated coverage", red_sel_moc.sky_fraction*41253*60**2)
            
        else:
            
            if red_sel_moc.sky_fraction > 0.95 *moc_red.sky_fraction:
                # selection is enough
                if verbose:
                    print('red selection is enough on iteration',i)
                break
            if mocs[i].intersection(moc_green).sky_fraction<mocs[i].sky_fraction*min(0.05,1/len(mocs)):
                # little intersection with green observations
                if verbose:
                    print(obs_table[i]['obs_id'],'little intersection with green observations')
                continue
            if  i!=0 and mocs[i].intersection(red_sel_moc).sky_fraction>0.98*mocs[i].sky_fraction:
                if keep_1duplicate == 1:
                    keep_1duplicate = 0
                else:
                    # nearly duplicate field
                    if verbose:
                        print(obs_table[i]['obs_id'],'nearly duplicate red')
                    keep_1duplicate = 1
                    continue
            if mocs[i].sky_fraction<red_sel_moc.sky_fraction*min(0.05,1/len(mocs)):
                # too small field
                if verbose:
                    print(obs_table[i]['obs_id'],'red field too small',i, mocs[i].sky_fraction*41253*3600)
                continue
            products = Observations.get_product_list(obs_table[i]['obsid'])
            if getcat and not(np.any(['cat' in p for p in products['productFilename']])):
                print(obs_table[i]['obs_id'],"without associated catalog")
                continue
            print("selected", obs_table[i]['obs_id'])
            red_sel_moc = red_sel_moc.union(mocs[i])
            print("updated coverage", red_sel_moc.sky_fraction*41253*60**2)
            red_sel[i] = 1
    
    # Then the green ones
    if getcat:
        green_sel = np.array([int(o['filters']==col_v and o['obs_id2'] in obs_table[red_sel==1]['obs_id2']) for o in obs_table])
    if getcat and simultaneous:
        print(sum(red_sel),'selected red observations')
        print(sum(green_sel),'selected green observations')
        for i, o in zip(range(len(obs_table)),obs_table):
            if green_sel[i]==0:
                continue
            if green_sel_moc is None:
                green_sel_moc = mocs[i]
            else:
                green_sel_moc = green_sel_moc.union(mocs[i])
                
    else:
        for i, o in zip(range(len(obs_table)),obs_table):
            if o['filters']=='F814W':
                continue
            
            if green_sel_moc is None:
                products = Observations.get_product_list(obs_table[i]['obsid'])
                if getcat and not(np.any(['cat' in p for p in products['productFilename']])):
                    print(obs_table[i]['obs_id'],"without associated catalog")
                    continue
                green_sel_moc = mocs[i]
                print("selected", obs_table[i]['obs_id'])
                green_sel[i] = 1  
                print("updated coverage", green_sel_moc.sky_fraction*41253*60**2)
            else:
                redgreen_sel_moc = green_sel_moc.intersection(red_sel_moc) 
                if 0 and (green_sel_moc.sky_fraction > 0.95 *moc_green.sky_fraction or redgreen_sel_moc.sky_fraction > 0.8 *moc_redgreen.sky_fraction):
                    # selection is enough
                    if verbose:
                        print('green selection is enough on iteration',i)
                    break
                if mocs[i].intersection(moc_red).sky_fraction<mocs[i].sky_fraction*min(0.05,1/len(mocs)):
                    # little intersection with red observations
                    if verbose:
                        print(obs_table[i]['obs_id'],'little intersection with red observations')
                    continue
                if i!=0 and mocs[i].intersection(green_sel_moc).sky_fraction>0.98*mocs[i].sky_fraction:
                    # nearly duplicate field
                    if verbose:
                        print(obs_table[i]['obs_id'],'nearly duplicate green')
                    continue
                if mocs[i].sky_fraction<green_sel_moc.sky_fraction*min(0.05,1/len(mocs)):
                    # too small field
                    if verbose:
                        print(obs_table[i]['obs_id'],'green field too small',i, mocs[i].sky_fraction*41253*3600)
                    continue
                products = Observations.get_product_list(obs_table[i]['obsid'])
                if getcat and not(np.any(['cat' in p for p in products['productFilename']])):
                    print(obs_table[i]['obs_id'],"without associated catalog")
                    continue
                print("selected", obs_table[i]['obs_id'])
                green_sel_moc = green_sel_moc.union(mocs[i])
                print("updated coverage", green_sel_moc.sky_fraction*41253*60**2)
                
                
                green_sel[i] = 1    
    
    if green_sel_moc is None:
        print("empty green MOC (3), skipping")
        continue
    
    redgreen_sel_moc = green_sel_moc.intersection(red_sel_moc)    
    obs_table['selected'] = green_sel+red_sel
    
    obs_table = obs_table[obs_table['selected']==1]
    ratest = ramin+(ramax-ramin)*np.random.random(10000)
    dectest = decmin+(decmax-decmin)*np.random.random(10000)
    imocred = red_sel_moc.contains(ratest*u.deg,dectest*u.deg)
    imocgreen = green_sel_moc.contains(ratest*u.deg,dectest*u.deg)
    plt.figure()
    plt.scatter(ratest[imocred],dectest[imocred],s=2,color='r',label='red',marker='o')
    plt.scatter(ratest[imocgreen],dectest[imocgreen],s=1,color='g',label='green',marker='_')
    plt.scatter([ 10.5333],[40.9169 ],s=50,label='M31-LRN2015')
    plt.legend(fontsize=10)
    plt.xlabel('RA')
    plt.ylabel('Dec')
    plt.gca().invert_xaxis()
    plt.savefig("../../OneDrive/Documents/figures/%s_MOC_hst_observations.pdf"%dirtarget)   
    
    
    if verbose:
        print('Selected red+green coverage:', redgreen_sel_moc.sky_fraction*3600*41253)    
        print(len(obs_table),'\n',obs_table['target_name','obsid','filters','t_exptime','decyear'][:20])
        print("M31-POS26" in obs_table['target_name'],
              "M31-POS41" in obs_table['target_name'])
    
    ### Print product list
    
    obsids = obs_table['obsid']
    products = Observations.get_product_list(obsids)
    flcproducts = products[np.array(['_flc' in p['productFilename'] for p in products])]
    if verbose:
        print(flcproducts['productFilename'])
    
    # Any DAOPHOT and Sextractor catalog?
    
    starcat = products[['cat' in p for p in products['productFilename']]]
    starcat = starcat[[not('total' in n or 'segment' in n or 'multiwave' in n) for n in starcat['productFilename']]]
    starcat = starcat[starcat['dataRights']=="PUBLIC"]
    if len(starcat)>0 and getcat: # not(os.path.isfile(dirtarget+'/catalog_'+dirtarget+'.fits'))
        if verbose:
            print("Found catalogs:\n",starcat['productFilename'])
            
            
        if not(os.path.exists(dirtarget)):
            os.makedirs(dirtarget)
        
        os.chdir(dirtarget)
        if os.path.isfile('catalog_'+dirtarget+'.fits'):
            os.replace('catalog_'+dirtarget+'.fits','catalog_'+dirtarget+'_old.fits')
            
        
        if not(os.path.isfile('catalog_'+dirtarget+'.fits')):
            for uri in starcat['dataURI']:
                fname = uri.split('/')[-1].split('filename=')[-1]
                if not(os.path.isfile(fname)):
                    try:
                        Observations.download_file(uri)
                    except:
                        r=requests.get('https://hla.stsci.edu/cgi-bin/'+uri.split('/')[-1])
                        with open(fname,'wb') as fl:
                            fl.write(r.content)
            Cat = None
            sexnames=' X-Center    Y-Center         RA         DEC        ID     TotMag     MagErr       CI    Flags     MagAp1     MagErr1    FluxAp1    FluxErr1     MagAp2     MagErr2    FluxAp2    FluxErr2        Bck    FlagsSE        pix        deg     MagIso     MagErrIso    FluxIso    FluxErrIso  MagIsoCor     MagErrIsoCor FluxIsoCor    FluxErrIsocor   Radius    MagAuto     MagErrAuto   FluxAuto    FluxErrAuto    MagBest     MagErrBest   FluxBest    FluxErrBest     countsTh        MagTh     countsPk        MagPk          nl0          nl1          nl2          nl3          nl4          nl5          nl6          nl7       Xmin       Ymin       Xmax       Ymax         X2      ErrX2         Y2      ErrY2         XY      ErrXY        CXX     ErrCXX        CYY     ErrCYY        CXY     ErrCXY        A     ErrA        B     ErrB    ThetaPix  ErrThetaPix    FWHM             W_X             W_Y            W_X2         W_ErrX2            W_Y2         W_ErrY2            W_XY         W_ErrXY           W_CXX        W_ErrCXX           W_CYY        W_ErrCYY           W_CXY        W_ErrCXY             W_A          W_ErrA             W_B          W_ErrB    W_Theta W_ErrTheta     W_FWHM    Elong    Ellip    Theta ErrTheta          RA1950         DEC1950    Theta1950 ErrTheta1950          RA_sky         DEC_sky    Theta_sky ErrTheta_sky Star/Gal       c108       c109       c110       c111       c112       c113       c114       c115       c116       c117       c118       c119       c120       c121       c122       c123       c124       c125       c126       c127       c128       c129       c130       c131       c132      X-Det      Y-Det       Back       Cnts       XCTE       YCTE Chip    Mag-CTE'.split()
            sexnames2=' X-Center    Y-Center         RA         DEC        ID     TotMag     MagErr       CI    Flags     MagAp1     MagErr1    FluxAp1    FluxErr1     MagAp2     MagErr2    FluxAp2    FluxErr2        Bck    FlagsSE        pix        deg     MagIso     MagErrIso    FluxIso    FluxErrIso  MagIsoCor     MagErrIsoCor FluxIsoCor    FluxErrIsocor   Radius    MagAuto     MagErrAuto   FluxAuto    FluxErrAuto    MagBest     MagErrBest   FluxBest    FluxErrBest     countsTh        MagTh     countsPk        MagPk          nl0          nl1          nl2          nl3          nl4          nl5          nl6          nl7       Xmin       Ymin       Xmax       Ymax         X2      ErrX2         Y2      ErrY2         XY      ErrXY        CXX     ErrCXX        CYY     ErrCYY        CXY     ErrCXY        A     ErrA        B     ErrB    ThetaPix  ErrThetaPix    FWHM             W_X             W_Y            W_X2         W_ErrX2            W_Y2         W_ErrY2            W_XY         W_ErrXY           W_CXX        W_ErrCXX           W_CYY        W_ErrCYY           W_CXY        W_ErrCXY             W_A          W_ErrA             W_B          W_ErrB    W_Theta W_ErrTheta     W_FWHM    Elong    Ellip    Theta ErrTheta          RA1950         DEC1950    Theta1950 ErrTheta1950          RA_sky         DEC_sky    Theta_sky ErrTheta_sky Star/Gal       c108       c109       c110       c111       c112       c113       c114       c115       c116       c117       c118       c119       c120       c121       c122       c123       c124       c125       c126       c127       c128       c129       c130       c131       c132'.split()
            daonames='X-Center   Y-Center         RA          DEC        ID   MagAp1   MagErr1   MagAp2   MagErr2     MSky    Stdev       Flux   TotMag   MagErr       CI    Flags      X-Det      Y-Det       Back       Cnts     XCTE     YCTE Chip  Mag-CTE'.split()
            daonames2='X-Center   Y-Center         RA          DEC        ID   MagAp1   MagErr1   MagAp2   MagErr2     MSky    Stdev       Flux   TotMag   MagErr       CI    Flags'.split()
            catnames = ['X-Center','Y-Center','RA','DEC','ID','MagAp1','MagErr1','MagAp2','MagErr2','CI','Flags','ref']
            
            for fl in starcat['productFilename']:
                if not(os.path.isfile(fl)):
                    continue
                if fl[-4:]=='ecsv':
                    cat = Table.read(fl, format='ascii.ecsv')
                    if len(cat)==0:
                        if verbose:
                            print('empty table, skipping')
                        continue
                    cat['ref'] = 'pnt-cat'
                else:
                    
                    cat = Table.read(fl, format='ascii.basic')
                    if len(cat)==0:
                        continue
                    if len(cat.colnames)==140:
                        cat.rename_columns(cat.colnames,sexnames)
                        cat['ref'] = 'sexphot'
                    elif len(cat.colnames)==132:
                        cat.rename_columns(cat.colnames,sexnames2)
                        cat['ref'] = 'sexphot'
                    elif len(cat.colnames)==24:
                        cat.rename_columns(cat.colnames,daonames)
                        cat['ref'] = 'daophot'
                    elif len(cat.colnames)==16:
                        cat.rename_columns(cat.colnames,daonames2)
                        cat['ref'] = 'daophot'
                    else:
                        if verbose:
                            print('failed to guess column names from column number, skip')
                        continue
                    cat.remove_columns([c for c in cat.colnames if not(c in catnames)])
                
                cat = cat[np.logical_and(abs(cat['RA']-(ramin+ramax)/2)<(ramax-ramin)/2,
                                         abs(cat['DEC']-(decmin+decmax)/2)<(decmax-decmin)/2)]
                if len(cat)==0:
                    if verbose:
                        print('no match with galaxy, skipping')
                    continue
                
                cat.meta = None
                filt = 'f'+fl.split('_f')[1].split('_')[0]
                if len(filt)!=5:
                    filt = 'f'+fl.split('_f')[2].split('_')[0]
                cat['filter'] = filt
                
                if Cat is None:
                    Cat = cat
                else:
                    Cat = vstack([cat,Cat])
            
            if Cat is not None:
                Cat.write('catalog_'+dirtarget+'.fits',overwrite=1)
        elif verbose:

            
            ('catalog already exists, skipping')
            
    elif len(obs_table_cat)>0:
        obsids = obs_table_cat['obsid']
        products = Observations.get_product_list(obsids)
        starcat = [p for p in products['productFilename'] if p[-3:]=='cat']
        if verbose:
            print('Newly found catalogs:\n',starcat)
        if len(starcat)>0 and not(os.path.isfile(dirtarget+'/catalog_'+dirtarget+'.fits')):
            if verbose:
                print("Found catalogs:\n",starcat)
            if not(os.path.exists(dirtarget)):
                os.makedirs(dirtarget)
            
            os.chdir(dirtarget)
            if not(os.path.isfile('catalog_'+dirtarget+'.fits')):
                for uri in starcat['dataURI']:
                    fname = uri.split('/')[-1].split('filename=')[-1]
                    if not(os.path.isfile(fname)):
                        try:
                            Observations.download_file(uri)
                        except:
                            r=requests.get('https://hla.stsci.edu/cgi-bin/'+uri.split('/')[-1])
                            with open(fname,'wb') as fl:
                                fl.write(r.content)
                Cat = None
                sexnames=' X-Center    Y-Center         RA         DEC        ID     TotMag     MagErr       CI    Flags     MagAp1     MagErr1    FluxAp1    FluxErr1     MagAp2     MagErr2    FluxAp2    FluxErr2        Bck    FlagsSE        pix        deg     MagIso     MagErrIso    FluxIso    FluxErrIso  MagIsoCor     MagErrIsoCor FluxIsoCor    FluxErrIsocor   Radius    MagAuto     MagErrAuto   FluxAuto    FluxErrAuto    MagBest     MagErrBest   FluxBest    FluxErrBest     countsTh        MagTh     countsPk        MagPk          nl0          nl1          nl2          nl3          nl4          nl5          nl6          nl7       Xmin       Ymin       Xmax       Ymax         X2      ErrX2         Y2      ErrY2         XY      ErrXY        CXX     ErrCXX        CYY     ErrCYY        CXY     ErrCXY        A     ErrA        B     ErrB    ThetaPix  ErrThetaPix    FWHM             W_X             W_Y            W_X2         W_ErrX2            W_Y2         W_ErrY2            W_XY         W_ErrXY           W_CXX        W_ErrCXX           W_CYY        W_ErrCYY           W_CXY        W_ErrCXY             W_A          W_ErrA             W_B          W_ErrB    W_Theta W_ErrTheta     W_FWHM    Elong    Ellip    Theta ErrTheta          RA1950         DEC1950    Theta1950 ErrTheta1950          RA_sky         DEC_sky    Theta_sky ErrTheta_sky Star/Gal       c108       c109       c110       c111       c112       c113       c114       c115       c116       c117       c118       c119       c120       c121       c122       c123       c124       c125       c126       c127       c128       c129       c130       c131       c132      X-Det      Y-Det       Back       Cnts       XCTE       YCTE Chip    Mag-CTE'.split()
                sexnames2=' X-Center    Y-Center         RA         DEC        ID     TotMag     MagErr       CI    Flags     MagAp1     MagErr1    FluxAp1    FluxErr1     MagAp2     MagErr2    FluxAp2    FluxErr2        Bck    FlagsSE        pix        deg     MagIso     MagErrIso    FluxIso    FluxErrIso  MagIsoCor     MagErrIsoCor FluxIsoCor    FluxErrIsocor   Radius    MagAuto     MagErrAuto   FluxAuto    FluxErrAuto    MagBest     MagErrBest   FluxBest    FluxErrBest     countsTh        MagTh     countsPk        MagPk          nl0          nl1          nl2          nl3          nl4          nl5          nl6          nl7       Xmin       Ymin       Xmax       Ymax         X2      ErrX2         Y2      ErrY2         XY      ErrXY        CXX     ErrCXX        CYY     ErrCYY        CXY     ErrCXY        A     ErrA        B     ErrB    ThetaPix  ErrThetaPix    FWHM             W_X             W_Y            W_X2         W_ErrX2            W_Y2         W_ErrY2            W_XY         W_ErrXY           W_CXX        W_ErrCXX           W_CYY        W_ErrCYY           W_CXY        W_ErrCXY             W_A          W_ErrA             W_B          W_ErrB    W_Theta W_ErrTheta     W_FWHM    Elong    Ellip    Theta ErrTheta          RA1950         DEC1950    Theta1950 ErrTheta1950          RA_sky         DEC_sky    Theta_sky ErrTheta_sky Star/Gal       c108       c109       c110       c111       c112       c113       c114       c115       c116       c117       c118       c119       c120       c121       c122       c123       c124       c125       c126       c127       c128       c129       c130       c131       c132'.split()
                daonames='X-Center   Y-Center         RA          DEC        ID   MagAp1   MagErr1   MagAp2   MagErr2     MSky    Stdev       Flux   TotMag   MagErr       CI    Flags      X-Det      Y-Det       Back       Cnts     XCTE     YCTE Chip  Mag-CTE'.split()
                daonames2='X-Center   Y-Center         RA          DEC        ID   MagAp1   MagErr1   MagAp2   MagErr2     MSky    Stdev       Flux   TotMag   MagErr       CI    Flags'.split()
                catnames = ['X-Center','Y-Center','RA','DEC','ID','MagAp1','MagErr1','MagAp2','MagErr2','CI','Flags','ref']
                
                for fl in starcat['productFilename']:
                    if not(os.path.isfile(fl)):
                        continue
                    if fl[-4:]=='ecsv':
                        cat = Table.read(fl, format='ascii.ecsv')
                        if len(cat)==0:
                            if verbose:
                                print('empty table, skipping')
                            continue
                        cat['ref'] = 'pnt-cat'
                    else:
                        
                        cat = Table.read(fl, format='ascii.basic')
                        if len(cat)==0:
                            continue
                        if len(cat.colnames)==140:
                            cat.rename_columns(cat.colnames,sexnames)
                            cat['ref'] = 'sexphot'
                        elif len(cat.colnames)==132:
                            cat.rename_columns(cat.colnames,sexnames2)
                            cat['ref'] = 'sexphot'
                        elif len(cat.colnames)==24:
                            cat.rename_columns(cat.colnames,daonames)
                            cat['ref'] = 'daophot'
                        elif len(cat.colnames)==16:
                            cat.rename_columns(cat.colnames,daonames2)
                            cat['ref'] = 'daophot'
                        else:
                            if verbose:
                                print('failed to guess column names from column number, skip')
                            continue
                        cat.remove_columns([c for c in cat.colnames if not(c in catnames)])
                    
                    cat = cat[np.logical_and(abs(cat['RA']-(ramin+ramax)/2)<(ramax-ramin)/2,
                                             abs(cat['DEC']-(decmin+decmax)/2)<(decmax-decmin)/2)]
                    if len(cat)==0:
                        if verbose:
                            print('no match with galaxy, skipping')
                        continue
                    
                    cat.meta = None
                    filt = 'f'+fl.split('_f')[1].split('_')[0]
                    if len(filt)!=5:
                        filt = 'f'+fl.split('_f')[2].split('_')[0]
                    cat['filter'] = filt
                    
                    
                    if Cat is None:
                        Cat = cat
                    else:
                        Cat = vstack([cat,Cat])
                
                if Cat is not None:
                    Cat.write('catalog_'+dirtarget+'.fits',overwrite=1)
            elif verbose:
                print('catalog already exists, skipping')

    os.chdir(r'C:\Users\hgtra\Downloads\MAST_catalogs')

    print(dirtarget+',%d,%d,%d,%f,%f,'%(len(starcat),len(obs_table),len(flcproducts),redgreen_sel_moc.sky_fraction*41253*3600,moc_redgreen.sky_fraction*41253*3600)+col_v,file=f)
    

f.close()

### Download data products for selected observations
# Works: https://hla.stsci.edu/cgi-bin/getdata.cgi?dataset=hst_06699_05_wfpc2_total_pc_drz.fits

if dl_products:
    dirtarget = target.replace(' ','')
    if not(os.path.exists(dirtarget)):
        os.makedirs(dirtarget)
    
    os.chdir(dirtarget)
    
    for o in obs_table:
        Observations.download_products(o['obsid'],
                                       productType=["SCIENCE", "PREVIEW"],
                                       extension="fits")
        
os.chdir(r'C:\Users\hgtra\OneDrive\Documents\Python Scripts')