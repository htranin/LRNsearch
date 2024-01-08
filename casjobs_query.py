# -*- coding: utf-8 -*-
"""
Created on Mon May 29 14:04:24 2023

@author: Hugo
"""

# Based on the notebook https://archive.stsci.edu/hst/hsc/help/api/hscv3_smc_api.html

#import astropy
from astropy.coordinates import SkyCoord
import time
#import sys
import os
#import requests
#import json
import numpy as np
#import matplotlib.pyplot as plt

#from PIL import Image
#from io import BytesIO

from astropy.table import Table

# set universal matplotlib parameters
#plt.rcParams.update({'font.size': 16})



tget = "NGC5204" # if "": whole catalog


HSCContext = "HSCv3"
# check that version of mastcasjobs is new enough
# we are using some features not in version 0.0.1
from pkg_resources import get_distribution
from packaging.version import Version as V

assert V(get_distribution("mastcasjobs").version) >= V('0.0.2'), """
A newer version of mastcasjobs is required.
Update mastcasjobs to current version using this command:
pip install --upgrade git+git://github.com/rlwastro/mastcasjobs@master
"""

import mastcasjobs

import getpass
if not os.environ.get('CASJOBS_USERID'):
    os.environ['CASJOBS_USERID'] = input('Enter Casjobs UserID:')
if not os.environ.get('CASJOBS_PW'):
    os.environ['CASJOBS_PW'] = getpass.getpass('Enter Casjobs password:')


def run_query(tget, verb=False):
    gals = Table.read("../../../Downloads/Catalogs/galaxy_compilation_Simbad_unique_1023.fits")
    gals = Table.read("../../../Downloads/Catalogs/galaxies_hecate_20Mpc_new.fits")
    #gals = Table.read("../../../Downloads/Catalogs/galaxy_hst_0723.fits")
    #gals = Table.read("sample_casjobs_reallyunique.fits")
    
    reduce = 0 # if reduce: query central 1 arcmin
    
    f = open('log_casjobs_query.txt','a')


    for gal in gals:
        if tget and gal['main_id'].replace("NAME","").replace(' ','')!=tget:
            continue
        f.close()
        f = open('log_cmd_master.txt','a')
        try:
            target = " ".join(gal['main_id'].replace('NAME','').split())
            if os.path.isfile('../Catalogs/HSCv3/'+target+"_data.fits") and 0:
                continue
            #target = "NGC 4395"
            #gal = gals[[target.replace(" ","")=="".join(g['main_id'].split()).replace("NAME","") for g in gals]]
            
            coord_g = SkyCoord.from_name(target)
            
            ra = coord_g.ra.degree
            dec = coord_g.dec.degree
            print(target, f'ra: {ra}\ndec: {dec}')
            
            DBtable = "HSC_galaxy"
            
            
            # drop table if it already exists
            jobs = mastcasjobs.MastCasJobs(context="MyDB")
            jobs.drop_table_if_exists(DBtable)
            
            
            
            # get main information
            ddec = gal['R1']/3600 # deg
            if reduce:
                ddec = 1./60
            dra = gal['R1']/np.cos(np.radians(gal['dec']))/3600
            
            
            if ddec>0.3:
                # tiling: nit**2 tiles
                nit = int(ddec//0.3+1)
                nit = 12
                print(gal['main_id'],nit**2,"tiles")
                for i in range(nit):
                    for j in range(nit):
                        DBtable = "HSC_galaxy_%d%d"%(i,j)
                        if abs(i-j)>=4 or i<6 or i==6 and j<=6:
                            continue
                        ra_min, ra_max = "%.4f"%(ra-dra+2*i*dra/nit), "%.4f"%(ra-dra+2*(i+1)*dra/nit),
                        dec_min, dec_max = "%.4f"%(dec-ddec+2*j*ddec/nit), "%.4f"%(dec+-ddec+2*(j+1)*ddec/nit)
                        # test: only central tile
                        if ra<float(ra_min) or ra>float(ra_max) or dec<float(dec_min) or dec>float(dec_max):
                            continue
                        print(ra_min,ra_max,dec_min,dec_max)
                        
                        jobs = mastcasjobs.MastCasJobs(context="MyDB")
                        jobs.drop_table_if_exists(DBtable)
                        
                        query = f"""
                    select 
                        p.MatchID, p.MatchRA, p.MatchDec, p.DSigma, p.AbsCorr,
                        p.Extinction, p.SpectrumFlag,
                        p.CI, p.CI_Sigma, 
                        c.W3_F336W,  c.W3_F336W_MAD,  c.W3_F336W_N,
                        c.A_F435W,  c.A_F435W_MAD,  c.A_F435W_N,
                        c.W3_F438W,  c.W3_F438W_MAD,  c.W3_F438W_N,
                        c.W2_F439W,  c.W2_F439W_MAD,  c.W2_F439W_N,
                        c.A_F475W,  c.A_F475W_MAD,  c.A_F475W_N,
                        c.W3_F475W, c.W3_F475W_MAD, c.W3_F475W_N,
                        c.A_F555W,  c.A_F555W_MAD,  c.A_F555W_N,
                        c.W3_F555W, c.W3_F555W_MAD, c.W3_F555W_N,
                        c.W2_F555W, c.W2_F555W_MAD, c.W2_F555W_N,
                        c.A_F606W,  c.A_F606W_MAD,  c.A_F606W_N,
                        c.W3_F606W, c.W3_F606W_MAD, c.W3_F606W_N,
                        c.W2_F606W, c.W2_F606W_MAD, c.W2_F606W_N,
                        c.A_F814W,  c.A_F814W_MAD,  c.A_F814W_N,
                        c.W3_F814W, c.W3_F814W_MAD, c.W3_F814W_N,
                        c.W2_F814W, c.W2_F814W_MAD, c.W2_F814W_N
                    into MyDB.{DBtable}
                    from dbo.fHtmCoverRegion('POLY J2000 {ra_min} {dec_min} {ra_min} {dec_max} {ra_max} {dec_max} {ra_max} {dec_min}') h
                    join SumPropMagAper2Cat p on p.htmid between h.HtmIDStart and h.HtmIDEnd
                    join SumMagAper2Cat c on c.MatchID=p.MatchID
                    where
                        coalesce(c.W3_F336W, c.A_F435W, c.W3_438W, c.W2_439W, c.W3_F475W, c.A_F475W, c.A_F555W, c.W3_F555W, c.W2_F555W, c.A_F606W, c.W3_F606W, c.W2_F606W) is not null and
                        coalesce(c.A_F814W, c.W3_F814W, c.W2_F814W) is not null and
                        p.MatchRA between {ra_min} and {ra_max} and
                        p.MatchDec between {dec_min} and {dec_max}
                    """
                    
                        
                        t0 = time.time()
                        results = jobs.quick(query, task_name="HSC galaxy", context=HSCContext)
                        if results[0][0]==0:
                            continue
                        print(f"Completed in {(time.time()-t0):.1f} sec")
                        print(target,"\n",results, file=f)
                        print(target,"\n",results)
                        
                        # fast retrieval using special MAST Casjobs service
                        try:
                            tab = jobs.get_table(DBtable, verbose=True)
                            tab.write("../../../Downloads/Catalogs/HSCv3/"+target+"_%d%d_data.fits"%(i,j),overwrite=1)
                            print("output successfully written")
                            jobs.drop_table_if_exists(DBtable)
                        except Exception as e:
                            print(e)
                            print("writing skipped")
                            continue
                        
                        #plt.scatter(tab['MatchRA'],tab['MatchDec'],alpha=0.5,c="k",s=1)
                    
            else:
                jobs = mastcasjobs.MastCasJobs(context="MyDB")
                jobs.drop_table_if_exists(DBtable)
                ra_min, ra_max = "%.4f"%(ra-dra), "%.4f"%(ra+dra),
                dec_min, dec_max = "%.4f"%(dec-ddec), "%.4f"%(dec+ddec)
    
                query = f"""
            select 
                p.MatchID, p.MatchRA, p.MatchDec, p.DSigma, p.AbsCorr,
                p.Extinction, p.SpectrumFlag,
                p.CI, p.CI_Sigma, 
                c.W3_F336W,  c.W3_F336W_MAD,  c.W3_F336W_N,
                c.A_F435W,  c.A_F435W_MAD,  c.A_F435W_N,
                c.W3_F438W,  c.W3_F438W_MAD,  c.W3_F438W_N,
                c.W2_F439W,  c.W2_F439W_MAD,  c.W2_F439W_N,
                c.A_F475W,  c.A_F475W_MAD,  c.A_F475W_N,
                c.W3_F475W, c.W3_F475W_MAD, c.W3_F475W_N,
                c.A_F555W,  c.A_F555W_MAD,  c.A_F555W_N,
                c.W3_F555W, c.W3_F555W_MAD, c.W3_F555W_N,
                c.W2_F555W, c.W2_F555W_MAD, c.W2_F555W_N,
                c.A_F606W,  c.A_F606W_MAD,  c.A_F606W_N,
                c.W3_F606W, c.W3_F606W_MAD, c.W3_F606W_N,
                c.W2_F606W, c.W2_F606W_MAD, c.W2_F606W_N,
                c.A_F814W,  c.A_F814W_MAD,  c.A_F814W_N,
                c.W3_F814W, c.W3_F814W_MAD, c.W3_F814W_N,
                c.W2_F814W, c.W2_F814W_MAD, c.W2_F814W_N
            into MyDB.{DBtable}
            from dbo.fHtmCoverRegion('POLY J2000 {ra_min} {dec_min} {ra_min} {dec_max} {ra_max} {dec_max} {ra_max} {dec_min}') h
            join SumPropMagAper2Cat p on p.htmid between h.HtmIDStart and h.HtmIDEnd
            join SumMagAper2Cat c on c.MatchID=p.MatchID
            where
                coalesce(c.W3_F336W, c.A_F435W, c.W3_F438W, c.W2_F439W, c.W3_F475W, c.A_F475W, c.A_F555W, c.W3_F555W, c.W2_F555W, c.A_F606W, c.W3_F606W, c.W2_F606W) is not null and
                coalesce(c.A_F814W, c.W3_F814W, c.W2_F814W) is not null and
                p.MatchRA between {ra_min} and {ra_max} and
                p.MatchDec between {dec_min} and {dec_max}
            """
            
            
                t0 = time.time()
                results = jobs.quick(query, task_name="HSC galaxy", context=HSCContext)
                if results[0][0]==0:
                    print(query)
                    continue
                    
                print(f"Completed in {(time.time()-t0):.1f} sec")
                print(target,"\n",results, file=f)
                print(target,"\n",results)
                
                # fast retrieval using special MAST Casjobs service
                tab = jobs.get_table(DBtable, verbose=True)
                if len(tab)==results[0][0]:
                    tab.write("../../../Downloads/Catalogs/HSCv3/"+target+"_data.fits",overwrite=1)
                else:
                    print("Failed to retrieve %s, please download file from CasJobs website (https://mastweb.stsci.edu/mcasjobs/output.aspx)"%(target+"_data.fits"))
                    
                    h = open("GALAXIES_TO_DOWNLOAD_ON_CASJOBS.txt","a")
                    h.write(target+"\n")
                    h.close()
                    return None
                
                #plt.scatter(tab['MatchRA'],tab['MatchDec'],alpha=0.5,c="k",s=1)
            if verb:
                return(query)
            else:
                return("")
        
        except Exception as e:
            print(target,"\n",e, file=f)
            print(target,"\n",e)
            continue


if __name__=="__main__":
    print(run_query(tget))