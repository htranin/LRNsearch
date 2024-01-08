# -*- coding: utf-8 -*-
"""
Created on Wed Jul  5 08:56:31 2023

@author: hgtra

Retrieves the MAST color image of the target and plots it in DS9
"""



from astroquery.mast import Observations
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
from mocpy import MOC
import os
import sys
from astropy.table import Table, vstack
import requests
import subprocess
from requests.exceptions import Timeout
import glob
from astropy.io import fits

if len(sys.argv)>1:
    target = sys.argv[1]
else:
    target = "NGC 7552"
    
plot = '' # 'candidates', 'all' or ''    

dirtarget = "../../../Downloads/MAST_catalogs/"+target.replace(' ','')
coord = SkyCoord.from_name(target)
ra,dec = coord.ra.deg, coord.dec.deg
#gals = Table.read('../../../Downloads/Catalogs/galaxy_compilation_Simbad_unique_0623.fits')
gals = Table.read('../../../Downloads/Catalogs/galaxies_hecate_20Mpc.fits')
gal = gals[np.array([target.replace(' ','')==name.replace(' ','') for name in gals['main_id']])]
r1 = gal['R1']/3600.


# Query the MAST database for >300s observations
obs_table = Observations.query_object(target, radius="%f deg"%r1)
obs_table = obs_table[obs_table['t_exptime']>300]

# find HST color images
obs_table = obs_table[['hst_' in l and 'red=' in l for l in obs_table['dataURL'] ]]

if len(obs_table)==0:
    print("no HST color image, exit")
    sys.exit()
else:
    print(len(obs_table),'color images')
    obs_table['image'] = ['https://hla.stsci.edu/cgi-bin/fitscut.cgi?'+url.split('cgi?')[-1].strip() for url in obs_table['dataURL']]
    

# add information on observations coverage
mocs = []
covs = []

for o in obs_table:
    m = None
    for s_region in o['s_region'].upper().replace(' J2000','').split('POLYGON'):
        if len(s_region)==0: #empty string
            continue
        vertices = np.array(s_region.split()).astype(float).reshape((-1,2))
        if m is None:
            m = MOC.from_polygon(vertices[:,0]*u.deg,vertices[:,1]*u.deg,max_depth=16)
        else:
            m = m.union(MOC.from_polygon(vertices[:,0]*u.deg,vertices[:,1]*u.deg,max_depth=16))
    mocs.append(m)
    covs.append(float("%.2e"%(m.sky_fraction*41253*3600))) #squared arcmin

obs_table['coverage'] = covs

# sort by coverage
obs_table = obs_table[np.argsort(obs_table['coverage'])[::-1]]

obs_table['img']=['hst'+url.split('=hst')[1][:-10]+'_'+url.split('=hst')[2][-18:-9]+'_'+url.split('size')[0][-14:-5] for url in obs_table['image']]

print(obs_table['target_name','img'])

# download first observation


if not(os.path.exists(dirtarget)):
    os.makedirs(dirtarget)

os.chdir(dirtarget)

l = glob.glob('*hst*.fits')
if len(l)>0:
    print('found',len(l),'images in local path:\n',l)
    for i in range(len(l)):
        fname = l[i]
        if 1:
            rinp = input('use %s... (type n for break)? '%(fname[:30]))
            if rinp=='n':
                break
            elif rinp!='y':
                continue
        
        r,g,b = fits.getdata(fname).reshape((3,-1))
        mins,maxs = [],[]
        for c in [r,g,b]:
            cmin,cmax = np.nanmin(c[::1000]),np.nanmax(c)/4
            maxs.append(str(cmax)[:7])
            h = np.histogram(c[::1000]-cmin+0.1,bins=np.geomspace(0.1,cmax-cmin+0.1,1000))
            cmin = h[1][h[0].argmax()]/2+cmin/2-0.1
            # new cmin = middlepoint of c_hmax and cmin
            mins.append(str(cmin)[:7])
            print('color new scale min max:',cmin,cmax)
            
        plot_fname = ''
        if plot=='candidates':
            plot_fname = "../../Catalogs/HSCv3/%s_candidates2.fits"%target
            if not(os.path.isfile(plot_fname)):
                plot_fname = "../../Catalogs/%s_candidates2.fits"%target
            filt_fname = "../filter_cand.txt"
            c_ra,c_dec = "MatchRA","MatchDec"
        elif plot=='all':
            plot_fname = "catalog_"+dirtarget+".fits"
            filt_fname = "../filter.txt"
            c_ra,c_dec = "RA","DEC"
        
        
        if plot and os.path.isfile(plot_fname):
            subprocess.call([r'C:\Users\hgtra\ds9.exe','-scale','log','-rgb',
                             '-red',fname,'-scale','limits',mins[0],maxs[0],
                             '-green',fname,'-cube','2','-scale','limits',mins[1],maxs[1],
                             '-blue',fname,'-cube','3','-scale','limits',mins[2],maxs[2],
                             '-rgb','lock','colorbar','yes',
                             '-cmap','7','0.09',
                             "-catalog","import","fits",plot_fname,
                             "-catalog","ra",c_ra,
                             "-catalog","dec",c_dec,
                             "-catalog","filter","load",filt_fname,
                             "-catalog","update"])
        else:
            subprocess.call([r'C:\Users\hgtra\ds9.exe','-scale','log','-rgb',
                             '-red',fname,'-scale','limits',mins[0],maxs[0],
                             '-green',fname,'-cube','2','-scale','limits',mins[1],maxs[1],
                             '-blue',fname,'-cube','3','-scale','limits',mins[2],maxs[2],
                             '-rgb','lock','colorbar','yes',
                             '-cmap','7','0.09'])
            
        os.chdir('..')
        sys.exit()


for i in range(len(obs_table)):
    url = obs_table['image'][i]
    fname = url[:-29].split('cgi?')[-1].replace('&amp;','_').replace('=','_')+'.fits'
    
    if 1:
        rinp = input('download (%s...) (type n for break)? '%(fname[:30]))
        if rinp=='n':
            break
        elif rinp!='y':
            continue
    try:
        r = requests.get(url,timeout=(30,60))
    except Timeout:
        rinp = input("request timed out after 20s, start a longer request? [y] ")
        if rinp=='n':
            continue
        r = requests.get(url)
    
    
    with open(fname,'wb') as fl:
        fl.write(r.content)
    print('file successfully written')
    
    r,g,b = fits.getdata(fname).reshape((3,-1))
    mins,maxs = [],[]
    for c in [r,g,b]:
        cmin,cmax = np.nanmin(c[::1000]),np.nanmax(c)/4
        maxs.append(str(cmax)[:7])
        h = np.histogram(c[::1000]-cmin+0.1,bins=np.geomspace(0.1,cmax-cmin+0.1,1000))
        cmin = h[1][h[0].argmax()]/2+cmin/2-0.1
        # new cmin = middlepoint of c_hmax and cmin
        mins.append(str(cmin)[:7])
        print('color new scale min max:',cmin,cmax)
        
    
    subprocess.call([r'C:\Users\hgtra\ds9.exe','-scale','log','-rgb',
                     '-red',fname,'-scale','limits',mins[0],maxs[0],
                     '-green',fname,'-cube','2','-scale','limits',mins[1],maxs[1],
                     '-blue',fname,'-cube','3','-scale','limits',mins[2],maxs[2],
                     '-rgb','lock','colorbar','yes',
                     '-cmap','7','0.09'])
    
os.chdir('..')
        

