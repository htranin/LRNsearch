# -*- coding: utf-8 -*-
"""
Created on Tue Jun  6 07:53:41 2023

@author: Hugo
"""

from astropy.table import Table, vstack
from astropy.coordinates import SkyCoord
import glob
from tqdm import tqdm
import numpy as np
import os
import matplotlib.pyplot as plt
import sys


plt.rcParams.update({'font.size':14})

cleaning = 1
csvtofits = 0
fileout = "YSG_candidates_HSC_281123.fits"
stack = 0

#l = glob.glob("newCMD/*candidates.fits")
os.chdir("../../../Downloads/Catalogs/HSCv3/")
l = glob.glob("*candidates*_mist.fits")

print(len(l),'lists of candidates')
#gals = Table.read('../galaxy_summary_1123.csv')
gals = Table.read('../galaxy_summary_new.csv')
#gals2 = Table.read("../../../Downloads/Catalogs/galaxy_compilation_Simbad_unique_1023.fits")
gals2 = Table.read("../galaxies_hecate_20Mpc_new.fits")

gnames2 = [' '.join(g.replace('NAME','').split()) for g in gals2['main_id']]
#gals.remove_column(gals.colnames[1])

#gals = gals[np.asarray(gals['mstarref'])>0]
gnames = [' '.join(g.replace('NAME','').split()) for g in gals['main_id']]
ind2 = [gnames2.index(n) if n in gnames2 else -1 for n in gnames]
if csvtofits:
    gals['metal'].name = 'logMH_Guseva'
    gals['Av'].name = 'A0'
    gals['source'] = 'HSCv3'
    gals['source'][['MAST' in c for c in gals['commentaire']]] = 'MAST'
    gals['pct'] = 50.
    gals['pct'][['pct=' in c for c in gals['commentaire']]] = [float(c.replace(',',' ').split('pct=')[1].split()[0]) for c in gals['commentaire'] if 'pct=' in c]
    gals['Av'] = [gals2[i2]['A0'] if i2!=-1 else np.nan for i2 in ind2]
    gals['dAv'] = 0.
    gals['dAv'][['dAv=' in c for c in gals['commentaire']]] = [float(c.replace(',',' ').split('dAv=')[1].split()[0]) for c in gals['commentaire'] if 'dAv=' in c]
    gals['metal'] = [gals2[i2]["logMH_Guseva"] if i2!=-1 else np.nan for i2 in ind2]
    gals['dmetal'] = 0.
    gals['dmetal'][['dmetal=' in c for c in gals['commentaire']]] = [float(c.replace(',',' ').split('dmetal=')[1].split()[0]) for c in gals['commentaire'] if 'dmetal=' in c]
    gals['forcedist'] = 0.
    gals['forcedist'][['forcedist=' in c for c in gals['commentaire']]] = [float(c.replace(',',' ').split('forcedist=')[1].split()[0]) for c in gals['commentaire'] if 'forcedist=' in c]

    #gals.write('../galaxy_summary_1123.fits',overwrite=True)
    #gals.write('../galaxy_summary_new.fits',overwrite=True)

gals = Table.read('../galaxy_summary_new.fits')
gnames = [' '.join(g.replace('NAME','').split()) for g in gals['main_id']]
ind2 = [gnames2.index(n) if n in gnames2 else -1 for n in gnames]

lnames = [i.split('_')[0] for i in l]
#g,i1,i2 = np.intersect1d(lnames,gnames,return_indices=True)
ind = [gnames.index(n) if n in gnames else -1 for n in lnames]

#ind = [i for i in range(len(gals)) if l[i].split('_')[0] in gnames and int("MAST" in )+int("MAST" in gnames)]

if stack:
    #l = l[10:]
    t = Table.read(l[0])
    gal = lnames[0]
    if len(gal)<23:
        gal = " "*(23-len(gal))+gal
    print("First galaxy:",gal)
    
    t.add_column(gal,name='galaxy',index=0)
    t.add_column(99.,name='distance',index=1)
    t.add_column(99.,name='Av',index=2)
    t.add_column(-99.,name='metal',index=3)
    t.add_column(99.,name='pctGM',index=4)
    t.add_column(99.,name='forcedist',index=5)
    t.add_column('        ',name='filter',index=6)
    t.add_column('HSCv3',name='source',index=6)
    
    
        
    t = t[:0] # keep only table structure, required to stack
    
    ncand = []
    for i in tqdm(range(0,len(l))):
        if ind[i]==-1:
            # galaxy not in selection
            continue
        if lnames[i]!="M 31" and int(gals[ind[i]]['source']=="MAST")-int('MAST' in l[i])!=0:
            # HSCv3 catalog for MAST galaxy or the other way around: skip
            continue
            
        cands = Table.read(l[i])
        if len(cands)==0:
            continue
        if not('x' in cands.colnames):
            print('skip',l[i])
            continue
        
        gal = lnames[i]
        cands['galaxy'] = gal
        cands['distance'] = gals[ind[i]]['Dist']
        r1,r2,pa = gals[ind[i]]['R1','R2','PA']
        ra1 = np.asarray(cands['MatchRA'])
        dec1 = np.asarray(cands['MatchDec'])
        ra2,dec2 = gals[ind[i]]['ra','dec']
        #r1,r2 = 1.26*r1, 1.26*r2
    
        c1 = SkyCoord(ra=ra1,dec=dec1,frame="icrs",unit="deg") # coordinates of TNS objects
        c2 = SkyCoord(ra=ra2,dec=dec2,frame="icrs",unit="deg") # coordinates of galaxies
        d2d = c2.separation(c1)
        posang = c2.position_angle(c1).deg-pa
    
        posmatch = (d2d.arcsec < r1*r2/np.sqrt(r1**2*np.sin(posang*np.pi/180)**2+r2**2*np.cos(posang*np.pi/180)**2))
        # compare each separation to the corresponding galaxy radius at position angle
        cands = cands[posmatch.reshape(-1)]
    
        
        if len(cands)==0:
            print("empty candidates / galaxy intersection, skip")
            continue
        
        cands['Av'] = gals[ind[i]]['Av']+gals[ind[i]]['dAv']
        cands['metal'] = gals[ind[i]]['metal']+gals[ind[i]]['dmetal']
        cands['pctGM'] = gals[ind[i]]['pct']
        cands['forcedist'] = gals[ind[i]]['forcedist']
        cands['source'] = gals[ind[i]]['source']
        
        if "W_" in l[i]:
            # add instrument and filter
            cands['filter'] = '_'.join(l[i].split('_')[-3:-1])
        else:
            cands['filter'] = '        '
        # flattening the catalog
        for c in cands.colnames:
            if len(cands[c].shape)>1:
                cands[c] = cands[c][:,0]
                
        ncand.append(len(cands))
        t = vstack((t,cands))
        
    
    t.write(fileout, overwrite=True)
    
else:
    t = Table.read(fileout)
    

u,nu = np.unique(t['galaxy'],return_index=True)
nu = np.array(list(nu)+[len(t)])
ncand = nu[1:]-nu[:-1]
    
plt.hist(ncand,bins=np.geomspace(1,1e4,50),color='C3')
    #plt.axvline(t['Pextin'].mean(),color='k',label='all galaxies')
plt.gca().set_xscale('log')
#plt.legend()
plt.xlabel('Number of candidates')
plt.ylabel('N galaxies')
plt.tight_layout()
plt.savefig('../hist_ncand.png')
#sys.exit()

if cleaning:

    # crossmatch with Ga#ia to remove stars / foreground objects
    filegaia = fileout.replace('.fits','_xmatch.fits')
    cmd = 'java -jar C:/Users/hgtra/OneDrive/Documents/stilts.jar cdsskymatch in="%s" ra="MatchRA" dec="MatchDec" cdstable="Gaia DR3 (Epoch 2016)" find=all radius=0.5 out="%s"'%(fileout,filegaia)
    print(cmd)
    os.system(cmd)
    cand_gaia = Table.read(filegaia)
    t['Gaia_star'] = np.nan
    pm_snr = np.maximum(abs(cand_gaia['pmRA']/cand_gaia['e_pmRA']),
                        abs(cand_gaia['pmDE']/cand_gaia['e_pmDE']))
    highpm = pm_snr > 4
    highps = np.logical_and(cand_gaia['PSS']>0.9999,np.logical_or(highpm,cand_gaia['distance']>1))
    ind_match = np.intersect1d(np.asarray(t['MatchID']),np.asarray(cand_gaia['MatchID']),return_indices=True)
    t['Gaia_star'][ind_match[1]] = 0 # non-nan: there is a Gaia source
    
    print("Found",len(cand_gaia),"Gaia counterparts")
    cand_gaia = cand_gaia[np.logical_or(highpm,highps)]
    print("Including",len(cand_gaia),"Gaia stars")
    
    ind_match = np.intersect1d(np.asarray(t['MatchID']),np.asarray(cand_gaia['MatchID']),return_indices=True)
    t['Gaia_star'][ind_match[1]] = cand_gaia[ind_match[2]]['Source_cds']
    # non-zero: Gaia source is a star
    t.write(fileout,overwrite=1)
    
    
    # crossmatch with Simbad
    filesimb = fileout.replace('.fits','_xmatch.fits')
    cmd = 'java -jar C:/Users/hgtra/OneDrive/Documents/stilts.jar cdsskymatch in="%s" ra="MatchRA" dec="MatchDec" cdstable="simbad" find=all radius=0.5 out="%s"'%(fileout,filesimb)
    print(cmd)
    os.system(cmd)
    cand_simb = Table.read(filesimb)
    # extended contaminants
    conta1 = np.logical_or(np.logical_or(np.logical_or(cand_simb['main_type']=='Cl*',
                                                       cand_simb['main_type']=='SNR'),
                                         np.logical_or(cand_simb['main_type']=='GlCl',
                                                       cand_simb['main_type']=='Assoc*')),
                           np.logical_or(np.logical_or(cand_simb['main_type']=='QSO',
                                                       cand_simb['main_type']=='SN'),
                                         np.logical_or(cand_simb['main_type']=='Galaxy',
                                                       cand_simb['main_type']=='LINER')))
    # point-like contaminants
    conta2 = np.logical_or(np.logical_or(np.logical_or(cand_simb['main_type']=='HMXB',
                                                       cand_simb['main_type']=='LMXB'),
                                         np.logical_or(cand_simb['main_type']=='deltaCep',
                                                       cand_simb['main_type']=='Cepheid')),
                           np.logical_or(np.logical_or(cand_simb['main_type']=='WR*',
                                                       cand_simb['main_type']=='RedSG*'),
                                         np.logical_or(cand_simb['main_type']=='BlueSG*',
                                                       cand_simb['main_type']=='Nova')))
    
    cand_simb = cand_simb[np.logical_or(conta1,conta2)]
    print("Found",len(cand_simb),"Simbad contaminants")
    t['Simbad_contam'] = "        "
    ind_match = np.intersect1d(np.asarray(t['MatchID']),np.asarray(cand_simb['MatchID']),return_indices=True)
    t['Simbad_contam'][ind_match[1]] = cand_simb[ind_match[2]]['main_type']
    t.write(fileout,overwrite=1)
    
    # crossmatch with Chandra
    filexray = fileout.replace('.fits','_xmatch.fits')
    csc_cat = r"C:\Users\hgtra\Downloads\Catalogs\CSC2.1_4XMMDR13_accurate1arcsec.fits"
    cmd = 'java -jar C:/Users/hgtra/OneDrive/Documents/stilts.jar tmatch2 matcher="skyerr" params="1" in1="%s" in2="%s" values1="MatchRA MatchDec 0" values2="ra dec poserr" find="best1" out="%s"'%(fileout,csc_cat,filexray)
    print(cmd)
    os.system(cmd)
    cand_xray = Table.read(filexray)
    print("Found",len(cand_xray),"X-ray counterparts")
    t['XRB_sep'] = np.nan
    ind_match = np.intersect1d(np.asarray(t['MatchID']),np.asarray(cand_xray['MatchID']),return_indices=True)
    t['XRB_sep'][ind_match[1]] = cand_xray[ind_match[2]]['Separation']
    
    plt.figure(figsize=(7,5))
    logFxFopt = np.log10(cand_xray[ind_match[2]]['flux'])+5.37+t['A_F814W'][ind_match[1]]/2.5
    i0 = np.logical_or(np.logical_or(np.isnan(np.asarray(cand_xray[ind_match[2]]['flux'])),
                                     np.isnan(np.asarray(t['A_F814W'][ind_match[1]]))),
                       cand_xray[ind_match[2]]['flux']==0)
    
    logFxFopt = logFxFopt[~i0]
    plt.hist(logFxFopt,bins=25,histtype='step',color='k')

    plt.xlabel("log(F_X/F_I)")
    plt.ylabel('Number of XRB candidates')
    
    plt.tight_layout()
    plt.legend()
    plt.savefig('C:/Users/hgtra/OneDrive/Documents/figures/hist_logFxFopt_YSG_xray.png')

    
    t.write(fileout,overwrite=1)
    
    
    # print info
    for typ in ['Cl*', 'SNR', 'GlCl', 'Assoc*', 'QSO', 'SN', 'Galaxy', 'LINER', 'HMXB',
                'LMXB', 'deltaCep', 'Cepheid', 'WR*', 'RedSG*', 'BlueSG*','Nova']:
        print("Number of",typ,np.sum(t['Simbad_contam']==typ))

# assign galaxy distances
#cand = cand[np.argsort(cand['galaxy'])]
#gals = Table.read('../../../Downloads/Catalogs/galaxy_hst_0723.fits')
#gals['name']=[" ".join(n.replace('NAME','').split()) for n in gals['main_id']]
#igal = np.intersect1d(cand['galaxy'],gals['name'],return_indices=True)
#igal1 = np.concatenate((igal[1],[len(cand)]))
#cand['Dist'] = np.nan
#for i in range(len(igal1)-1):
#    cand['Dist'][igal1[i]:igal1[i+1]] = gals[igal[2][i]]['Dist']
#cand.write(fileout,overwrite=1)

# plot results

#sys.exit()

plt.figure(figsize=(7,5))
plt.hist(cand_xray['Separation']*cand_xray['poserr'],bins=50,histtype="step",label="CSC2")
plt.hist(cand_gaia['angDist'],bins=50,histtype="step",label="Gaia")
plt.hist(cand_simb['angDist'],bins=50,histtype="step",label="Simbad")

plt.xlabel("Separation (arcsec)")
plt.ylabel('Number of counterparts')
plt.gca().set_yscale('log')
plt.tight_layout()
plt.legend()
plt.savefig('C:/Users/hgtra/OneDrive/Documents/figures/hist_sep_contaminants.png')

gal,ugal = np.unique(t['galaxy'],return_index=1)
ugal = np.concatenate((ugal,[len(t)]))
npergal = ugal[1:]-ugal[:-1]
plt.figure(figsize=(7,5))
plt.hist(npergal,bins=np.geomspace(0.8,1000,20),color='0.5')
plt.gca().set_xscale('log')
plt.xlabel('Number of candidates per galaxy')
plt.ylabel('Number of galaxies')
plt.savefig('C:/Users/hgtra/OneDrive/Documents/figures/hist_npergal_cand2.png')

plt.figure(figsize=(7,5))
plt.hist(t['distance'],bins=np.geomspace(0.3,16,20),color='0.5')
plt.gca().set_xscale('log')
plt.xlabel('Distance (Mpc)')
plt.ylabel('Number of candidates')
plt.savefig('C:/Users/hgtra/OneDrive/Documents/figures/hist_dist_cand2.png')

plt.figure(figsize=(7,5))
plt.hist(t['A_F814W']-5*np.log10(t['distance']*1e6)+5,bins=np.linspace(-13,1,20),color='0.5')
plt.xlim((1,-13))
plt.ylim((0,1.3*plt.ylim()[1]))
plt.xlabel('Absolute magnitude A_F814W')
plt.ylabel('Number of candidates')
plt.savefig('C:/Users/hgtra/OneDrive/Documents/figures/hist_F814W_abs_cand2.png')

plt.figure(figsize=(7,5))
gaia = ~np.isnan(t['Gaia_star'])
simbad = np.array([c.strip()!="" for c in t['Simbad_contam']])
plt.hist(t['A_F814W'],bins=np.linspace(16,26,20),color='0.5',label='YSG candidates')
plt.hist((t['A_F814W'])[gaia],bins=np.linspace(16,26,20),color='C6',label='Gaia contaminants')
plt.hist((t['A_F814W'])[simbad],bins=np.linspace(16,26,20),color='C3',label='Simbad contaminants')
plt.xlim((26,16))
plt.ylim((0,1.1*plt.ylim()[1]))
plt.legend()
plt.xlabel('Observed magnitude A_F814W')
plt.ylabel('Number of candidates')
plt.savefig('C:/Users/hgtra/OneDrive/Documents/figures/hist_F814W_cand2.png')

plt.figure(figsize=(5,5))
hc = np.histogram((t['A_F814W'])[np.logical_or(gaia,simbad)],bins=np.linspace(16,26,21))[0]
h0 = np.histogram((t['A_F814W']),bins=np.linspace(16,26,21))[0]
plt.errorbar(np.linspace(16.25,25.75,20)[h0>3],(hc/h0)[h0>3],yerr=(np.sqrt(hc)/h0)[h0>3],color='k')
plt.xlabel('Observed magnitude A_F814W')
plt.ylabel('Fraction of contaminants')
plt.xlim((26,16))
plt.tight_layout()
plt.savefig('C:/Users/hgtra/OneDrive/Documents/figures/fcont_F814W_cand2.png')


# CMD of candidates with data for a large galaxy as background
print('Galaxy with maximum number of candidates (%d):'%(npergal.max()),gal[npergal.argmax()])
gmax = np.argsort(npergal)[-14] # index of maximum, fine-tunable. Ex -1, -2, -14
filename = "C:/Users/hgtra/Downloads/Catalogs/HSCv3/%s_data.fits"%(gal[gmax])
data = Table.read(filename)
i = np.argmax([len(np.where(~np.isnan(data['A_F475W']))[0]),
               len(np.where(~np.isnan(data['A_F555W']))[0]),
               len(np.where(~np.isnan(data['A_F606W']))[0]),
               len(np.where(~np.isnan(data['W3_F475W']))[0]),
               len(np.where(~np.isnan(data['W3_F555W']))[0]),
               len(np.where(~np.isnan(data['W3_F606W']))[0]),
               len(np.where(~np.isnan(data['W2_F555W']))[0]),
               len(np.where(~np.isnan(data['W2_F606W']))[0])])

col_v = ['A_F475W','A_F555W','A_F606W',
         'W3_F475W','W3_F555W','W3_F606W',
         'W2_F555W','W2_F606W'][i]

dist = t[ugal[gmax]]['distance']
y = data[col_v]
x = y-data['A_F814W']
y += -5*np.log10(dist*1e6)+5 # absolute magnitude
x[np.abs(x)>10]=np.nan
data = data[~np.isnan(x)]
y = y[~np.isnan(x)]
x = x[~np.isnan(x)]

X = t[col_v]-t['A_F814W']
Y = t[col_v]-5*np.log10(t['distance']*1e6)+5

    
plt.figure(figsize=(10,7))
plt.scatter(x,y,s=1,alpha=0.05)
plt.scatter(x[0],y[0],s=10,alpha=0.2,label=gal[gmax],c="C0")
plt.scatter(X,Y,s=10,alpha=0.2,label='YSG candidates (all galaxies)')
plt.scatter(X[gaia],Y[gaia],s=50,alpha=0.5,label='Gaia contaminants', c="C6", marker='*')
plt.scatter(X[simbad],Y[simbad],s=30,alpha=0.5,label='Simbad contaminants', c="C3", marker='s')
if plt.xlim()[0]<-2:
    plt.xlim((-2,plt.xlim()[1]))
if plt.xlim()[1]>3.5:
    plt.xlim((plt.xlim()[0],3.5))
plt.ylim((plt.ylim()[1],plt.ylim()[0]))

plt.xlabel(col_v+' - A_F814W')
plt.ylabel(col_v+' (absolute magnitude)')
plt.legend()

plt.savefig('C:/Users/hgtra/OneDrive/Documents/figures/CMD_all_candidates_%s_2.png'%col_v)

ibright=np.logical_and(abs(t['A_F475W']-t['A_F814W']-0.45)<0.15,t['A_F475W']-5*np.log10(t['distance']*1e6)+5<-7)

t['contam'] = len(t)*[0]
t['contam'][np.logical_or(simbad,gaia)] = 1
print("Bright candidates F475W:\n",t[ibright]['A_F475W','galaxy','MatchRA','MatchDec','contam'])



