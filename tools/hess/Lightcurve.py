#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import json
import os
import shutil
import sys

try:
    import numpy as np
    _numpy_available = True
except ImportError:
    _numpy_available = False

try:
    from oda_api.json import CustomJSONEncoder
except ImportError:
    from json import JSONEncoder as CustomJSONEncoder

_galaxy_wd = os.getcwd()


# In[1]:


import astropy.units as u
from astropy.coordinates import SkyCoord

import matplotlib.pyplot as plt
from gammapy.data import DataStore
from gammapy.makers import MapDatasetMaker
from gammapy.makers.utils import make_theta_squared_table
from gammapy.maps import Map, MapAxis, WcsGeom
from astropy.io import fits
import numpy as np
from numpy import pi,cos,sin,sqrt,log10
import os
from astropy import wcs
from astropy.io import fits
from oda_api.data_products import PictureProduct
from oda_api.data_products import ODAAstropyTable
from oda_api.data_products import ImageDataProduct
from astropy.time import Time


# In[11]:


hess_data="gammapy-datasets/1.1/hess-dl3-dr1/"
if not(os.path.exists(hess_data)):
    get_ipython().system('gammapy download datasets')
    
data_store = DataStore.from_dir(hess_data)


# if(os.path.exists('hess_dl3_dr1.tar.gz')==False):
#     !wget https://zenodo.org/record/1421099/files/hess_dl3_dr1.tar.gz
#     !tar -zxvf hess_dl3_dr1.tar.gz

# In[62]:


src_name='Crab' #http://odahub.io/ontology#AstrophysicalObject
RA=83.628700  # http://odahub.io/ontology#PointOfInterestRA
DEC = 22.014700 # http://odahub.io/ontology#PointOfInterestDEC
T1='2003-10-09T13:16:00.0'# http://odahub.io/ontology#StartTime
T2='2005-10-10T13:16:00.0' # http://odahub.io/ontology#EndTime
Radius=2.5  #http://odahub.io/ontology#AngleDegrees
R_s=0.2     #http://odahub.io/ontology#AngleDegrees
Emin=100.    #http://odahub.io/ontology#Energy_GeV
Emax=10000. #http://odahub.io/ontology#Energy_GeV
NTbins=10 # http://odahub.io/ontology#Integer


# In[ ]:


with open('inputs.json', 'r') as fd:
    inp_dic = json.load(fd)
if '_data_product' in inp_dic.keys():
    inp_pdic = inp_dic['_data_product']
else:
    inp_pdic = inp_dic

for vn, vv in inp_pdic.items():
    if vn != '_selector':
        globals()[vn] = type(globals()[vn])(vv)


# In[63]:


T1=Time(T1, format='isot', scale='utc').mjd
T2=Time(T2, format='isot', scale='utc').mjd

dates1=data_store.obs_table['DATE-OBS']
dates2=data_store.obs_table['DATE-END']
times1=data_store.obs_table['TIME-OBS']
times2=data_store.obs_table['TIME-END']
OBSIDs=data_store.obs_table['OBS_ID']
Tstart=[]
Tstop=[]
for i in range(len(dates1)):
    Tstart.append(Time(dates1[i]+'T'+times1[i], format='isot', scale='utc').mjd)
    Tstop.append(Time(dates2[i]+'T'+times2[i], format='isot', scale='utc').mjd)
    
RA_pnts=np.array(data_store.obs_table['RA_PNT'])
DEC_pnts=np.array(data_store.obs_table['DEC_PNT'])
Tstart=np.array(Tstart)


# In[64]:


Coords_s=SkyCoord(RA,DEC,unit='degree')
COORDS_pnts=SkyCoord(RA_pnts,DEC_pnts,unit='degree')
seps=COORDS_pnts.separation(Coords_s).deg


# In[65]:


mask=np.where((seps<Radius) & (Tstart>T1) & (Tstop<T2))[0]
obs_ids=OBSIDs[mask]
Tbegs=Tstart[mask]
if(len(obs_ids)==0):
    message='No data found'
    raise RuntimeError('No data found')
obs_ids


# In[66]:


Tbins=np.linspace(T1,T2,NTbins+1)
Tmin=Tbins[:-1]
Tmax=Tbins[1:]
Tmean=(Tmin+Tmax)/2.
Tbins


# In[67]:


OBSlist=[]
for obs in obs_ids:
    OBSlist.append(hess_data+'/data/hess_dl3_dr1_obs_id_0'+str(obs)+'.fits.gz')
OBSlist


# In[68]:


flux=np.zeros(NTbins)
flux_err=np.zeros(NTbins)
flux_b=np.zeros(NTbins)
flux_b_err=np.zeros(NTbins)
Expos=np.zeros(NTbins)
for count,f in enumerate(OBSlist):
    hdul=fits.open(f)
    RA_pnt=hdul[1].header['RA_PNT']
    DEC_pnt=hdul[1].header['DEC_PNT']
    Texp=hdul[1].header['LIVETIME']
    Trun_start=hdul[1].header['TSTART']
    dRA=RA-RA_pnt
    dDEC=DEC-DEC_pnt
    RA_b=RA_pnt-dRA
    DEC_b=DEC_pnt-dDEC
    Coords_b=SkyCoord(RA_b,DEC_b,unit='degree')    
    Coords_pnt=SkyCoord(RA_pnt,DEC_pnt,unit='degree')
    dist=Coords_pnt.separation(Coords_s).deg

    ev=hdul['EVENTS'].data
    ev_ra=ev['RA']
    ev_dec=ev['DEC']
    ev_en=ev['ENERGY']
    ev_time=(ev['TIME']-Trun_start)/86400.+Tbegs[count]
    print(ev_time[0])
    ev_coords=SkyCoord(ev_ra,ev_dec,unit='degree')
    sep_s=ev_coords.separation(Coords_s).deg
    sep_b=ev_coords.separation(Coords_b).deg

    hdu=hdul['AEFF'].data
    EEmin=hdu['ENERG_LO'][0]
    EEmax=hdu['ENERG_HI'][0]
    EE=sqrt(EEmin*EEmax)
    EEbins=np.concatenate((EEmin,[EEmax[-1]]))
    AA=hdu['EFFAREA'][0]+1e-10
    Thmin=hdu['THETA_LO'][0]
    Thmax=hdu['THETA_HI'][0]   
    ind=np.argmin((Thmin-dist)**2)
    AA=AA[ind]*Texp*1e4
    mask=np.where((sep_s<R_s))
    cts1=np.histogram2d(ev_time[mask],ev_en[mask],bins=[Tbins,EEbins])[0]
    mask=(sep_b<R_s)
    cts2=np.histogram2d(ev_time[mask],ev_en[mask],bins=[Tbins,EEbins])[0]
    src=cts1-cts2
    src_err=sqrt(cts1+cts2)
    flux+=np.sum(src/AA,axis=1)
    flux_err+=np.sum(src_err/AA,axis=1)
    flux_b+=np.sum(cts2/AA,axis=1)
    flux_b_err+=np.sum(sqrt(cts2)/AA,axis=1)
    hdul.close()


# In[69]:


if(message==''):
    plt.errorbar(Tmean,flux,yerr=flux_err,xerr=[Tmean-Tmin,Tmax-Tmean],linestyle='none',label='source')
    plt.errorbar(Tmean,flux_b,yerr=flux_b_err,xerr=[Tmean-Tmin,Tmax-Tmean],linestyle='none',label='background')
    plt.xlabel('Time, MJD')
    plt.ylabel('Flux, cts/cm$^2$s')
    plt.yscale('log')
    ymin=min(min(flux-flux_err),min(flux_b-flux_b_err))
    ymax=max(max(flux+flux_err),max(flux_b+flux_b_err))    
    plt.ylim(ymin/2.,2*ymax)
    plt.legend(loc='lower left')
    plt.savefig('Lightcurve.png',format='png')


# In[70]:


if (message==''):
    bin_image = PictureProduct.from_file('Lightcurve.png')
    from astropy.table import Table
    data=[Tmean,Tmin,Tmax,flux,flux_err,flux_b,flux_b_err]
    names=('Tmean[MJD]','Tmin[MJD]','Tmax[MJD]','Flux[counts/cm2s]','Flux_error[counts/cm2s]','Background[counts/cm2s]','Background_error[counts/cm2s]')
    lc = ODAAstropyTable(Table(data, names = names))


# In[71]:


picture = bin_image # http://odahub.io/ontology#ODAPictureProduct
lightcurve_astropy_table = lc # http://odahub.io/ontology#ODAAstropyTable


# In[ ]:





# In[ ]:


_simple_outs, _oda_outs = [], []
_galaxy_meta_data = {}
_oda_outs.append(('out_Lightcurve_picture', 'picture_galaxy.output', picture))
_oda_outs.append(('out_Lightcurve_lightcurve_astropy_table', 'lightcurve_astropy_table_galaxy.output', lightcurve_astropy_table))

for _outn, _outfn, _outv in _oda_outs:
    _galaxy_outfile_name = os.path.join(_galaxy_wd, _outfn)
    if isinstance(_outv, str) and os.path.isfile(_outv):
        shutil.move(_outv, _galaxy_outfile_name)
        _galaxy_meta_data[_outn] = {'ext': '_sniff_'}
    elif getattr(_outv, "write_fits_file", None):
        _outv.write_fits_file(_galaxy_outfile_name)
        _galaxy_meta_data[_outn] = {'ext': 'fits'}
    elif getattr(_outv, "write_file", None):
        _outv.write_file(_galaxy_outfile_name)
        _galaxy_meta_data[_outn] = {'ext': '_sniff_'}
    else:
        with open(_galaxy_outfile_name, 'w') as fd:
            json.dump(_outv, fd, cls=CustomJSONEncoder)
        _galaxy_meta_data[_outn] = {'ext': 'json'}

for _outn, _outfn, _outv in _simple_outs:
    _galaxy_outfile_name = os.path.join(_galaxy_wd, _outfn)
    if isinstance(_outv, str) and os.path.isfile(_outv):
        shutil.move(_outv, _galaxy_outfile_name)
        _galaxy_meta_data[_outn] = {'ext': '_sniff_'}
    elif _numpy_available and isinstance(_outv, np.ndarray):
        with open(_galaxy_outfile_name, 'wb') as fd:
            np.savez(fd, _outv)
        _galaxy_meta_data[_outn] = {'ext': 'npz'}
    else:
        with open(_galaxy_outfile_name, 'w') as fd:
            json.dump(_outv, fd)
        _galaxy_meta_data[_outn] = {'ext': 'expression.json'}

with open(os.path.join(_galaxy_wd, 'galaxy.json'), 'w') as fd:
    json.dump(_galaxy_meta_data, fd)
print("*** Job finished successfully ***")

