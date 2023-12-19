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
from astropy.io import fits
import numpy as np
from numpy import pi,cos,sin,sqrt,log10
import os

from oda_api.data_products import PictureProduct
from oda_api.data_products import ODAAstropyTable
import os
from astropy.time import Time
if(os.path.exists('hess_dl3_dr1.tar.gz')==False):
    get_ipython().system('wget https://zenodo.org/record/1421099/files/hess_dl3_dr1.tar.gz')
    get_ipython().system('tar -zxvf hess_dl3_dr1.tar.gz')


# In[2]:


src_name='Crab' #http://odahub.io/ontology#AstrophysicalObject
RA = 83.628700  # http://odahub.io/ontology#PointOfInterestRA
DEC = 22.014700 # http://odahub.io/ontology#PointOfInterestDEC
T1='2000-10-09T13:16:00.0'# http://odahub.io/ontology#StartTime
T2='2022-10-10T13:16:00.0' # http://odahub.io/ontology#EndTime
Radius=2.5  #http://odahub.io/ontology#AngleDegrees
R_s=0.2     #http://odahub.io/ontology#AngleDegrees

Emin=100.    #http://odahub.io/ontology#Energy_GeV
Emax=10000. #http://odahub.io/ontology#Energy_GeV
NEbins=20 # http://odahub.io/ontology#Integer


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


# In[3]:


Emin=Emin/1e3
Emax=Emax/1e3
Ebins=np.logspace(log10(Emin),log10(Emax),NEbins+1)
Ebins
Emin=Ebins[:-1]
Emax=Ebins[1:]
Emean=sqrt(Emin*Emax)
lgEmean=log10(Emean)


# In[4]:


T1=Time(T1, format='isot', scale='utc').mjd
T2=Time(T2, format='isot', scale='utc').mjd
message=''
RA_pnts=[]
DEC_pnts=[]
DL3_files=[]
OBSIDs=[]
Tstart=[]
Tstop=[]
flist=os.listdir('data')
for f in flist:
    if(f[-7:]=='fits.gz'):
        DL3_files.append(f)
        OBSIDs.append(int(f[20:26]))
        hdul=fits.open('data/'+f)
        RA_pnts.append(float(hdul[1].header['RA_PNT']))
        DEC_pnts.append(float(hdul[1].header['DEC_PNT']))   
        Tstart.append(Time(hdul[1].header['DATE-OBS']+'T'+hdul[1].header['TIME-OBS'], format='isot', scale='utc').mjd)
        Tstop.append(Time(hdul[1].header['DATE-END']+'T'+hdul[1].header['TIME-END'], format='isot', scale='utc').mjd)
        hdul.close()


# In[5]:


Coords_s=SkyCoord(RA,DEC,unit='degree')
COORDS_pnts=SkyCoord(RA_pnts,DEC_pnts,unit='degree')
seps=COORDS_pnts.separation(Coords_s).deg


# In[6]:


mask=np.where((seps<Radius) & (Tstart>T1) & (Tstop<T2))[0]
OBSlist=[]
for i in mask:
    OBSlist.append(DL3_files[i])
if(len(OBSlist)==0):
    message='No data found'
    raise RuntimeError('No data found')
message


# In[9]:


cts_s=np.zeros(NEbins)
cts_b=np.zeros(NEbins)
Expos=np.zeros(NEbins)
for f in OBSlist:
    hdul=fits.open('data/'+f)
    RA_pnt=hdul[1].header['RA_PNT']
    DEC_pnt=hdul[1].header['DEC_PNT']
    Texp=hdul[1].header['LIVETIME']
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
    ev_time=ev['TIME'] 
    ev_coords=SkyCoord(ev_ra,ev_dec,unit='degree')
    sep_s=ev_coords.separation(Coords_s).deg
    sep_b=ev_coords.separation(Coords_b).deg

    hdu=hdul['AEFF'].data
    EEmin=hdu['ENERG_LO'][0]
    EEmax=hdu['ENERG_HI'][0]
    lgEE=log10(sqrt(EEmin*EEmax))
    lgAA=log10(hdu['EFFAREA'][0]+1e-10)
    Thmin=hdu['THETA_LO'][0]
    Thmax=hdu['THETA_HI'][0]   
    ind=np.argmin((Thmin-dist)**2)
    Expos+=10**(np.interp(lgEmean,lgEE,lgAA[ind]))*Texp
    mask=(sep_s<R_s)
    cts_s+=np.histogram(ev_en[mask],bins=Ebins)[0]
    mask=(sep_b<R_s)
    cts_b+=np.histogram(ev_en[mask],bins=Ebins)[0]
    hdul.close()


# In[10]:


flux=(cts_s-cts_b)/(Emax-Emin)*Emax*Emin/(Expos*1e4)
flux_err=sqrt(cts_s+cts_b)/(Emax-Emin)*Emax*Emin/(Expos*1e4)
plt.errorbar(Emean,flux,yerr=flux_err,xerr=[Emean-Emin,Emax-Emean])
plt.xscale('log')
plt.yscale('log')
plt.xlabel('$E$, TeV')
plt.ylabel('$E^2 dN/dE$, erg/(cm$^2$s)')
plt.savefig('Spectrum.png',format='png')


# In[11]:


bin_image = PictureProduct.from_file('Spectrum.png')
from astropy.table import Table
data=[Emean,Emin,Emax,flux,flux_err,cts_s,cts_b,Expos*1e4]
names=('Emean[TeV]','Emin[TeV]','Emax[TeV]','Flux[TeV/cm2s]','Flux_error[TeV/cm2s]','Cts_s','Cts_b','Exposure[cm2s]')
spec = ODAAstropyTable(Table(data, names = names))


# In[12]:


picture_png = bin_image # http://odahub.io/ontology#ODAPictureProduct
spectrum_astropy_table = spec # http://odahub.io/ontology#ODAAstropyTable


# In[ ]:





# In[ ]:


_simple_outs, _oda_outs = [], []
_galaxy_meta_data = {}
_oda_outs.append(('out_Spectrum_picture_png', 'picture_png_galaxy.output', picture_png))
_oda_outs.append(('out_Spectrum_spectrum_astropy_table', 'spectrum_astropy_table_galaxy.output', spectrum_astropy_table))

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

