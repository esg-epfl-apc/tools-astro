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
from numpy import pi,cos,sin,sqrt
import os
from astropy import wcs
from astropy.io import fits
from oda_api.data_products import PictureProduct
from oda_api.data_products import ImageDataProduct
from astropy.time import Time


# In[2]:


if(os.path.exists('hess_dl3_dr1.tar.gz')==False):
    get_ipython().system('wget https://zenodo.org/record/1421099/files/hess_dl3_dr1.tar.gz')
    get_ipython().system('tar -zxvf hess_dl3_dr1.tar.gz')


# In[3]:


src_name='Crab' #http://odahub.io/ontology#AstrophysicalObject
RA=83.628700  # http://odahub.io/ontology#PointOfInterestRA
DEC = 22.014700 # http://odahub.io/ontology#PointOfInterestDEC
T1='2000-10-09T13:16:00.0'# http://odahub.io/ontology#StartTime
T2='2022-10-10T13:16:00.0' # http://odahub.io/ontology#EndTime
Radius=2.5  #http://odahub.io/ontology#AngleDegrees
pixsize=0.1 #http://odahub.io/ontology#AngleDegrees ; oda:label "Pixel size"
Emin=100.    #http://odahub.io/ontology#Energy_GeV
Emax=10000. #http://odahub.io/ontology#Energy_GeV


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


# In[8]:


cdec=cos(DEC*pi/180.)
Npix=int(4*Radius/pixsize)+1
RA_bins=np.linspace(RA-Radius/cdec,RA+Radius/cdec,Npix+1)
DEC_bins=np.linspace(DEC-Radius,DEC+Radius,Npix+1)
image=np.zeros((Npix,Npix))
for f in OBSlist:
    hdul=fits.open('data/'+f)
    ev=hdul['EVENTS'].data
    ev_ra=ev['RA']
    ev_dec=ev['DEC']
    ev_en=ev['ENERGY']
    ev_time=ev['TIME']    
    h=np.histogram2d(ev_ra,ev_dec,bins=[RA_bins,DEC_bins])
    image+=h[0]
    hdul.close()

plt.imshow(np.flip(image,axis=1),extent=(RA_bins[-1],RA_bins[0],DEC_bins[0],DEC_bins[-1]),origin='lower')
plt.colorbar()

plt.xlabel('RA, degrees')
plt.ylabel('DEC,degrees')
plt.savefig('Image.png',format='png')


# In[9]:


# Create a new WCS object.  The number of axes must be set
# from the start
w = wcs.WCS(naxis=2)

# Set up an "Airy's zenithal" projection
# Vector properties may be set with Python lists, or Numpy arrays
w.wcs.crpix = [Npix/2., Npix/2.]
w.wcs.cdelt = np.array([pixsize/cdec, pixsize])
w.wcs.crval = [RA, DEC]
w.wcs.ctype = ["RA---AIR", "DEC--AIR"]
w.wcs.set_pv([(2, 1, 45.0)])

# Now, write out the WCS object as a FITS header
header = w.to_header()

# header is an astropy.io.fits.Header object.  We can use it to create a new
# PrimaryHDU and write it to a file.
hdu = fits.PrimaryHDU(image,header=header)
hdu.writeto('Image.fits',overwrite=True)
hdu=fits.open('Image.fits')
im=hdu[0].data
from astropy.wcs import WCS
wcs = WCS(hdu[0].header)
plt.subplot(projection=wcs)
plt.imshow(im,  origin='lower')
plt.grid(color='white', ls='solid')
plt.xlabel('RA')
plt.ylabel('Dec')


# In[10]:


bin_image = PictureProduct.from_file('Image.png')
fits_image=ImageDataProduct.from_fits_file('Image.fits')


# In[11]:


picture = bin_image # http://odahub.io/ontology#ODAPictureProduct
image = fits_image # http://odahub.io/ontology#Image


# In[ ]:





# In[ ]:


_simple_outs, _oda_outs = [], []
_galaxy_meta_data = {}
_oda_outs.append(('out_Image_picture', 'picture_galaxy.output', picture))
_oda_outs.append(('out_Image_image', 'image_galaxy.output', image))

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

