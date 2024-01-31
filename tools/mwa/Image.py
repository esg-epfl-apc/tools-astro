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


from astropy.io import fits
import astropy.units as u
import matplotlib.pyplot as plt
from astropy.wcs import WCS
from astroquery.skyview import SkyView
from oda_api.data_products import PictureProduct
from oda_api.data_products import ImageDataProduct


# In[2]:


src_name='1ES 0229+200' #http://odahub.io/ontology#AstrophysicalObject
RA=38.202562  # http://odahub.io/ontology#PointOfInterestRA
DEC =20.288191 # http://odahub.io/ontology#PointOfInterestDEC
T1='2000-10-09T13:16:00.0'# http://odahub.io/ontology#StartTime
T2='2022-10-10T13:16:00.0' # http://odahub.io/ontology#EndTime
Radius=1.  #http://odahub.io/ontology#AngleDegrees
pixsize=0.01 #http://odahub.io/ontology#AngleDegrees
Frequency='GLEAM 170-231 MHz'    # http://odahub.io/ontology#String ; oda:allowed_value "GLEAM 72-103 MHz","GLEAM 103-134 MHz","GLEAM 139-170 MHz","GLEAM 170-231 MHz"


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


pixels=int(2*Radius/pixsize)+1
Radius*=u.deg
pos=str(RA)+', '+str(DEC)
pixels


# In[4]:


hdul=SkyView.get_images(position=pos,
                       survey=[Frequency],pixels=pixels,radius=Radius)  


# In[5]:


hdu=hdul[0]
hdu[0].header
wcs = WCS(hdu[0].header)


# In[6]:


image=hdu[0].data


# In[7]:


ax=plt.subplot(projection=wcs)
im=ax.imshow(image,origin='lower')
ax.coords.grid(True, color='white', ls='solid')
plt.colorbar(im,label='Jy/beam')
plt.savefig('Image.png',format='png',bbox_inches='tight')


# In[8]:


hdu.writeto('Image.fits',overwrite=True)
bin_image = PictureProduct.from_file('Image.png')
fits_image=ImageDataProduct.from_fits_file('Image.fits')


# In[9]:


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

