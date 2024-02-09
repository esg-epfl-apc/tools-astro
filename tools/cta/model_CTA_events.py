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


from pathlib import Path
import numpy as np
import astropy.units as u
from astropy.coordinates import Angle, SkyCoord
from astropy.io import fits
from astropy.time import Time
from regions import CircleSkyRegion
import matplotlib.pyplot as plt
from IPython.display import display
from gammapy.data import (
    DataStore,
    FixedPointingInfo,
    Observation,
    observatory_locations,
)
from gammapy.datasets import MapDataset, MapDatasetEventSampler
from gammapy.irf import load_irf_dict_from_file
from gammapy.makers import MapDatasetMaker
from gammapy.maps import MapAxis, WcsGeom, Map, RegionGeom, WcsMap
from gammapy.modeling.models import (
    ExpDecayTemporalModel,
    FoVBackgroundModel,
    Models,
    PointSpatialModel,
    PowerLawNormSpectralModel,
    PowerLawSpectralModel,
    SkyModel,
    TemplateSpatialModel,
)

from gammapy.data import PointingMode
from gammapy.irf import load_cta_irfs

from oda_api.data_products import PictureProduct
from oda_api.data_products import LightCurveDataProduct
from oda_api.data_products import ODAAstropyTable
from oda_api.data_products import NumpyDataProduct
from astropy.table import Table 

from astropy.nddata import StdDevUncertainty
from astropy import units as u

from numpy import log10,cos,pi,sqrt
from numpy.random import randint

from numpy import log,exp


# In[20]:


get_ipython().run_cell_magic('bash', '', 'rm -r IRFS | echo "Ok"\nmkdir IRFS\ncd IRFS\nwget https://zenodo.org/records/5499840/files/cta-prod5-zenodo-fitsonly-v0.1.zip\nunzip cta-prod5-zenodo-fitsonly-v0.1.zip\ncd fits\nfor fn in *.gz ; do tar -zxvf $fn; done \n')


# In[18]:


#We simulate point source in wobble observaiton, 
#0.4 degree off-axis

src_name='Mrk 421' #http://odahub.io/ontology#AstrophysicalObject
RA=166.113809  # http://odahub.io/ontology#PointOfInterestRA
DEC =38.208833 # http://odahub.io/ontology#PointOfInterestDEC
T1='2000-10-09T13:16:00.0'# http://odahub.io/ontology#StartTime
T2='2022-10-10T13:16:00.0' # http://odahub.io/ontology#EndTime

#Exposure time in hours
Texp=50. #http://odahub.io/ontology#TimeIntervalHours

#Source redshift
z=0.1 #http://odahub.io/ontology#Float

#Source flux normalisaiton F0 in 1/(TeV cm2 s) at reference energy E0 
F0=4e-13 # http://odahub.io/ontology#Float
E0=1. # http://odahub.io/ontology#Energy_TeV
Gamma=1.75 #http://odahub.io/ontology#Float

#source extension in degrees
sigma=0. #http://odahub.io/ontology#Float

#Pointing RA and DEC
#pnt_RA=166.613809   # http://odahub.io/ontology#RightAscensionDegrees
#pnt_DEC = 38.208833 # http://odahub.io/ontology#DeclinationDegrees


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


# In[19]:


#We use "cube" model filled with counts
#this is experimentally established version of 
#conversion of the flux normalisation to model intrinsic units
norm_ps=4.67*F0/1e-10

# For point sources, we will simulate in a square of 1 degree around the source
# with pixels 0.02 degree
ROI=1.
pixsize=0.02
if (sigma<0.02):
    sigma=0.02

source=SkyCoord(RA,DEC,unit='degree')
cdec=cos(DEC*pi/180.)
pnt_RA=RA-0.4/cdec
pnt_DEC=DEC
pnt=SkyCoord(pnt_RA,pnt_DEC,unit='degree')


# In[22]:


#Here we read the IRFs
hdul=fits.open('IRFS/fits/Prod5-North-20deg-AverageAz-4LSTs09MSTs.18000s-v0.1.fits.gz')
Aeff=hdul['EFFECTIVE AREA'].data
th_min=Aeff['THETA_LO']
th_min
E_MIN=Aeff['ENERG_LO'][0][:-2]
E_MAX=Aeff['ENERG_HI'][0][:-2]
E=sqrt(E_MIN*E_MAX)
Ebins=np.concatenate((E_MIN,[E_MAX[-1]]))
n_dec=log10(E_MAX)-log10(E_MIN)
n_Ebins=len(E_MIN)
A=Aeff['EFFAREA'][0,0][:-2]


# In[23]:


d=np.genfromtxt('Franceschini17.txt')
en=d[:,0]
tau_001=d[:,1]
tau_003=d[:,2]
tau_01=d[:,3]
tau_03=d[:,4]
tau_05=d[:,5]
tau_10=d[:,6]
tau_15=d[:,7]
tau_20=d[:,8]
tau_30=d[:,9]

zgrid=np.array([0.01,0.03,0.1,0.3,0.5,1.0,1.5,2.0,3.0])
mask=(z-zgrid)>0
ind=len(zgrid[mask])-1
print(ind)
increment=(z-zgrid[ind])/(zgrid[ind+1]-zgrid[ind])
tau=d[:,1:]


tau_ss=(tau[:,ind+1]-tau[:,ind])*increment+tau[:,ind]
tau_s=np.interp(log(E),log(en),tau_ss)
plt.scatter(E,tau_s)
#plt.plot(E,tau_s)
plt.scatter(en,tau_ss)
if z<=0.01:
    tau_s=0.*tau_s
    
def tau(x):
    return np.interp(log(x),log(en),tau_ss)


# In[24]:


#Just in case, point source is a narrow Gaussian
def Gauss(x, y, A, x0, y0, sigma_x, sigma_y):
    return A*np.exp(-(x-x0)**2/(2*sigma_x**2) -(y-y0)**2/(2*sigma_y**2))


# In[25]:


cdec=cos(DEC*pi/180)
Npix=int(2*ROI/pixsize)+1
src_pixx=int(Npix/2)
print(src_pixx)
print(Npix)

x = np.linspace(RA-ROI/cdec,RA+ROI/cdec,Npix)
y = np.linspace(DEC-ROI,DEC+ROI,Npix)
x0 = RA
y0 = DEC
A = 1
Xg, Yg = np.meshgrid(x, y)

# point source will be a Gaussian of o.02 degrees
Z_ps = Gauss(Xg, Yg, A, x0, y0, sigma/cdec, sigma)

cube=np.zeros((n_Ebins,Npix,Npix))
for i in range(n_Ebins):
    cube[i]=norm_ps*Z_ps*(E[i])**(-(Gamma-1))*exp(-tau_s[i])
    


# In[26]:


hdul=fits.open('3d.fits',mode='update')
hdul[0].header['CRVAL1']=RA
hdul[0].header['CRVAL2']=DEC
hdul[0].header['CDELT2']=pixsize
hdul[0].header['CDELT1']=-pixsize/cdec
hdul[0].header['NAXIS1']=Npix
hdul[0].header['NAXIS2']=Npix
hdul[0].header['CRPIX1']=Npix/2.
hdul[0].header['CRPIX2']=Npix/2.
hdul[0].header['NAXIS3']=n_Ebins
hdul[0].header['WCSSHAPE']=str((Npix,Npix,n_Ebins))
hdul[1].data['E_MIN']=E_MIN
hdul[1].data['E_MAX']=E_MAX
hdul[1].data['ENERGY']=E
hdul[0].data=cube
hdul.close()


# In[27]:


m_cube=Map.read('3d.fits')
hdul=fits.open('3d.fits')
cube=hdul[0].data
ebins=hdul[1].data
E_MIN=ebins['E_MIN']
E_MAX=ebins['E_MAX']
E_MEAN=ebins['ENERGY']
n_bins_per_decade=int(len(E_MEAN)/(log10(E_MAX[-1])-log10(E_MIN[0])))
dth=hdul[0].header['CDELT2']
Npix_ph=hdul[0].header['NAXIS1']
Npix_th=hdul[0].header['NAXIS2']
Nebins=hdul[0].header['NAXIS3']
src_ph=hdul[0].header['CRPIX1']
src_th=hdul[0].header['CRPIX2']
src_ra=hdul[0].header['CRVAL1']
src_dec=hdul[0].header['CRVAL2']
source=SkyCoord(src_ra,src_dec,unit='degree')
m_cube


# In[28]:


# telescope is pointing at a fixed position in ICRS for the observation
pointing = FixedPointingInfo(
    fixed_icrs=pnt, mode=PointingMode.POINTING
)
livetime = Texp * u.hr
location = observatory_locations["cta_north"]

# irfs = load_irf_dict_from_file(path / irf_filename)
filename = "data/Prod5-North-20deg-AverageAz-4LSTs09MSTs.180000s-v0.1.fits.gz"
#irfs = load_cta_irfs(filename)
irfs = load_irf_dict_from_file(filename)
irfs


# In[29]:


map0 = m_cube.slice_by_idx({"energy": 0})
map0.plot()
plt.show()


# In[30]:


energy_axis_true = MapAxis.from_energy_bounds(
    str(E_MIN[0])+" TeV", str(E_MAX[-1])+" TeV", nbin=n_bins_per_decade, per_decade=True, name="energy_true"
)

#energy_axis = MapAxis.from_energy_bounds("1 TeV", "10 TeV", nbin=selected_n_bins_per_decade, per_decade=True)
energy_axis = MapAxis.from_energy_bounds(
    str(E_MIN[0])+" TeV", str(E_MAX[-1])+" TeV", nbin=n_bins_per_decade, per_decade=True
)

geom = WcsGeom.create(
    skydir=pointing.fixed_icrs,
    width=(2*ROI, 2*ROI),
    binsz=pixsize,
    frame="icrs",
    axes=[energy_axis],
)

migra_axis = MapAxis.from_bounds(-2, 2, nbin=150, node_type="edges", name="migra")


# In[31]:


def GetBinSpectralModel(E, bins_per_decade=20, norm=1):
    amplitude=1e-12 * u.Unit("cm-2 s-1") * norm
    from gammapy.modeling.models import (
        GaussianSpectralModel,
    )
    sigma = (10**(1/bins_per_decade)-1) * E
    return GaussianSpectralModel(mean=E, sigma=sigma, amplitude=amplitude)


# In[32]:


spec = m_cube.get_spectrum()
spec.plot()


# In[33]:


energy_bins = m_cube.geom.axes['energy'].center
len(energy_bins), float(np.min(energy_bins)/u.TeV), float(np.max(energy_bins)/u.TeV)


# In[34]:


int_bin_flux = spec.data.flatten()


# In[35]:


observation = Observation.create(
    obs_id=1001,
    pointing=pointing,
    livetime=livetime,
    irfs=irfs,
    location=location,
)
print(observation)


# In[36]:


empty = MapDataset.create(
    geom,
    energy_axis_true=energy_axis_true,
    migra_axis=migra_axis,
    name="my-dataset",
)
maker = MapDatasetMaker(selection=["exposure", "background", "psf", "edisp"])
dataset = maker.run(empty, observation)


# In[37]:


bin_models = []
for i, (flux, E) in enumerate(zip(int_bin_flux, energy_bins)):
    if flux == 0:
        continue
    spectral_model_delta = GetBinSpectralModel(E, norm=flux) # normalizing here
    spacial_template_model = TemplateSpatialModel(m_cube.slice_by_idx({"energy": i}), filename=f'cube_bin{i}.fit', normalize=True)
    sky_bin_model = SkyModel(
        spectral_model=spectral_model_delta,
        spatial_model=spacial_template_model,
        name=f"bin_{i}",
    )
    bin_models.append(sky_bin_model)


# In[38]:


bkg_model = FoVBackgroundModel(dataset_name="my-dataset")
models = Models(bin_models + [bkg_model])


# In[39]:


dataset.models = models


# In[ ]:


sampler = MapDatasetEventSampler(random_state=0)
events = sampler.run(dataset, observation)


# In[ ]:


E=events.energy/u.TeV
ras=events.radec.ra.deg
decs=events.radec.dec.deg
#plt.hist(E,bins=np.logspace(-2,2,41))

mask=events.table['MC_ID'] > 0
plt.hist(E[mask],bins=np.logspace(-2,2,41),alpha=0.5,label='source')
mask=events.table['MC_ID'] == 0
plt.hist(E[mask],bins=np.logspace(-2,2,41),alpha=0.5,label='background')


plt.xscale('log')
plt.yscale('log')
plt.legend(loc='upper right')
plt.savefig('event_spectrum.png',format='png')


# In[ ]:


mask=np.where(E>1)
print(len(E[mask]),len(ras[mask]),len(decs[mask]))

plt.hist2d(ras[mask],decs[mask],bins=[np.linspace(pnt_RA-ROI/cdec,pnt_RA+ROI/cdec,Npix),np.linspace(DEC-ROI,DEC+ROI,Npix)])
plt.colorbar()


# In[ ]:


print(f"Save events ...") 
primary_hdu = fits.PrimaryHDU()
hdu_evt = fits.BinTableHDU(events.table)
hdu_gti = fits.BinTableHDU(dataset.gti.table, name="GTI")
hdu_all = fits.HDUList([primary_hdu, hdu_evt, hdu_gti])
hdu_all.writeto(f"./events.fits", overwrite=True)
####################


# In[ ]:


hdul=fits.open('events.fits')
T_exp=hdul['EVENTS'].header['ONTIME']
events=hdul['EVENTS'].data


# In[ ]:


coords_s=SkyCoord(RA,DEC,unit='degree')
RA_bkg=pnt_RA-(RA-pnt_RA)
DEC_bkg=pnt_DEC-(DEC-pnt_DEC)
coords_b=SkyCoord(RA_bkg,DEC_bkg,unit='degree')

Es=events['ENERGY']
ras=events['RA']
decs=events['DEC']
coords=SkyCoord(ras,decs,unit='degree')
seps_s=coords.separation(coords_s).deg
seps_b=coords.separation(coords_b).deg
coords_s=SkyCoord(RA,DEC,unit='degree')


# In[ ]:


hdul=fits.open('Prod5-North-20deg-AverageAz-4LSTs09MSTs.180000s-v0.1.fits.gz')
Aeff=hdul['EFFECTIVE AREA'].data
th_min=Aeff['THETA_LO']
th_min
Emin_irf=Aeff['ENERG_LO'][0]
Emax_irf=Aeff['ENERG_HI'][0]
E_irf=sqrt(Emin_irf*Emax_irf)
Emin_irf
Ebins_irf=np.concatenate((Emin_irf,[Emax_irf[-1]]))
A=Aeff['EFFAREA'][0,0]


# In[ ]:


th_cut=0.3
mask=(seps_s<th_cut)
E_s=Es[mask]
h1=plt.hist(E_s,bins=Ebins_irf,alpha=0.5)
mask=(seps_b<th_cut)
E_b=Es[mask]
h2=plt.hist(E_b,bins=Ebins_irf,alpha=0.5)
plt.xscale('log')
plt.yscale('log')
Ns=h1[0]
Nb=h2[0]
Src=Ns-Nb
Src_err=sqrt(Ns+Nb)
print(Src,Src_err)


# In[ ]:


EE=E_irf[4:-3]
Flux=Src/(Emax_irf-Emin_irf)*E_irf**2/(A*1e4)/T_exp
Flux_err=Src_err/(Emax_irf-Emin_irf)*E_irf**2/(A*1e4)/T_exp
Flux_err=Flux_err+100*(Flux_err==0)
FFlux=Flux[4:-3]
FFlux_err=Flux_err[4:-3]

from scipy.optimize import curve_fit
def PL(x,Norm):
    return Norm*(x/E0)**(2-Gamma)*exp(-tau(x))


plt.errorbar(EE,FFlux,yerr=FFlux_err)

popt, pcov = curve_fit(PL, EE, FFlux, sigma=FFlux_err)
F0_best=popt[0]
y=F0_best*(EE/E0)**(2-Gamma)*exp(-tau(EE))
plt.yscale('log')
plt.xscale('log')
F0_err=sqrt(pcov[0,0])
SN=F0_best/F0_err
plt.plot(EE,y,label='S/N='+str(SN))
plt.legend(loc='upper right')
plt.ylim(1e-14,1e-9)
plt.savefig('spectrum.png',format='png')


# In[ ]:


#theta2 plot
thbin=0.1
nb=int(1/thbin**2)
print(nb)
bins=np.linspace(0,1,nb)
h1=plt.hist(seps_s**2,bins=bins,alpha=0.5)
h2=plt.hist(seps_b**2,bins=bins,alpha=0.5)
plt.axvline(th_cut**2)
plt.xlim(0,0.3)

cts_s=sum(h1[0][:10])
cts_b=sum(h2[0][:10])
SN=(cts_s-cts_b)/sqrt(cts_b)
print((cts_s-cts_b)/sqrt(cts_b))
plt.text(0.15,max(h1[0][:10]),'S/N='+str(SN))

plt.savefig('theta2.png',format='png')


# In[ ]:


events_fits = NumpyDataProduct.from_fits_file('events.fits')
spectrum_png = PictureProduct.from_file('spectrum.png')
theta2_png = PictureProduct.from_file('theta2.png')


# In[ ]:


theta2 = theta2_png # http://odahub.io/ontology#ODAPictureProduct
spectrum = spectrum_png # http://odahub.io/ontology#ODAPictureProduct
events_fits = events_fits # https://odahub.io/ontology/#Spectrum


# In[ ]:





# In[ ]:





# In[ ]:


_simple_outs, _oda_outs = [], []
_galaxy_meta_data = {}
_oda_outs.append(('out_model_CTA_events_theta2', 'theta2_galaxy.output', theta2))
_oda_outs.append(('out_model_CTA_events_spectrum', 'spectrum_galaxy.output', spectrum))
_simple_outs.append(('out_model_CTA_events_events_fits', 'events_fits_galaxy.output', events_fits))

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

