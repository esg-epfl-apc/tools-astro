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


import numpy as np
import matplotlib.pyplot as plt
from numpy import log10,pi,log,exp
from scipy.interpolate import interp1d
import glob

from oda_api.data_products import PictureProduct
from oda_api.data_products import ODAAstropyTable


# In[2]:


get_ipython().system('pwd')


# In[3]:


get_ipython().system('git lfs install ; git lfs pull')


# In[4]:


get_ipython().run_cell_magic('bash', '', './install.sh\n')


# In[15]:


Gamma=2.4 # http://odahub.io/ontology#Float
Emax=1e15 # http://odahub.io/ontology#Float
Ecut=1e14 # http://odahub.io/ontology#Float
B=1e4     # http://odahub.io/ontology#Float
source_size_cm=3.0856e+13  # http://odahub.io/ontology#Float
background_norm_mode='density_cm3' # http://odahub.io/ontology#String ; oda:allowed_value "absolute","density_cm3","energy_density_eV_cm3"
corona_backround_norm_cm3=4.18876e15  # http://odahub.io/ontology#Float
disk_background_norm_cm3=1.39654e17   # http://odahub.io/ontology#Float
min_steps=100 #http://odahub.io/ontology#Integer ; oda:lower_limit 10 ; oda:upper_limit 10000


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


# In[16]:


background_norm_mode_index=["absolute","density_cm3","energy_density_eV_cm3"].index(background_norm_mode)
distance_Mpc = source_size_cm/3.0856775807e+24
electron_distance_Mpc = 100*distance_Mpc # electrons propagate longer
min_step_Mpc = distance_Mpc/min_steps
electron_min_step_Mpc = electron_distance_Mpc/min_steps


# In[17]:


fname='pp_a'+str(Gamma)+'_'+str(f"{Ecut:.2e}")+'_'+str(B)
fname='run'
get_ipython().system('cp propagation/pp_a2.4_E1e14.xsw propagation/{fname}.xsw')


# In[18]:


get_ipython().system('python propagation/replacePar CommonAlpha {str(Gamma)} propagation/{fname}.xsw')
get_ipython().system('python propagation/getPar CommonAlpha propagation/{fname}.xsw')

get_ipython().system('python propagation/replacePar Emax {str(Emax)} propagation/{fname}.xsw')
get_ipython().system('python propagation/getPar Emax propagation/{fname}.xsw')

Ecut_rel=Ecut/Emax
get_ipython().system('python propagation/replacePar injSpectraHigherCutoff {str(Ecut_rel)} propagation/{fname}.xsw')
get_ipython().system('python propagation/getPar injSpectraHigherCutoff propagation/{fname}.xsw')


get_ipython().system('python propagation/replacePar B0 {str(B)} propagation/step??_template.xsw')
get_ipython().system('python propagation/getPar B0 propagation/step2g_template.xsw')

get_ipython().system('python propagation/replacePar iroMult {str(disk_background_norm_cm3)} propagation/step??_template.xsw')
get_ipython().system('python propagation/getPar iroMult propagation/step2g_template.xsw')

get_ipython().system('python propagation/replacePar IROBackgroundNormMode {str(background_norm_mode_index)} propagation/step??_template.xsw')
get_ipython().system('python propagation/getPar IROBackgroundNormMode propagation/step2g_template.xsw')

get_ipython().system('python propagation/replacePar CustomBackgroundNormMode {str(background_norm_mode_index)} propagation/step??_template.xsw')
get_ipython().system('python propagation/getPar CustomBackgroundNormMode propagation/step2g_template.xsw')

get_ipython().system('python propagation/replacePar L_short_distance_test {str(distance_Mpc)} propagation/step?g_template.xsw')
get_ipython().system('python propagation/getPar L_short_distance_test propagation/step2g_template.xsw')

get_ipython().system('python propagation/replacePar L_short_distance_test {str(electron_distance_Mpc)} propagation/step3e_template.xsw')
get_ipython().system('python propagation/getPar L_short_distance_test propagation/step3e_template.xsw')

get_ipython().system('python propagation/replacePar microStep {str(min_step_Mpc)} propagation/step?g_template.xsw')
get_ipython().system('python propagation/getPar microStep propagation/step2g_template.xsw')

get_ipython().system('python propagation/replacePar microStep {str(electron_min_step_Mpc)} propagation/step3e_template.xsw')
get_ipython().system('python propagation/getPar microStep propagation/step3e_template.xsw')


# In[9]:


get_ipython().run_cell_magic('bash', '', 'cd propagation\n./run_all.sh run.xsw\n')


# In[10]:


get_ipython().system("cat 'propagation/results/run_step4g/fin/tot'")


# In[11]:


d=np.genfromtxt('propagation/results/run_step4g/fin/tot')
E=d[:,0]
gam=d[:,1]
nu=d[:,2]
elec=d[:,3]
prot=d[:,4]
plt.plot(E,gam)
plt.plot(E,nu)
plt.xscale('log')
plt.yscale('log')
plt.ylim(1e-20,1e-8)


# In[12]:


d=np.genfromtxt('NGC_1068_contour.csv')
gam=d[:,0]
f=d[:,1]
x_nu=np.logspace(3,5,10)
ymax=np.zeros(len(x_nu))
ymin=np.ones(len(x_nu))
for i in range(len(gam)):
    y=3*f[i]*1e-11*(x_nu/1e3)**(2-gam[i])
    ymax=np.maximum(y,ymax)
    ymin=np.minimum(y,ymin)
plt.fill_between(x_nu,ymin,ymax,color='blue',alpha=0.5)

d=np.genfromtxt('SED_1068.csv')
nu=d[:,0]
nufnu=d[:,1]
en=2*pi*6.6e-16*nu/1e9
efe=nufnu*1e-23/1.6
plt.scatter(en,efe)

d=np.genfromtxt('Fermi_1068.csv')
ee=d[:,0]
ff=d[:,1]
ff_max=d[:,2]
ff_min=d[:,3]
plt.errorbar(ee,ff,yerr=[ff-ff_min,ff_max-ff])

A=7e43/1.6e-12 /1e6 /log(10/2.) #1/eV/s
FX0=A*1e6/1e12/4/pi/(16.3*3e24)**2
print(FX0)
x=[1e-6,1e-4]
y=[FX0,FX0]
plt.plot(x,y,color='black',linestyle='dashed')
A=4e43/1.6e-12 /1e6 /log(10/2.) #1/eV/s
FX0=A*1e6/1e12/4/pi/(16.3*3e24)**2
A=14e43/1.6e-12 /1e6 /log(10/2.) #1/eV/s
FX1=A*1e6/1e12/4/pi/(16.3*3e24)**2
x=[1e-6,1e-4,1e-4,1e-6]
y=[FX1,FX1,FX0,FX0]
plt.fill(x,y,alpha=0.3,color='black',linestyle='none')

spec_file = 'propagation/results/run_step4g/fin/tot'
d=np.genfromtxt(spec_file)
E=d[:,0]*1e-9
gam=d[:,1]     # for original format d[:,3]
nu=d[:,2]           # for original format (d[:,4]+d[:,5]+d[:,7]+d[:,8])
f_nu = interp1d(np.log(E), nu , kind='linear', fill_value=0.)
interp_nu = f_nu(np.log(x_nu))
N = 0.8*np.min(ymax/interp_nu)
plt.plot(E,N*gam*(E>5.e-9), label='$\\gamma$',color='blue',alpha=0.3)
plt.plot(E,N*gam*(E>5.e-7), color='black')
plt.plot(E,N*nu, label='$\\nu$')
backgr = np.genfromtxt('propagation/results/run_step4g/fin/backgr')
Eb = backgr[:,0] * 1e-9
FE2 = backgr[:,0] * backgr[:,1] * 1e-12 * 3e10/4/np.pi*(1e13/4.6285163711e+25)**2
plt.plot(Eb, 20*FE2, label='backgr3')

plt.xscale('log')
plt.yscale('log')
plt.ylim(6e-14,3e-8)
plt.xlim(3e-11,1e6)
plt.tick_params(labelsize=12)
plt.xlabel(r'$E$, GeV',fontsize=12)
plt.ylabel(r'$EF_E$, TeV/(cm$^2$s)',fontsize=12)
plt.savefig('Spectrum.png',format='png',bbox_inches='tight')


# In[13]:


bin_image = PictureProduct.from_file('Spectrum.png')
from astropy.table import Table
data=[E,gam,nu]
names=('E[eV]','E2 dN/dE gamma [TeV/cm2s]','E2 dN/dE nu[ TeV/cm2 s]')
spectrum = ODAAstropyTable(Table(data, names = names))


# In[14]:


picture = bin_image # http://odahub.io/ontology#ODAPictureProduct
spectrum_astropy_table = spectrum # http://odahub.io/ontology#ODAAstropyTable


# In[ ]:





# In[ ]:


_simple_outs, _oda_outs = [], []
_galaxy_meta_data = {}
_oda_outs.append(('out_Simulate_picture', 'picture_galaxy.output', picture))
_oda_outs.append(('out_Simulate_spectrum_astropy_table', 'spectrum_astropy_table_galaxy.output', spectrum_astropy_table))

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
