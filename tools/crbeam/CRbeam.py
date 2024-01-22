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


# In[ ]:


src_name='NGC 1365' #http://odahub.io/ontology#AstrophysicalObject

z_start=0  #http://odahub.io/ontology#Float
Npart=2000 #http://odahub.io/ontology#Integer ; oda:lower_limit 1 ; oda:upper_limit 100000
particle_type='gamma' # http://odahub.io/ontology#String ; oda:allowed_value "gamma","electron","proton"
Emax=30   #http://odahub.io/ontology#Energy_TeV
Emin=0.01   #http://odahub.io/ontology#Energy_TeV
EminSource=1.   #http://odahub.io/ontology#Energy_TeV
Gamma=2. #http://odahub.io/ontology#Float
EGMF_fG=10 #http://odahub.io/ontology#Float
lmaxEGMF_Mpc=5 #http://odahub.io/ontology#Float
jet_half_size=5.  #http://odahub.io/ontology#degree
jet_direction=0.  #http://odahub.io/ontology#degree
psf=1.  #http://odahub.io/ontology#degree
EBL='Franceschini 2017' # http://odahub.io/ontology#String ; oda:allowed_value "Franceschini 2017","Stecker 2016 lower limit","Stecker 2016 upper limit","Inoue 2012 Baseline","Inoue 2012 lower limit","Inoue 2012 upper limit"


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


# In[ ]:


from oda_api.data_products import PictureProduct
from oda_api.data_products import LightCurveDataProduct
from oda_api.data_products import ODAAstropyTable
from oda_api.data_products import NumpyDataProduct
from astropy.table import Table

from specutils import Spectrum1D
from astropy.nddata import StdDevUncertainty
from astropy import units as u
from utils import find_redshift


# In[ ]:


if z_start <= 0:
    z_start = find_redshift(src_name)


# In[ ]:


from crbeam import CRbeam
import numpy as np
from matplotlib import pyplot as plt
from spec_plot import plot_spec
import subprocess
import light_curve as lc

EGMF=EGMF_fG*1e-15

# internal units are eV
Emax *= 1e12
Emin *= 1e12
EminSource *= 1e12

background = {
    "Franceschini 2017" : 12,
    "Stecker 2016 lower limit" : 10,
    "Stecker 2016 upper limit" : 11,
    "Inoue 2012 Baseline" : 3,
    "Inoue 2012 lower limit" : 4,
    "Inoue 2012 upper limit": 5
    }[EBL]

prog = CRbeam(z=z_start, nparticles=Npart, primary=particle_type, emax=Emax, emin=Emin, emin_source=EminSource, EGMF=EGMF, lmaxEGMF=lmaxEGMF_Mpc, background=background)
cmd = prog.command
cmd


# In[ ]:


# Here is one way to call CRbeam without global cache
# prog.run(force_overwrite=False)
# Here we call CRbeam with global cache
# prog.run_cached(overwrite_local_cache=True)
# to clear s3 cache run the following command
# prog.remove_cache()


# In[ ]:


# To report the progress we will split running CRbeam into 10 parts
n_steps=10
# Initialize multistep simulation
data_exists = not prog.start_multistage_run(overwrite_local_cache=True, n_steps=n_steps)
proceed = not data_exists


# In[ ]:


prog.p.nparticles


# In[ ]:


proceed


# In[ ]:


prog.output_dir
# next we run every step in a separate cell


# In[ ]:


if proceed:
    for step in range(n_steps):
        print(f"running simulation {100 * step // n_steps}%", prog.output_dir)
        proceed = prog.simulation_step()
        # todo: report progress using rest API
    print(f"running simulation 100%", prog.output_dir)


# In[ ]:


assert not proceed, "must be completed before this cell"
if not data_exists:
    prog.end_multistep_run()


# In[ ]:


def adjust_weights(mc_file, power):
# converting weights to mimic required injection spectrum power
    header = ''
    with open(mc_file, 'rt') as lines:
        for line in lines:
            if len(line) > 0 and line[0] == '#':
                header += line[1:].strip() + '\n'
            else:
                break
    weight_col = 2
    E_src_col = 12
    data = np.loadtxt(mc_file)
    weight = data[:,weight_col-1]
    E_src = data[:, E_src_col-1]
    orig_power = prog.power_law  # CRBeam is always called with fixed power=1 to optimize cache usage 
    weight *= np.power(E_src/Emax, -(power-orig_power))
    output_file = f'{mc_file}_p{power}'
    np.savetxt(output_file, data, header=header.strip(), fmt='%.6g')
    return output_file


# In[ ]:


# this code will not execute program since files already exist on server
# prog.run(force_overwrite=False)


# In[ ]:


# one can upload cache explicitely by call
# prog.upload_cache()


# In[ ]:


print(prog.output_path)
get_ipython().system('ls {prog.output_path}')


# In[ ]:


mc_file = prog.output_path + "/photon"


# In[ ]:


# how the data looks like
get_ipython().system("cat {mc_file} | awk 'NR<=5'")


# In[ ]:


if Emax!=EminSource:
    mc_file = adjust_weights(mc_file, Gamma)


# In[ ]:


# how the data looks like
get_ipython().system("cat {mc_file} | awk 'NR<=5'")


# In[ ]:


#! . ./makeSpecE2.sh {mc_file}
# the code above will not work on MMODA as of Sep 30 2023
# here is workaround
subprocess.run(['bash', 'makeSpecE2.sh', mc_file])


# In[ ]:


# rotating the beam

if EGMF > 0:
    from mc_rotate import mc_rotate
    mc_rotated_file = mc_rotate(mc_file, jet_half_size, jet_direction, psf_deg=psf)
else:
    mc_rotated_file = mc_file
mc_rotated_file


# In[ ]:


# calculating the energy spectrum
#! . ./makeSpecE2.sh {mc_rotated_file}
subprocess.run(['bash', 'makeSpecE2.sh', mc_rotated_file])


# In[ ]:


# how the rotated data looks like
get_ipython().system("cat {mc_rotated_file} | awk 'NR<=5'")


# In[ ]:


# reading the source distance in Mpc from the data file
T_Mpc = lc.get_distance_Mpc(mc_file)
T_Mpc


# In[ ]:


# building the spectrum figure for total flux
spec_file = mc_file + '.spec'
spec_fig = plot_spec(spec_file, title='spectrum', Emin=Emin, Emax=Emax, ext='png', show=False, logscale=True)
spec_image = PictureProduct.from_file(spec_fig)


# In[ ]:


# building the spectrum figure for the flux within PSF
spec_rotated_file = mc_rotated_file+ '.spec'
spec_rotated_fig  = plot_spec(spec_rotated_file, title='spectrum', Emin=Emin, Emax=Emax, ext='png', show=False, logscale=True)
spec_rotated_image = PictureProduct.from_file(spec_rotated_fig)


# In[ ]:


spec_rotated_file


# In[ ]:


lc_params = dict(
    logscale=True,
    max_t=-99, # 8760000,
    n_points=30,
    psf=psf,
    min_n_particles=10,
    cut_0=True
)

# building the light curve figure
if EGMF > 0:
    light_curve_fig = lc.make_plot(mc_rotated_file, **lc_params)
    light_curve_image = PictureProduct.from_file(light_curve_fig)
else:
    # to avoid possible problems with absent figure
    light_curve_image = PictureProduct.from_file(spec_fig)


# In[ ]:


l_curve = None
if EGMF > 0:
    delay, weights = lc.get_counts(mc_rotated_file, **lc_params)
    t, f, N = lc.light_curve(delay, weights, **lc_params)
    l_curve = LightCurveDataProduct.from_arrays(
        times = t*3600., # t is in hours
        fluxes = f,
        errors = f/np.sqrt(N),
        time_format='unix'
    )


# In[ ]:


import seaborn as sns
if EGMF > 0:
    data = np.loadtxt(mc_rotated_file)
    non_zero_delays = data[:, 5]>0
    n_non_zero_delays = np.sum(non_zero_delays)
    print(f'{n_non_zero_delays/len(data)*100:g} % of {len(data)} events have delay > 0')
    if n_non_zero_delays > 0:
        delay_years = data[non_zero_delays, 5] * T_Mpc*3.2637977241e+06 # delay in years
        weight = data[non_zero_delays, 1] 
        #sns.distplot(np.log10(delay_years), bins=4, hist_kws={'weights':weight})
        sns.histplot(x=delay_years, weights=weight)


# In[ ]:


d=np.genfromtxt(spec_file)
ee=d[:,0]
ff=d[:,1]
ff_err=ff/np.sqrt(d[:,2])
plt.errorbar(ee,ff,yerr=ff_err,label='total flux')
d=np.genfromtxt(spec_rotated_file)
ee1=d[:,0]
ff1=d[:,1]
ff_err1=ff1/np.sqrt(d[:,2])
plt.errorbar(ee1,ff1,yerr=ff_err1,label='PSF flux')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Energy, eV')
plt.ylabel('$E^2dN/dE$, arbitrary units')
plt.legend(loc='lower right')
plt.savefig('Spectrum.png',format='png',bbox_inches='tight')

bin_image = PictureProduct.from_file('Spectrum.png')

data=[ee,ff,ff_err]
names=('Energy','Total_flux','Total_flux_err')
table = ODAAstropyTable(Table(data, names = names))
data=[ee,ff,ff_err]
names=('Energy','PSF_flux','PSF_flux_err')
table1 = ODAAstropyTable(Table(data, names = names))


# In[ ]:


flux_unit = u.eV/u.cm/u.cm/u.s/u.sr
spec = Spectrum1D(spectral_axis=ee*u.eV, flux=ff*flux_unit, uncertainty=StdDevUncertainty(ff_err))
spec.write('spec.fits', overwrite=True)
dp = NumpyDataProduct.from_fits_file('spec.fits')
spec_rotated = Spectrum1D(spectral_axis=ee1*u.eV, flux=ff1*flux_unit, uncertainty=StdDevUncertainty(ff_err1))
spec_rotated.write('spec_rotated.fits', overwrite=True)
dp_rotated = NumpyDataProduct.from_fits_file('spec_rotated.fits')


# In[ ]:


spectrum_png = bin_image # http://odahub.io/ontology#ODAPictureProduct
light_curve_png = light_curve_image # http://odahub.io/ontology#ODAPictureProduct
total_spectrum_table = table # http://odahub.io/ontology#ODAAstropyTable
psf_spectrum_table = table1 # http://odahub.io/ontology#ODAAstropyTable
lc_result = l_curve # http://odahub.io/ontology#LightCurve
spectrum = dp # https://odahub.io/ontology/#Spectrum
spectrum_rotated = dp_rotated # https://odahub.io/ontology/#Spectrum


# In[ ]:


spec_png = spec_image # http://odahub.io/ontology#ODAPictureProduct
spec_rotated_png = spec_rotated_image # http://odahub.io/ontology#ODAPictureProduct


# In[ ]:





# In[ ]:


_simple_outs, _oda_outs = [], []
_galaxy_meta_data = {}
_oda_outs.append(('out_CRbeam_spectrum_png', 'spectrum_png_galaxy.output', spectrum_png))
_oda_outs.append(('out_CRbeam_light_curve_png', 'light_curve_png_galaxy.output', light_curve_png))
_oda_outs.append(('out_CRbeam_total_spectrum_table', 'total_spectrum_table_galaxy.output', total_spectrum_table))
_oda_outs.append(('out_CRbeam_psf_spectrum_table', 'psf_spectrum_table_galaxy.output', psf_spectrum_table))
_oda_outs.append(('out_CRbeam_lc_result', 'lc_result_galaxy.output', lc_result))
_simple_outs.append(('out_CRbeam_spectrum', 'spectrum_galaxy.output', spectrum))
_simple_outs.append(('out_CRbeam_spectrum_rotated', 'spectrum_rotated_galaxy.output', spectrum_rotated))

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

