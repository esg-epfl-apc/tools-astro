import numpy as np
import sys
import os

import matplotlib.pyplot as plt

def plot_spec(spec_file, title='', Emin=0, Emax=0, ext='pdf', show=True, logscale=True):
    data = np.loadtxt(spec_file)
    if Emin > 0:
        data = data[Emin <= data[:,0]]
    else:
        Emin = 10**np.floor(np.log10(np.min(data[:,0])))
    if Emax > 0:
        data = data[data[:,0] <= Emax]
    else:
        Emax = 10**np.ceil(np.log10(np.max(data[:,0])))
        
    ymin = np.min(data[:,1] - data[:,1]/np.sqrt(data[:,2]))
    ymax = np.max(data[:,1] + data[:,1]/np.sqrt(data[:,2]))
        
    if title:
        plt.title(title)
    if logscale:
        plt.xscale('log')
        plt.yscale('log')
        if ymin <= 0:
            ymin = data[:,1] - data[:,1]/np.sqrt(data[:,2])
            ymin = ymin[ymin>0]
            if len(ymin) > 0:
                ymin = np.min(ymin)
            else:
                ymin = 0
        if ymin > 0:
            ymin = 10**np.floor(np.log10(ymin))
            ymax = 10**np.ceil(np.log10(ymax))
            plt.ylim(ymin, ymax)
    else:
        plt.ylim(ymin-0.1*(ymax-ymin), ymax+0.1*(ymax-ymin))
            
    
    plt.xlim(Emin, Emax)
        
    #plt.plot(data[:,2])
    plt.errorbar(data[:,0], data[:,1], yerr=data[:,1]/np.sqrt(data[:,2]), fmt='-o')
    out_file = spec_file + '.' + ext
    plt.savefig(out_file)
    if show:
        plt.show()
    return out_file