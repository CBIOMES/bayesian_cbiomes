import netCDF4 as nc4
import numpy as np
from dateutil.parser import parse
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

'''
This script regrids the data in SeaFlow_SizeDist.csv into a format
compatible with the matrix population model. The parameters v_min,
m and delta_v_inv are set below, m and delta_v_inv can also be 
specified as input parameters. If create_plots is set to true, both
the old and new size grid are plotted along with the old and regridded 
data set. The regridded data is saved in NetCDF format.
'''

def read_csv(fname):
    with open(fname) as f:

        header = f.readline()
        header_parts=header.strip().split(',')
        
        v_info = np.array([float(x.strip('(]" ')) for x in header_parts[2:]])
        
        data = {'t_min':[], 'par':[], 'data':[]}
        
        for i,line in enumerate(f):
            line_parts = line.strip().split(',')
            
            date = parse(line_parts[0])
            if i == 0:
                date0 = date
                #m = len(line_parts)-2
           
            data['t_min'].append((date-date0).total_seconds()//60)
            data['par'].append(float(line_parts[1]))
            data['data'].append([float(x) for x in line_parts[2:]])

    for k in data:
        data[k] = np.array(data[k])
    
    data['v_edges'] = np.empty(len(v_info)//2+1)
    data['v_edges'][:-1] = v_info[::2]
    data['v_edges'][-1] = v_info[-1]

    return data

data = read_csv('SeaFlow_SizeDist.csv')
create_plots = True

#
# specify v
#

v_min = 0.0135
#m = 20
#delta_v_inv = 5
m = 15
delta_v_inv = 5
import sys
if len(sys.argv) > 1:
    m = int(sys.argv[1])
    if len(sys.argv) > 2:
        delta_v_inv = int(sys.argv[2])

delta_v = 1.0/delta_v_inv
v = v_min * 2**(np.arange(m+1)*delta_v) # to get m intervals, we need m+1 edges

if create_plots:
    fig,ax = plt.subplots()
    ax.plot(data['v_edges'], np.ones_like(data['v_edges']), marker='s', label='data v edges')
    ax.plot(v, np.zeros_like(v), marker='o', label='model v edges')
    ax.legend()
    ax.set(ylim=(-1,2), yticks=[])
    #plt.show()

# data['data'] is num_t x num_v

def regrid(v0, w0, v1, permit_smallgrid=False):
    # assert that old grid is contained in new one
    assert v1[0] <= v0[0]
    if not permit_smallgrid:
        assert v1[-1] >= v0[-1]
    # check sizes
    if v0.size != w0.shape[0]+1:
        raise ValueError('Size of v0 ({}) does not match shape of w0 {}.'.format(v0.size, w0.shape))

    w1 = np.zeros((v1.size-1,w0.shape[1]))
    #v_data = data['v_edges']

    i1 = 1 # first right edge
    for i0 in range(1,v0.size):
        if v0[i0] < v1[i1]:
            if v0[i0-1] >= v1[i1-1]:
                w1[i1-1,:] += w0[i0-1,:]
            else:
                raise RuntimeError('This should not happen.')
        elif i1+1 == v1.size:
            # can only happen for permit_smallgrid
            a = (v1[i1]-v0[i0-1])/(v0[i0]-v0[i0-1])
            # left side
            w1[i1-1,:] += w0[i0-1,:]*a
            # right side is not covered by v1
        elif v0[i0] < v1[i1+1]:
            a = (v1[i1]-v0[i0-1])/(v0[i0]-v0[i0-1])
            # left side
            w1[i1-1,:] += w0[i0-1,:]*a
            # right side
            w1[i1,:] += w0[i0-1,:]*(1-a)
            i1 += 1
        else:
            raise NotImplementedError('Case not yet covered, chose a v_min smaller or equal to smallest value in data grid.')
    
    if v1[-1] >= v0[-1]:
        assert np.all(np.abs(np.sum(w0,axis=0)-np.sum(w1,axis=0))<1e-9)
    else:
        print('maximum loss of data due to reduction in grid coverage: {}'.format(np.max(np.abs(np.sum(w0,axis=0)-np.sum(w1,axis=0)))))
    return w1

w = regrid(data['v_edges'], data['data'].T, v, permit_smallgrid=True)

def add_colorbar(ax, **cbarargs):
    axins_cbar = inset_axes(ax, width='3%', height='90%', loc=5, bbox_to_anchor=(0.05,0.0,1,1), bbox_transform=ax.transAxes)
    mpl.colorbar.ColorbarBase(axins_cbar, orientation='vertical', **cbarargs)

if create_plots:
    fig,axs = plt.subplots(3,1,sharex=True,figsize=(12,10))
    ax = axs[0]
    ax.plot(data['t_min'], data['par'], color='gold')
    ax.set(ylabel='PAR (µmol photons m$^{-2}$ s$^{-1}$)')
    ax = axs[1]
    pc = ax.pcolormesh(data['t_min'], data['v_edges'], data['data'].T/np.sum(data['data'],axis=1), rasterized=True)
    add_colorbar(ax, label='size class proportion', norm=pc.norm, cmap=pc.cmap)
    ax.set(ylabel='size ($\mu$m$^3$)')
    ax = axs[2]
    pc = ax.pcolormesh(data['t_min'], v, w/np.sum(w,axis=0), rasterized=True)
    add_colorbar(ax, label='size class proportion', norm=pc.norm, cmap=pc.cmap)
    ax.set(xlabel='minutes since start', ylabel='size ($\mu$m$^3$)', ylim=axs[1].get_ylim())
    ax.text(0.01, 0.98, 'm = {}, delta_v_inv = {}'.format(m, delta_v_inv), ha='left', va='top', transform=ax.transAxes, size=16)
    
    #fig.patch.set_alpha(0.0)
    #fname = 'SeaFlow_SizeDist_regrid-{}-{}.pdf'.format(m, delta_v_inv)
    #fig.savefig(fname, bbox_inches='tight')

#
# write to netCDF
#

fname = 'SeaFlow_SizeDist_regrid-{}-{}.nc'.format(m, delta_v_inv)
with nc4.Dataset(fname,'w') as nc:
    nc.createDimension('time', w.shape[1])
    nc.createDimension('size', w.shape[0])
    nc.createDimension('size_bounds', w.shape[0]+1)

    nc.createVariable('time', int, ('time',), fill_value=False)
    nc.variables['time'][:] = data['t_min']
    nc.variables['time'].units = 'minutes since start of experiment'

    nc.createVariable('size', float, ('size',), fill_value=False)
    nc.variables['size'][:] = v[:-1]
    nc.variables['size'].units = 'um3'
    nc.variables['size'].long_name = 'the lower bound for each size class (used as the size of each size class in the model)'

    nc.createVariable('size_bounds', float, ('size_bounds',), fill_value=False)
    nc.variables['size_bounds'][:] = v
    nc.variables['size_bounds'].units = 'um3'
    nc.variables['size_bounds'].long_name = 'the lower and upper bounds for each size class'

    nc.createVariable('w_obs', float, ('size','time'), fill_value=False)
    nc.variables['w_obs'][:] = w/np.sum(w,axis=0)
    nc.variables['w_obs'].units = 'proportion'
    nc.variables['w_obs'].long_name = 'proportion of the population in each size class'

    nc.createVariable('PAR', float, ('time',), fill_value=False)
    nc.variables['PAR'][:] = data['par']
    nc.variables['PAR'].units = 'umol photons/m2/s'

    nc.createVariable('m', int, fill_value=False)
    nc.variables['m'][:] = m
    nc.variables['m'].long_name = 'the number of size classes'

    nc.createVariable('delta_v_inv', int, fill_value=False)
    nc.variables['delta_v_inv'][:] = delta_v_inv
    nc.variables['delta_v_inv'].long_name = 'the inverse of delta_v, spacifying the logarithmic scaling of size classes'

    nc.createVariable('v_min', float, fill_value=False)
    nc.variables['v_min'][:] = v_min
    nc.variables['v_min'].long_name = 'the lower bound of the smallest size class (used as the size of the smallest size class in the model)'

if create_plots:
    plt.show()
