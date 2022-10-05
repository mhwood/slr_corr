
import os
import netCDF4 as nc4
import numpy as np
import matplotlib.pyplot as plt
import argparse
from datetime import datetime, timedelta

def read_etan_field(run_dir):

    if 'diagsEXF.0000014600.t001.nc' in os.listdir(run_dir + '/mnc_0001'):
        file_1 = run_dir + '/mnc_0001/diagsEXF.0000014600.t001.nc'
        file_2 = run_dir + '/mnc_0002/diagsEXF.0000014600.t002.nc'
        grid_file_1 = run_dir + '/mnc_0001/grid.t001.nc'
        grid_file_2 = run_dir + '/mnc_0002/grid.t002.nc'
    else:
        file_1 = run_dir + '/mnc_0001/diagsEXF.0000014600.t002.nc'
        file_2 = run_dir + '/mnc_0002/diagsEXF.0000014600.t001.nc'
        grid_file_1 = run_dir + '/mnc_0001/grid.t002.nc'
        grid_file_2 = run_dir + '/mnc_0002/grid.t001.nc'

    ##############################################################
    # read in the EtaN grids
    ds = nc4.Dataset(file_1)
    X1 = ds.variables['X'][:]
    var1 = ds.variables['EXFempmr'][:, :, :, :]
    ds.close()
    var1 = var1[:, 0, :, :]

    ds = nc4.Dataset(file_2)
    X2 = ds.variables['X'][:]
    var2 = ds.variables['EXFempmr'][:, :, :, :]
    ds.close()

    var2 = var2[:, 0, :, :]

    if np.max(X1) > np.max(X2):
        var_grid = np.concatenate([var2, var1], axis=2)
    else:
        var_grid = np.concatenate([var1, var2], axis=2)

    ##############################################################
    # read in the area grids
    ds = nc4.Dataset(grid_file_1)
    X1 = ds.variables['X'][:]
    rA_var1 = ds.variables['rA'][:, :]
    XC_var1 = ds.variables['XC'][:, :]
    YC_var1 = ds.variables['YC'][:, :]
    ds.close()

    ds = nc4.Dataset(grid_file_2)
    X2 = ds.variables['X'][:]
    rA_var2 = ds.variables['rA'][:, :]
    XC_var2 = ds.variables['XC'][:, :]
    YC_var2 = ds.variables['YC'][:, :]
    ds.close()

    if np.max(X1) > np.max(X2):
        rA_grid = np.concatenate([rA_var2, rA_var1], axis=1)
        XC_grid = np.concatenate([XC_var2, XC_var1], axis=1)
        YC_grid = np.concatenate([YC_var2, YC_var1], axis=1)
    else:
        rA_grid = np.concatenate([rA_var1, rA_var2], axis=1)
        XC_grid = np.concatenate([XC_var1, XC_var2], axis=1)
        YC_grid = np.concatenate([YC_var1, YC_var2], axis=1)

    return(XC_grid,YC_grid,var_grid)

def calculate_EtaN_timeseries(EtaN):

    startdate = datetime(1991,12,27)
    start_dec_yr = 1991+(361/365.0)

    timeseries = np.zeros((np.shape(EtaN)[0],2))
    non_zero_indices = EtaN[0,:,:]!=0
    for i in range(np.shape(EtaN)[0]):
        dec_yr = start_dec_yr+(i/365.25)
        timeseries[i,0] = dec_yr
        step_field = EtaN[i,:,:]
        timeseries[i,1] = np.mean(step_field[non_zero_indices])

    timeseries[:,1] -= timeseries[0,1]

    return(timeseries)

def plot_EtaN_timeseries(config_dir,EtaN_slr_uncorrected_timeseries):

    fig = plt.figure(figsize=(9,4))

    plt.plot(EtaN_slr_uncorrected_timeseries[:,0],EtaN_slr_uncorrected_timeseries[:,1])

    if 'plots' not in os.listdir(config_dir):
        os.mkdir(config_dir+'/plots')

    plt.grid()

    plt.savefig(config_dir+'/plots/EXF_timeseries_comparison.png')
    plt.close(fig)

def plot_etan_comparison(config_dir):

    model_1 = 'precip_balanced'
    run_dir = config_dir+'/run_'+model_1
    XC, YC, EtaN_slr_uncorrected = read_etan_field(run_dir)
    EtaN_slr_uncorrected_timeseries = calculate_EtaN_timeseries(EtaN_slr_uncorrected)

    plot_EtaN_timeseries(config_dir, EtaN_slr_uncorrected_timeseries)






if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where run dirs are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    plot_etan_comparison(config_dir)