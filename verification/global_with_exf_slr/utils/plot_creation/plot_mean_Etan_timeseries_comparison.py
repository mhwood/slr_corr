
import os
import netCDF4 as nc4
import numpy as np
import matplotlib.pyplot as plt
import argparse
from datetime import datetime, timedelta

def read_etan_field(run_dir,spinup=False,n_proc=2):

    if spinup:
        timestep = '00000'
    else:
        timestep = '14600'

    if n_proc ==1:
        file_1 = run_dir + '/mnc_0001/diagSurf.00000'+timestep+'.t001.nc'
        file_2 = run_dir + '/mnc_0002/diagSurf.00000'+timestep+'.t002.nc'
        grid_file_1 = run_dir + '/mnc_0001/grid.t001.nc'
        grid_file_2 = run_dir + '/mnc_0002/grid.t002.nc'

        ##############################################################
        # read in the EtaN grids
        ds = nc4.Dataset(file_1)
        var_grid = ds.variables['ETAN'][:, :, :, :]
        ds.close()

        var_grid = var_grid[:, 0, :, :]

        ##############################################################
        # read in the area grids
        ds = nc4.Dataset(grid_file_1)
        XC_grid = ds.variables['XC'][:, :]
        YC_grid = ds.variables['YC'][:, :]
        rA_grid = ds.variables['rA'][:, :]
        mask_grid = ds.variables['HFacC'][:, :, :]
        ds.close()

        mask_grid = mask_grid[0,:,:]


    if n_proc ==2:
        if 'diagSurf.00000'+timestep+'.t001.nc' in os.listdir(run_dir + '/mnc_0001'):
            file_1 = run_dir + '/mnc_0001/diagSurf.00000'+timestep+'.t001.nc'
            file_2 = run_dir + '/mnc_0002/diagSurf.00000'+timestep+'.t002.nc'
            grid_file_1 = run_dir + '/mnc_0001/grid.t001.nc'
            grid_file_2 = run_dir + '/mnc_0002/grid.t002.nc'
        else:
            file_1 = run_dir + '/mnc_0001/diagSurf.00000'+timestep+'.t002.nc'
            file_2 = run_dir + '/mnc_0002/diagSurf.00000'+timestep+'.t001.nc'
            grid_file_1 = run_dir + '/mnc_0001/grid.t002.nc'
            grid_file_2 = run_dir + '/mnc_0002/grid.t001.nc'

        ##############################################################
        # read in the EtaN grids
        ds = nc4.Dataset(file_1)
        X1 = ds.variables['X'][:]
        var1 = ds.variables['ETAN'][:, :, :, :]
        ds.close()
        var1 = var1[:, 0, :, :]

        ds = nc4.Dataset(file_2)
        X2 = ds.variables['X'][:]
        var2 = ds.variables['ETAN'][:, :, :, :]
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
        mask_var1 = ds.variables['HFacC'][:,:,:]
        rA_var1 = ds.variables['rA'][:, :]
        XC_var1 = ds.variables['XC'][:, :]
        YC_var1 = ds.variables['YC'][:, :]
        ds.close()
        mask_var1 = mask_var1[0,:,:]

        ds = nc4.Dataset(grid_file_2)
        X2 = ds.variables['X'][:]
        mask_var2 = ds.variables['HFacC'][:, : ,:]
        rA_var2 = ds.variables['rA'][:, :]
        XC_var2 = ds.variables['XC'][:, :]
        YC_var2 = ds.variables['YC'][:, :]
        ds.close()
        mask_var2 = mask_var2[0, :, :]

        if np.max(X1) > np.max(X2):
            rA_grid = np.concatenate([rA_var2, rA_var1], axis=1)
            XC_grid = np.concatenate([XC_var2, XC_var1], axis=1)
            YC_grid = np.concatenate([YC_var2, YC_var1], axis=1)
            mask_grid = np.concatenate([mask_var2, mask_var1], axis=1)
        else:
            rA_grid = np.concatenate([rA_var1, rA_var2], axis=1)
            XC_grid = np.concatenate([XC_var1, XC_var2], axis=1)
            YC_grid = np.concatenate([YC_var1, YC_var2], axis=1)
            mask_grid = np.concatenate([mask_var1, mask_var2], axis=1)

    mask_grid[mask_grid > 0] = 1
    mask_grid = mask_grid.astype(int)

    return(XC_grid,YC_grid,var_grid, rA_grid, mask_grid)

def calculate_EtaN_timeseries(EtaN,area_grid, mask_grid, start_dec_yr):

    timeseries = np.zeros((np.shape(EtaN)[0],2))
    non_zero_indices = EtaN[0,:,:]!=0
    for i in range(np.shape(EtaN)[0]):
        dec_yr = start_dec_yr+(i/365.25)
        timeseries[i,0] = dec_yr
        # step_field = EtaN[i,:,:]
        timeseries[i,1] = np.sum(area_grid*EtaN[i,:,:]*mask_grid)/np.sum(area_grid*mask_grid)

    # timeseries[:,1] -= timeseries[0,1]

    return(timeseries)

def read_observation_timeseries(config_dir,start_dec_yr):

    obs = np.fromfile(config_dir+'/data/slr_observations.bin','>f4')
    time = np.arange(len(obs))*(1/365.25)+start_dec_yr
    observation_timeseries = np.column_stack([time,obs])

    return(observation_timeseries)

def plot_EtaN_timeseries(config_dir,observation_timeseries,EtaN_slr_uncorrected_timeseries,
                         EtaN_slr_corrected_timeseries,EtaN_slr_spinup):

    fig = plt.figure(figsize=(9,4))

    # plt.plot(EtaN_slr_spinup[:, 0], EtaN_slr_spinup[:, 1], label='spinup')
    plt.plot(EtaN_slr_uncorrected_timeseries[:,0],EtaN_slr_uncorrected_timeseries[:,1],label='Without slr_corr')
    plt.plot(observation_timeseries[:, 0], observation_timeseries[:, 1], label='"Observations"', alpha=0.4)
    plt.plot(EtaN_slr_corrected_timeseries[:, 0], EtaN_slr_corrected_timeseries[:, 1], label='With slr_corr',linewidth=2)


    plt.legend(ncol=3)

    if 'plots' not in os.listdir(config_dir):
        os.mkdir(config_dir+'/plots')

    plt.text(1992, 0.00075, 'precip field tuned online\nto match "observations" provided',color= 'green')

    plt.text(1996,0,'Results from simple toy model, similar to\nthe global_with_exf verification experiment')

    plt.ylabel('Mean EtaN (m)')

    # plt.gca().set_xlim([1985,2003])
    # plt.gca().set_ylim([-0.14, -0.1])

    plt.grid(linestyle='--',alpha=0.4)

    plt.title('Demo of MITgcm pkg "slr_corr"')

    plt.savefig(config_dir+'/plots/EtaN_timeseries_comparison.png')
    plt.close(fig)

def plot_etan_comparison(config_dir):

    # model_2 = 'precip_balanced'
    # run_dir = config_dir + '/run_' + model_2
    # XC, YC, EtaN_precip_balanced = read_etan_field(run_dir)
    # EtaN_precip_balanced_timeseries = calculate_EtaN_timeseries(EtaN_precip_balanced)

    model_1 = 'slr_uncorrected'
    run_dir = config_dir + '/run_' + model_1
    XC, YC, EtaN_slr_uncorrected, area_grid, mask_grid = read_etan_field(run_dir,n_proc=1)
    start_dec_yr = 1991 + (361 / 365.0)
    EtaN_slr_uncorrected_timeseries = calculate_EtaN_timeseries(EtaN_slr_uncorrected, area_grid, mask_grid, start_dec_yr)

    model_2 = 'spinup'
    run_dir = config_dir+'/run_'+model_2
    XC, YC, EtaN_spinup, area_grid, mask_grid = read_etan_field(run_dir,spinup=True,n_proc=1)
    start_dec_yr = 1972
    EtaN_slr_spinup = calculate_EtaN_timeseries(EtaN_spinup, area_grid, mask_grid, start_dec_yr)

    model_3 = 'slr_corrected'
    run_dir = config_dir + '/run_' + model_3
    XC, YC, EtaN_corrected_timeseries, area_grid, mask_grid = read_etan_field(run_dir, n_proc=1)
    start_dec_yr = 1991 + (361 / 365.0)
    EtaN_slr_corrected_timeseries = calculate_EtaN_timeseries(EtaN_corrected_timeseries, area_grid, mask_grid, start_dec_yr)

    start_dec_yr = 1991 + (361 / 365.0)
    observation_timeseries = read_observation_timeseries(config_dir,start_dec_yr)

    plot_EtaN_timeseries(config_dir, observation_timeseries, EtaN_slr_uncorrected_timeseries,EtaN_slr_corrected_timeseries, EtaN_slr_spinup)







if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where run dirs are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    plot_etan_comparison(config_dir)