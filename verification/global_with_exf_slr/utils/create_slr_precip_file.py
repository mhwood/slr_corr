import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from scipy.interpolate import griddata
import shutil
import argparse

def read_balanced_precip_and_evap_files(config_dir):

    precip_file = os.path.join(config_dir,'data','ecco_precip_balanced.bin')
    precip_grid = np.fromfile(precip_file,'>f4')
    precip_grid = np.reshape(precip_grid,(12,40,90))

    evap_file = os.path.join(config_dir, 'data', 'ecco_evap.bin')
    evap_grid = np.fromfile(evap_file, '>f4')
    evap_grid = np.reshape(evap_grid, (12, 40, 90))

    lon = np.arange(2,360,4)
    lat = np.arange(-78,80,4)
    Lon, Lat = np.meshgrid(lon,lat)

    return(Lon, Lat, precip_grid, evap_grid)

def read_area_and_mask_grids(grid_run_dir,n_proc = 2):

    if n_proc == 1:
        grid_file_1 = grid_run_dir + '/mnc_0001/grid.t001.nc'

        ##############################################################
        # read in the area grids
        ds = nc4.Dataset(grid_file_1)
        rA_grid = ds.variables['rA'][:, :]
        mask_grid = ds.variables['HFacC'][:, :, :]
        ds.close()
        mask_grid = mask_grid[0, :, :]

    if n_proc == 2:
        if 'grid.t001.nc' in os.listdir(grid_run_dir + '/mnc_0001'):
            grid_file_1 = grid_run_dir + '/mnc_0001/grid.t001.nc'
            grid_file_2 = grid_run_dir + '/mnc_0002/grid.t002.nc'
        else:
            grid_file_1 = grid_run_dir + '/mnc_0001/grid.t002.nc'
            grid_file_2 = grid_run_dir + '/mnc_0002/grid.t001.nc'

        ##############################################################
        # read in the area grids
        ds = nc4.Dataset(grid_file_1)
        X1 = ds.variables['X'][:]
        mask_var1 = ds.variables['HFacC'][:,:,:]
        rA_var1 = ds.variables['rA'][:, :]
        ds.close()
        mask_var1 = mask_var1[0,:,:]

        ds = nc4.Dataset(grid_file_2)
        X2 = ds.variables['X'][:]
        mask_var2 = ds.variables['HFacC'][:, : ,:]
        rA_var2 = ds.variables['rA'][:, :]
        ds.close()
        mask_var2 = mask_var2[0, :, :]

        if np.max(X1) > np.max(X2):
            rA_grid = np.concatenate([rA_var2, rA_var1], axis=1)
            mask_grid = np.concatenate([mask_var2, mask_var1], axis=1)
        else:
            rA_grid = np.concatenate([rA_var1, rA_var2], axis=1)
            mask_grid = np.concatenate([mask_var1, mask_var2], axis=1)

    mask_grid[mask_grid > 0] = 1
    mask_grid = mask_grid.astype(int)

    return(rA_grid, mask_grid)

def adjust_ecco_grids_for_slr_accumulation(ecco_precip_grid, area_grid, mask_grid):

    slr_per_year = 0.0036 # m (3.6 mm/yr)
    slr_per_month = slr_per_year/12 # m
    seconds_in_month = 24 * 60 * 60 * 30 #s

    total_wet_area = np.sum(area_grid*mask_grid) #m2
    total_months = 12

    ecco_slr_precip_grid = np.zeros_like(ecco_precip_grid)

    extra_volume_in_month = total_wet_area*slr_per_month
    rate_correction = extra_volume_in_month / (total_wet_area * seconds_in_month)

    print('SLR in month (m):',extra_volume_in_month/total_wet_area)
    print('Extra volume per month (m3): ', extra_volume_in_month)
    print('Volumetric rate in entire domain (m3/s)', extra_volume_in_month/seconds_in_month)
    print('Total wet area (m2): ',total_wet_area)
    print('Mean rate per wet cell (m/s): ',rate_correction)

    for month in range(total_months):
        # if month == 0:
        #     C = plt.imshow(ecco_precip_grid[month,:,:],origin='lower')
        #     plt.colorbar(C)
        #     plt.show()
        month_grid = np.copy(ecco_precip_grid[month,:,:])
        month_grid[mask_grid>0] += rate_correction
        ecco_slr_precip_grid[month, :, :] = month_grid

    return(ecco_slr_precip_grid)

def get_accumulation_timeseries(precip_grid, evap_grid, area_grid, total_wet_area):
    total_accumulation = np.zeros(12, )
    seconds_in_month = 24 * 60 * 60 * 30
    for i in range(12):
        net_rate = precip_grid[i, :, :] - evap_grid[i, :, :]
        volume_accumulation = np.sum(net_rate*area_grid*seconds_in_month)
        accumulation = volume_accumulation/total_wet_area
        if i == 0:
            total_accumulation[i] = accumulation
        else:
            total_accumulation[i] = total_accumulation[i - 1] + accumulation
    return(total_accumulation)

def create_slr_file(config_dir):

    print(' - Reading the balanced precip file')
    Lon, Lat, ecco_precip_grid, ecco_evap_grid = read_balanced_precip_and_evap_files(config_dir)

    print(' - Reading the area grid from the spinup model')
    grid_run_dir = config_dir+'/run_spinup'
    area_grid, mask_grid = read_area_and_mask_grids(grid_run_dir,n_proc = 1)

    ecco_slr_precip_grid = adjust_ecco_grids_for_slr_accumulation(ecco_precip_grid, area_grid, mask_grid)

    # this manual correction adds a trend of 3 mm in 365 days
    # ecco_slr_precip_grid = np.copy(ecco_precip_grid)
    # ecco_slr_precip_grid[ecco_precip_grid != 0] += 3.3e-9

    # total_wet_area = np.sum(mask_grid*area_grid)
    # balanced_accumulation = get_accumulation_timeseries(ecco_precip_grid, ecco_evap_grid, area_grid, total_wet_area)
    # slr_accumulation = get_accumulation_timeseries(ecco_slr_precip_grid, ecco_evap_grid, area_grid, total_wet_area)
    # plt.plot(balanced_accumulation)
    # plt.plot(slr_accumulation)
    # plt.show()

    # for i in range(1,np.shape(ecco_evap_grid)[0]):
    #     plt.subplot(1,2,1)
    #     C = plt.imshow(ecco_evap_grid[i,:,:],origin='lower',vmin=0,vmax=1e-7)
    #     plt.colorbar(C)
    #     plt.title('Evap grid (month = '+str(i)+')')
    #
    #     plt.subplot(1, 2, 2)
    #     C = plt.imshow(ecco_precip_grid[i, :, :], origin='lower',vmin=0,vmax=1e-7)
    #     plt.colorbar(C)
    #     plt.title('Precip grid (month = ' + str(i) + ')')
    #
    #     plt.show()

    data_output_dir = os.path.join(config_dir,'data')
    precip_output_file = os.path.join(data_output_dir,'ecco_precip_slr.bin')
    ecco_slr_precip_grid.ravel('C').astype('>f4').tofile(precip_output_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the configuration is stored (../).", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    create_slr_file(config_dir)

