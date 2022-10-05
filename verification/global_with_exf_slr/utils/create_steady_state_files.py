import os
import numpy as np
import matplotlib.pyplot as plt
import netCDF4 as nc4
from scipy.interpolate import griddata
import shutil
import argparse

def read_ncep_emp_file(mitgcm_dir):

    emp_file = os.path.join(mitgcm_dir,'verification','tutorial_global_oce_latlon','input','ncep_emp.bin')
    emp_grid = np.fromfile(emp_file,'>f4')
    emp_grid = np.reshape(emp_grid,(12,40,90))

    lon = np.arange(2,360,4)
    lat = np.arange(-78,80,4)
    Lon, Lat = np.meshgrid(lon,lat)

    return(Lon, Lat, emp_grid)

def read_area_and_mask_grids(grid_run_dir,n_proc = 2):

    if n_proc == 1:
        file_1 = grid_run_dir + '/mnc_0001/diagSurf.00000' + timestep + '.t001.nc'
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

def read_ecco_precip_evap_grids(fw_flux_folder, Lon, Lat, emp_grid):

    ecco_evap_grid = np.zeros_like(emp_grid)
    ecco_precip_grid = np.zeros_like(emp_grid)

    for month in range(1,13):
        file_name = 'OCEAN_AND_ICE_SURFACE_FW_FLUX_mon_mean_2000-'+'{:02d}'.format(month)+'_ECCO_V4r4_latlon_0p50deg.nc'
        ds = nc4.Dataset(fw_flux_folder+'/'+file_name)
        if month==1:
            ecco_lat = ds.variables['latitude'][:]
            ecco_lon = ds.variables['longitude'][:]
            ecco_lon[ecco_lon<0] += 360
            ecco_Lon, ecco_Lat = np.meshgrid(ecco_lon, ecco_lat)
            points = np.column_stack([ecco_Lon.ravel(),ecco_Lat.ravel()])

        evap = ds.variables['EXFevap'][:,:]
        precip = ds.variables['EXFpreci'][:, :] + ds.variables['EXFroff'][:, :]

        indices = evap.ravel() < 1
        points_subset = points[indices,:]
        evap_subset = evap.ravel()
        evap_subset = evap_subset[indices]
        precip_subset = precip.ravel()
        precip_subset = precip_subset[indices]

        interp_evap = griddata(points_subset,evap_subset,(Lon,Lat),method = 'nearest')
        interp_evap[emp_grid[0,:,:]==0]=0
        ecco_evap_grid[month-1,:,:] = interp_evap

        interp_precip = griddata(points_subset, precip_subset, (Lon, Lat),method = 'nearest')
        interp_precip[emp_grid[0, :, :] == 0] = 0
        ecco_precip_grid[month - 1, :, :] = interp_precip

    return(ecco_evap_grid,ecco_precip_grid)

def adjust_ecco_grids_for_zero_net_annual_accumulation(ecco_evap_grid, ecco_precip_grid, area_grid, mask_grid):

    ecco_evap_grid[ecco_evap_grid < 0] = 0

    total_months = 12
    total_wet_area = np.sum(area_grid*mask_grid)

    print('  - Iterating to balance the precip grid')
    for i in range(10):
        print('    - Iteration '+str(i))
        # calculate total flux (m3/s)
        total_precip_volume_flux = 0
        total_evap_volume_flux = 0
        for month in range(total_months):
            total_precip_volume_flux += np.sum(ecco_precip_grid[month,:,:] * area_grid * mask_grid)
            total_evap_volume_flux += np.sum(ecco_evap_grid[month, :, :] * area_grid * mask_grid)
        total_imbalance_flux = total_precip_volume_flux - total_evap_volume_flux

        print('      - Total flux imbalance before correction: ' + str(total_imbalance_flux) +' m3/s')

        correction = total_imbalance_flux/(total_wet_area*total_months)
        print('      - Correction: '+str(correction)+' m3/s')

        for month in range(total_months):
            new_month_precip_grid = np.copy(ecco_precip_grid[month,:,:])
            new_month_precip_grid[mask_grid>0] -= correction
            ecco_precip_grid[month, :, :] = new_month_precip_grid
        ecco_precip_grid[ecco_precip_grid < 0] = 0

    return(ecco_evap_grid,ecco_precip_grid)

def copy_tutorial_files(mitgcm_dir,data_output_dir):

    for file_name in ['lev_s.bin','lev_sss.bin','lev_t.bin','lev_sst.bin','trenberth_taux.bin','trenberth_tauy.bin','bathymetry.bin']:
        if file_name not in os.listdir(data_output_dir):
            src_file = os.path.join(mitgcm_dir,'verification','tutorial_global_oce_latlon','input',file_name)
            dst_file = os.path.join(data_output_dir,file_name)
            shutil.copyfile(src_file,dst_file)

def create_steady_state_files(mitgcm_dir,data_input_dir,data_output_dir):

    print(' - Reading the NCEP E-P file from the tutorial')
    Lon, Lat, emp_grid = read_ncep_emp_file(mitgcm_dir)

    print(' - Reading the area and mask grids')
    grid_run_dir = mitgcm_dir+'/configurations/precip_slr_corr/run_spinup'
    area_grid, mask_grid = read_area_and_mask_grids(grid_run_dir,n_proc=2)

    print(' - Reading the ECCO evap and precip files')
    ecco_evap_grid, ecco_precip_grid = read_ecco_precip_evap_grids(data_input_dir, Lon, Lat, emp_grid)

    ecco_evap_grid, ecco_precip_grid = adjust_ecco_grids_for_zero_net_annual_accumulation(ecco_evap_grid,ecco_precip_grid,area_grid, mask_grid)

    print(' - Sanity check: ')
    print('    - Nonzero precip values: '+str(np.sum(ecco_precip_grid<0)))
    print('    - Nonzero evap values: ' + str(np.sum(ecco_evap_grid < 0)))

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

    evap_output_file = os.path.join(data_output_dir,'ecco_evap.bin')
    ecco_evap_grid.ravel('C').astype('>f4').tofile(evap_output_file)

    precip_output_file = os.path.join(data_output_dir,'ecco_precip_balanced.bin')
    ecco_precip_grid.ravel('C').astype('>f4').tofile(precip_output_file)

    copy_tutorial_files(mitgcm_dir, data_output_dir)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--mitgcm_dir", action="store",
                        help="The directory where the MITgcm clone is stored.", dest="mitgcm_dir",
                        type=str, required=True)

    parser.add_argument("-i", "--data_input_dir", action="store",
                        help="The directory where ECCO FE_FLUX files for year 2000 are stored.", dest="data_input_dir",
                        type=str, required=True)

    parser.add_argument("-o", "--data_output_dir", action="store",
                        help="The directory where output will be stored (data).", dest="data_output_dir",
                        type=str, required=True)

    args = parser.parse_args()
    mitgcm_dir = args.mitgcm_dir
    data_input_dir = args.data_input_dir
    data_output_dir = args.data_output_dir

    create_steady_state_files(mitgcm_dir,data_input_dir,data_output_dir)

