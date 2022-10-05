
import os
import netCDF4 as nc4
import numpy as np
import matplotlib.pyplot as plt
import argparse

def read_and_plot_grid_from_mnc_files(config_dir):

    run_dir = config_dir+'/run_spinup'

    if 'diagSurf.0000000000.t001.nc' in os.listdir(run_dir+'/mnc_0001'):
        file_1 = run_dir+'/mnc_0001/diagSurf.0000000000.t001.nc'
        file_2 = run_dir + '/mnc_0002/diagSurf.0000000000.t002.nc'
        grid_file_1 = run_dir + '/mnc_0001/grid.t001.nc'
        grid_file_2 = run_dir + '/mnc_0002/grid.t002.nc'
    else:
        file_1 = run_dir + '/mnc_0001/diagSurf.0000000000.t002.nc'
        file_2 = run_dir + '/mnc_0002/diagSurf.0000000000.t001.nc'
        grid_file_1 = run_dir + '/mnc_0001/grid.t002.nc'
        grid_file_2 = run_dir + '/mnc_0002/grid.t001.nc'

    ##############################################################
    # read in the EtaN grids
    ds = nc4.Dataset(file_1)
    X1 = ds.variables['X'][:]
    var1 = ds.variables['ETAN'][:, :, :, :]
    ds.close()
    var1 = var1[:,0,:,:]

    ds = nc4.Dataset(file_2)
    X2 = ds.variables['X'][:]
    var2 = ds.variables['ETAN'][:, :, :, :]
    ds.close()

    var2 = var2[:, 0, :, :]

    if np.max(X1)>np.max(X2):
        var_grid = np.concatenate([var2,var1],axis = 2)
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

    dry_area = np.sum(rA_grid[var_grid[-1,:,:]==0])
    wet_area = np.sum(rA_grid[var_grid[-1, :, :] != 0])
    total_area = np.sum(rA_grid)

    print(dry_area,total_area,dry_area/total_area)

    # plt.imshow(var[-1,:,:],origin='lower')
    #     # plt.show()

    time = np.arange(np.shape(var_grid)[0])
    time = time/365

    total_eta = np.sum(var_grid,axis=1)
    total_eta = np.sum(total_eta, axis=1)
    mean_eta = total_eta/total_area

    fig = plt.figure(figsize=(9,9))

    plt.subplot(3,1,1)
    plt.title('Mean EtaN (m)')
    plt.plot(time,mean_eta)

    plt.subplot(3,2,3)
    plt.pcolormesh(XC_grid,YC_grid,var_grid[365*5,:,:]-var_grid[365*4,:,:], vmin=-0.1, vmax = 0.1, cmap='seismic')
    plt.title('Year 5 - Year 4 (m)')

    plt.subplot(3, 2, 4)
    plt.pcolormesh(XC_grid, YC_grid,var_grid[365*10,:,:]-var_grid[365*9,:,:], vmin=-0.1, vmax = 0.1, cmap='seismic')
    plt.title('Year 10 - Year 9 (m)')

    plt.subplot(3, 2, 5)
    plt.pcolormesh(XC_grid, YC_grid,var_grid[365*15,:,:]-var_grid[365*14,:,:], vmin=-0.1, vmax = 0.1, cmap='seismic')
    plt.title('Year 15 - Year 14 (m)')

    plt.subplot(3, 2, 6)
    plt.pcolormesh(XC_grid, YC_grid,var_grid[365*20-1,:,:]-var_grid[365*19,:,:], vmin=-0.1, vmax = 0.1, cmap='seismic')
    plt.title('Year 20 - Year 19 (m)')

    if 'plots' not in os.listdir(config_dir):
        os.mkdir(config_dir+'/plots')

    plt.savefig(config_dir+'/plots/EtaN_spinup.png')
    plt.close(fig)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument("-d", "--config_dir", action="store",
                        help="The directory where the L1, L1, and L3 configurations are stored.", dest="config_dir",
                        type=str, required=True)

    args = parser.parse_args()
    config_dir = args.config_dir

    read_and_plot_grid_from_mnc_files(config_dir)