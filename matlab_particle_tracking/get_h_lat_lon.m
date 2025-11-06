%%   grabs h and lat_rho & lon_rho from nc data file
      %sid = input('use 1] depth avg. 2] surface   ');
      sid = 1;
      if sid == 1,
       fn='../../../data/gbe_0001_his_vars_test_30days_uvbar_untilt.nc'; % file name for model output (velocities)
       else 
	fn = '../../../data/gbe_0001_his_vars_30days_u_v_untilt.nc';
       end 
	fnz='../../../data/gbe_0001_his_vars_test_30days_zeta.nc'; % file name for model output (velocities)

%%t=ncread(fn,'ocean_time');
otu = ncreadatt(fn, 'ocean_time', 'units');
dnuvs = datenum(str2num(otu(15:18)), str2num(otu(20:21)), str2num(otu(23:24)), str2num(otu(26:27)), str2num(otu(29:30)), str2num(otu(32:33)));
lon_rho = ncread(fn,'lon_rho');
[xi_rho, eta_rho] = size(lon_rho);
lat_rho = ncread(fn,'lat_rho');
h = ncread(fn, 'h');
[Nr, Er] = lltoutm(lat_rho,lon_rho,23);
h_bot = 24-h;

outfile = ['./bathymetry_from_gbe_0001_his_vars.nc'];

nccreate(outfile, 'dnuvs', 'Dimensions', {'ndnuvs', 1}, 'Datatype', 'double');
ncwrite(outfile, 'dnuvs', double(dnuvs));
ncwriteatt(outfile, 'dnuvs', 'long_name', 'Start Datenum for data run');
ncwriteatt(outfile, 'dnuvs', 'units', 'Days (Datenum)');
ncwriteatt(outfile,'/', 'velocity_filename',fn);	

ncwriteatt(outfile, '/', 'source', 'get_h_lat_lon.m');

nccreate(outfile, 'lon_rho', 'Dimensions', {'xi_rho', xi_rho, 'eta_rho', eta_rho}, 'Datatype', 'double');
ncwriteatt(outfile, 'lon_rho', 'units', 'deg E');
ncwriteatt(outfile, 'lon_rho', 'long_name', 'longitude at RHO-points');
ncwriteatt(outfile, 'lon_rho', 'type', 'data');
ncwrite(outfile, 'lon_rho', lon_rho);

nccreate(outfile, 'lat_rho', 'Dimensions', {'xi_rho', xi_rho, 'eta_rho', eta_rho}, 'Datatype', 'double');
ncwriteatt(outfile, 'lat_rho', 'units', 'deg N');
ncwriteatt(outfile, 'lat_rho', 'long_name', 'latitude at RHO-points');
ncwriteatt(outfile, 'lat_rho', 'type', 'data');
ncwrite(outfile, 'lat_rho', lat_rho);

nccreate(outfile, 'h', 'Dimensions', {'xi_rho', xi_rho, 'eta_rho', eta_rho}, 'Datatype', 'double');
ncwriteatt(outfile, 'h', 'units', 'meter');
ncwriteatt(outfile, 'h', 'grid', 'grid');
ncwriteatt(outfile, 'h', 'long_name', 'grid_bathymetry at RHO-points');
ncwriteatt(outfile, 'h', 'type', 'data');
ncwrite(outfile, 'h', h);

nccreate(outfile, 'h_bot', 'Dimensions', {'xi_rho', xi_rho, 'eta_rho', eta_rho}, 'Datatype', 'double');
ncwriteatt(outfile, 'h_bot', 'units', 'meter');
ncwriteatt(outfile, 'h_bot', 'long_name', 'bathymetric elevation NAVD88 at RHO-points');
ncwriteatt(outfile, 'h_bot', 'type', 'data');
ncwrite(outfile, 'h_bot', h_bot);



