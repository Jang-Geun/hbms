function [gdone] = write_netcdf_append_func(outfile, outpath, T, Xll, Yll, Xcart, Ycart, start, pid);
%%function [gdone] = write_netcdf_append_func(outfile, outpath, T, Xll, Yll, Xcart, Ycart, start, pid);
%
%	writes particle tracking data to outfile and appends as needed.
%	will not overwrite file if start == 1 and the outfile exists.
%
%	inputs:		outfile	= netcdf filename to store data, e.g., outfile = 'traj.nc'
%			outpath = path to output ncfile
%			T	= time vector (in seconds)
%			Xll	= longitude of particles (deg E)
%			Yll	= latitude of particles (deg N)
%			Xcart = cartesian coords of particles
%			Ycart = cartesian coords of particles
%			start	= time index to start writing (1 for first write)
%			pid	= 0] dont print to screen  1] print info to screen
%
%	outputs		gdone	= 0] outfile already exists, return without writing  1] wrote file ok
%

if (0),
	outpath = '/home/field/achilds/great_bay_modeling/image_code/';
	outfile = 'N160ktest_1_traj.mat';
	T =  1:.5:5;
	tv = 1:length(T);
	Xll = [(101:110)' (201:210)'];
	Yll = [(301:310)' (401:410)'];
	Xcart = [(101:110)' (201:210)'];
	Ycart = [(301:310)' (401:410)'];
	start = 1;
end;

[nt, npts] = size(Xll);
test = isempty(Xcart);

gdone = 0;
outfile = [outpath outfile];

%% check to see if file exists when start == 1
if start == 1,
	fid = fopen(outfile, 'r');
	if fid > 0,
		fclose(fid);
		fprintf(1, 'outfile %s already exists.  wont overwrite.\n', outfile);
		return;
	end;
end;

if (start == 1),
	nccreate(outfile, 'T', 'Dimensions', {'time', Inf}, 'Datatype', 'double', 'Format', '64bit');
	ncwriteatt(outfile, 'T', 'long_name', 'time');
	ncwriteatt(outfile, 'T', 'units', 'seconds since start of run');
	ncwriteatt(outfile, 'T', 'time_zone', 'none');

	nccreate(outfile, 'Xll', 'Dimensions', {'npts', npts, 'time', nt}, 'Datatype', 'double');
	ncwriteatt(outfile, 'Xll', 'units', 'deg E');
	ncwriteatt(outfile, 'Xll', 'long_name', 'longitude deg E');
	ncwriteatt(outfile, 'Xll', 'type', 'data');

	nccreate(outfile, 'Yll', 'Dimensions', {'npts', npts, 'time', nt}, 'Datatype', 'double');
	ncwriteatt(outfile, 'Yll', 'units', 'deg N');
	ncwriteatt(outfile, 'Yll', 'long_name', 'latitude deg N');
	ncwriteatt(outfile, 'Yll', 'type', 'data');
        if test == 0 
		nccreate(outfile, 'Xcart', 'Dimensions', {'npts', npts, 'time', nt}, 'Datatype', 'double');
		ncwriteatt(outfile, 'Xcart', 'units', 'meters');
		ncwriteatt(outfile, 'Xcart', 'long_name', 'X-dir cartesian coords in ROMS');
		ncwriteatt(outfile, 'Xcart', 'type', 'data');
	
		nccreate(outfile, 'Ycart', 'Dimensions', {'npts', npts, 'time', nt}, 'Datatype', 'double');
		ncwriteatt(outfile, 'Ycart', 'units', 'meters');
		ncwriteatt(outfile, 'Ycart', 'long_name', 'Y-dir cartesian coords in ROMS');
		ncwriteatt(outfile, 'Ycart', 'type', 'data');
        end

	ncwriteatt(outfile, '/', 'source', 'interp_Npts_func');
end;
ncwrite(outfile, 'T', double(T), start);
ncwrite(outfile, 'Xll', double(Xll'), [1 start]);
ncwrite(outfile, 'Yll', double(Yll'), [1 start]);
if test == 0 
	ncwrite(outfile, 'Xcart', double(Xcart'), [1 start]);
	ncwrite(outfile, 'Ycart', double(Ycart'), [1 start]);
end

if (pid), fprintf(1, 'outfile %s successfully written.\n', outfile); end;
gdone = 1;
return;
