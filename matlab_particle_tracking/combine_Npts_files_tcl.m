%basenameout = ['N160k2_A_minus0_30day2']; %change for each file 
basenameout = input('input basename:   [no path, without the _#.nc] ', 's');
path = ['../image_code/'];
fname1 = [path basenameout '_' '1_traj.nc'];
[T1, Xll1, Yll1, A1] = read_netcdf_func(fname1, 1, 2, 1, 1, inf, 1);	
deltaT = diff(T1)*24*60; %% in minutes
finfo = ncinfo(fname1);
[nt, np, ss] = finfo.Dimensions.Length;
pid = 0;
%%wsize = 100;
wsize = 49;
%%nf = 8; %chnage for each run 
aline = input('list of files to combine, e.g., 1 2 3 4 ..  ', 's'); 
[nfa, nf] = sscanf(aline, '%d');
Xcart = [];
Ycart = [];

%%tstride = 5; %read every time steps 
fprintf(1, 'Time interval is %f minutes.\n', deltaT);
Tinterval = input('Enter output time interval (minutes)  ');
tstride = round(Tinterval/deltaT);
fprintf(1, 'tstride = %d\n', tstride);

pstride = 1;
pind1 = 1;
pcount = inf;
nloop = floor(nt/(wsize*tstride));
nlo = nt - nloop*wsize*tstride;
tcount = wsize;
tend = 1 - tstride;
fnout = [basenameout '_all_tstride_' num2str(tstride) '_traj.nc'];
fprintf(1, 'output file is:  %s%s\n', path, fnout);

for k = 1:nloop
    tind1 = tend + tstride;
    tend = tind1 + tstride*(wsize-1);
    Xlla = [];
    Ylla = [];
    Ta = [];
    for m = 1:nf
	n = nfa(m);
	fname = [path basenameout '_' num2str(n) '_traj.nc'];
        if k ==1 ,fprintf(1, 'loading %s \n', fname);  end
	[T, Xll, Yll, A] = read_netcdf_func(fname, tind1, tcount, tstride, pind1, pcount, pstride);	
	Xlla = [Xlla Xll'];
        Ylla = [Ylla Yll'];
	if n == 1
	   Ta = [Ta; T(:)];
	end
    end
   T = Ta;
   Xll = Xlla;
   Yll = Ylla;
   start = 1+(k-1)*wsize;
   [gdone] = write_netcdf_append_func(fnout, path, T, Xll, Yll, Xcart, Ycart, start, pid);  
   if k == 1
	nccreate([path fnout], 'A', 'Dimensions', {'nA', 1}, 'Datatype', 'single');
        ncwrite([path fnout], 'A', single(A));
        ncwriteatt([path fnout], 'A', 'long_name', 'Dispersion Coeff');
        ncwriteatt([path fnout], 'A', 'units', 'm^2/s');
   end

end

tcount = floor(nlo/tstride);
if tcount > 0 
	tind1 = tend + tstride;
	 Xlla = [];
         Ylla = [];
         Ta = [];

	for n = 1:nf 
		[T, Xll, Yll, A] = read_netcdf_func(fname, tind1, tcount, tstride, pind1, pcount, pstride);
		Xlla = [Xlla Xll'];
		Ylla = [Ylla Yll'];	
		if n == 1
			Ta = [Ta; T(:)];
		end
		T = Ta;
		Xll = Xlla;
		Yll = Ylla;
	end
	start = 1+(k)*wsize;
        [gdone] = write_netcdf_append_func(fnout, path, T, Xll, Yll, Xcart, Ycart, start, pid);
	
end


