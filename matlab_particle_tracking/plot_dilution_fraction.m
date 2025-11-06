
pid = 0;

lid = input('Location? 1] Near Site 2] West of Bridge   3]  Center Harbor  ');
%%lid = 1;
if lid == 1
        lat_o = 42.888324366;
        lon_o = -70.825482967;
	fname0 = 'N_50000_circ_site_30sec_A_1.000_cont_feed_all_tstride_10_traj.nc';
	%%fname0 = 'N_1000_circ_site_30sec_A_1.000_cont_feed_all_tstride_1_traj.nc';
end
if lid == 2
        lat_o = 42.8970592218;
        lon_o = -70.8167521264;
	fname0 = 'N_50000_circ_bridge_30sec_A_1.000_cont_feed_all_tstride_10_traj.nc';
end
if lid == 3
        lat_o = 42.8988680105;
        lon_o = -70.8210248385;
	fname0 = 'N_50000_circ_center_30sec_A_1.000_cont_feed_all_tstride_10_traj.nc';
end
enmm = [348.0 355 4746.5 4755.0]*1000;

[No, Eo]= lltoutm(lat_o, lon_o, 23);

fnamepath = '../image_code/';
fname = [fnamepath fname0];
movid = input('make gif movie?  0] no   1] yes  ');
if movid == 1,
	s1 = findstr(fname0, '.nc');
	vnamepath = '../gifmovies_concentration/';
	vname = [vnamepath fname0(1:s1) 'gif'];
end;

tsid = input('write out time series at start point to file?  0] no  1] yes   ');
if tsid == 1,
	s1 = findstr(fname0, '.nc');
	tsnamepath = '../dilution_time_series/';
	tsname = [tsnamepath fname0(1:s1-1) '_ts_start_pt.mat'];
	tsnameplot = [tsnamepath fname0(1:s1-1) '_ts_start_pt.jpg'];
end;
	


finfo = ncinfo(fname);
[ntimes, npts, ss] = finfo.Dimensions.Length;
[T1, X1, Y1,A1] = read_netcdf_func(fname, 1, 2, 1, 1, inf, 1);
deltaT = diff(T1)*24*60;  %% in minutes
Tend = T1(1) + (deltaT/24/60)*ntimes;

fprintf(1, 'start time:  %s\n', datestr(T1(1)));
fprintf(1, '  end time:  %s\n', datestr(Tend));
fprintf(1, '    deltaT:  %f minutes\n', deltaT);
Tinterval = input('enter plotting time inteval (minutes):  ');
tstride = round(Tinterval/deltaT);
%%tstride = 20; %%1;
fprintf(1, '   tstride:  %d\n', tstride);
saveinterval = input('save frame interval (integer number):  [0 = dont save, 1] every frame] ');
if (saveinterval)
	s1 = findstr(fname0, '.nc');
	pnamepath = ['../frames/'];
	pnamebase = [pnamepath fname0(1:s1-1)];
end;
	
pstride = 1;
pind1 = 1;
pcount = inf;
%%wsize = 100;
%%tcount = 100;
tcount = 49;
tend = 1 - tstride;
%%nloop = floor(ntimes/(wsize*tstride));
%%nlo = ntimes - nloop*wsize*tstride;
nloop = floor(ntimes/(tcount*tstride));
nlo = ntimes - nloop*tcount*tstride;

%%wid = 100;
%%ht = 100;
wid = 25;
ht = 25;
nw = (enmm(2)-enmm(1));
nh = (enmm(4)-enmm(3));

%%Ev = [enmm(1) + wid/2 : wid : enmm(2)];
%%Nv = [enmm(3) + ht/2 : ht : enmm(4)];
Ev = [enmm(1) : wid : enmm(2)];
Nv = [enmm(3) : ht : enmm(4)];
Eedges = [enmm(1)-wid/2 : wid : enmm(2)+wid/2];
Nedges = [enmm(3)-ht/2  : ht  : enmm(4)+ht/2];

Eo_edges = [Eo-wid/2 : wid : Eo+wid/2];
No_edges = [No-ht/2 : wid : No+ht/2];

[Eaa, Naa] = meshgrid(Ev, Nv);
[nE, nN] = size(Eaa);

DT = delaunayTriangulation([Eaa(:) Naa(:)]);
 
dq = 30;
Eq = [enmm(1) : dq: enmm(2)];
Nq = [enmm(3) : dq: enmm(4)];
[Eqa, Nqa] = meshgrid(Eq, Nq);


figid = 11;
imid = 1;
[did] = plot_hampton_image_func(figid, 1, imid);
hl = colorbar('peer', gca);
set(get(hl, 'ylabel'), 'String', 'Dilution Fraction');
set(hl, 'fontsize', 16);
set(gca, 'colorscale', 'log');
colormap('jet');

if movid == 1,
	axis(enmm/1000)
	fprintf(1, 'resize figure for gif movie.  Hit return when done.\n');
	pause;
end;

framecount = 0;
elapsedtime = 0;
if exist('pc'), clear pc; end;

Tdilution = zeros(tcount*nloop, 1);
Xdilution = zeros(tcount*nloop, 1);
ndil = 0;

for i=1:nloop

	tind1 = tend + tstride;
	tend = tind1 + tstride*(tcount-1);
	%%tend = tind1 + tstride*(wsize-1);
	[T, X, Y,A] = read_netcdf_func(fname, tind1, tcount, tstride, pind1, pcount, pstride);
	X = X';
	Y = Y';
	T = T(:);

	[N, E] = lltoutm(Y, X, 23);

 
	for kk= 1:tcount

		ndil = ndil + 1;
		Xdilution(ndil) = (histcounts2(E(kk,:), N(kk,:), Eo_edges, No_edges)')/npts;
		Tdilution(ndil) = T(kk);

		elapsedtime = (T(kk)-T1(1))*24;  %% in hours

		if (pid), fprintf(1, 'before histcounts\n'); end;
		C = histcounts2(E(kk,:), N(kk,:), Eedges, Nedges)';
		if (pid), fprintf(1, 'after histcounts\n'); end;

		if (pid), fprintf(1, 'before my_interp_func\n'); end;
		Cq = my_interp_func(DT, C(:), Eqa, Nqa, 0)/npts;
		Cq = reshape(Cq, size(Eqa));
		if (pid), fprintf(1, 'after my_interp_func\n'); end;

		df = find(Cq == 0);
		if numel(df), Cq(df) = NaN; end;

		if exist('pc'), delete(pc); end;
		pc = pcolor(Eqa/1000, Nqa/1000, Cq);		
		shading('flat');
		%%shading('interp'); %% doesnt work with colorscale log
		caxis([1e-4 1e-2]);
		title([datestr(T(kk)) '  Elapsed time: ' num2str(elapsedtime, '%6.2f') ' hrs'], 'fontsize', 24);
		xlabel('Easting (km)', 'fontsize', 18);
		ylabel('Northing (km)', 'fontsize', 18);
		set(gca, 'tickdir', 'out', 'fontsize', 14);

		if movid == 1,
			drawnow
			frame = getframe(figid);
			im = frame2im(frame);
			[imind, cm] = rgb2ind(im,256);
			if (i==1 & kk==1),
				imwrite(imind, cm, vname, 'gif', 'Loopcount', inf);
			else
				imwrite(imind, cm, vname, 'gif', 'WriteMode', 'append', 'DelayTime', 0.001);
			end;
		end;

		framecount = framecount + 1;
		if (saveinterval),
			if (mod(framecount, saveinterval)==1),
				pname = [pnamebase '_' num2str(elapsedtime, '%6.2f') 'hrs.jpg'];
				fprintf(1, 'writing frame to %s\n', pname);
				eval(['print -djpeg ' pname]);
			end;
		end;

		if (pid), 
			fprintf(1, 'hit return to go to next step\n'); 
			pause;
		else
			pause(0.01);
		end;
	end;
end;

figure(figid+1);
clf;
semilogy((Tdilution-Tdilution(1))*24, Xdilution, 'ko-');
%%datetick
xlabel('Elapsed Time (hrs)', 'fontsize', 18);
ylabel('Dilution Fraction', 'fontsize', 18);
title('Dilution Fraction from Start Point', 'fontsize', 24);
set(gca, 'tickdir', 'out', 'fontsize', 14);

%% write out dilution time series at start point
if (tsid),
	fprintf(1, 'saveing ts figure to %s\n', tsnameplot);
	eval(['print -djpeg ' tsnameplot]);
	fprintf(1, 'writing ts at start point to file: %s\n', tsname);
	eval(['save ' tsname ' Tdilution Xdilution']);
end;

				

