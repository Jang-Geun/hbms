input1 = [ input('enter file name (e.g. gb_traj_N_times_A ):    ', 's') ];
%%fid = input('File type: 0] .nc 1] .m   ');
fid = 0;
if fid
   infile = [input1 '.mat'];
else 
   infile = [input1 '.nc'];
end
%%movid = input('movid as 0] dont export 1] mpeg-4 2] gif   ');
movid = 2;
if movid == 2, fprintf(1, 'exporting to gif movie\n'); end;
%%sid = input('use 1] depth avg. 2] surf  ');
sid = 1;
if sid == 1, fprintf(1, 'using depth-averaged currents\n'); end;
if sid == 2, fprintf(1, 'using surface currents\n'); end;

if fid
   load([inpath infile]);
else
    inpath = '../image_code/';
    fname = [inpath infile];
    finfo = ncinfo(fname);
    [nt,npts,ss] = finfo.Dimensions.Length;
    [T1, X1, Y1,A1] = read_netcdf_func(fname, 1, 2, 1, 1, inf, 1);

    ntimes = nt;
    %%wsize = 100;
    wsize = 49;
    tstride = 1;
    pstride = 1;
    pind1 = 1;
    pcount = inf;
    tcount = wsize;
    tend = 1 - tstride;
    nloop = floor(ntimes/(wsize*tstride));
    nlo = ntimes - nloop*wsize*tstride;
end
%this loads A, R_o, T, X, fn, Y 

%infn = input('read and plot velocities, too?  1=y  0=no   ');
infn = 0;
if infn == 0, fprintf(1, 'not plotting velocities\n'); end;
if infn == 1, fprintf(1, 'plotting velocities\n'); end;
if (infn),
	% file name for model output (velocities)
       if sid == 1,
          fn = '../../../data/gbe_0001_his_vars_test_30days_uvbar.nc';
       else
	  fn = '../../../data/gbe_0001_his_vars_30days_u_v_surface.nc'
       end
       scale = 5;
end;

%%yid = input('Which year to plot eelgrass?  0] dont plot   1] 2018  2] 2016  3] 2015  4] 1996  5] 1986   ');
yid=0;
if yid > 0
     llid = 0;
else 
     %%llid = input('Plot as 0] N and E  1] lat lon   ');
     llid = 0;
end;
%%dateid = input('Units of time read in  0] seconds 1] datenum   ');
dateid = 1;
%%plid = input('Plot drifter trajectories? 0] no 1] yes   ');
plid = 0;
if plid == 0, fprintf(1, 'not plotting drifters\n'); end;
if plid == 1, fprintf(1, 'plotting drifters\n'); end;
if fid 
   X = Xll';
   Y = Yll';
   ntimes = size(X,1);
   N_pts = size(X,2);
   scale = 3;
end
pid0 = 0;
if movid == 1,
   mpath = ['../mpegmovies/'];
   mname = [input1 '_mpegmovies.mp4'];
   mout = [mpath mname];
   v = VideoWriter(mout, 'MPEG-4');
   open(v);
end
if movid == 2,
   vpath = ['../gifmovies/'];
   vname = [input1 '_gifmovies.gif'];
   vout = [vpath vname];
end

%visualization
%plot great bay image
figid = 11;
if llid== 0
    %%[did] = plot_greatbay_image_func(figid,1);
	imid = 1;
    [did] = plot_hampton_image_func(figid, 1, imid);
    ax = [340 365 4760 4787];
    xlabel('Eastings (km)', 'fontsize', 18);
    ylabel('Northings (km)', 'fontsize', 18);
else
        [did] = plot_greatbay_image_func_latlon(figid,1);
      %  ax = [-70.91 -70.83 43.04 43.12];
end

fprintf(1,'hit return after resizing figure');
pause

elapsedtime = 0;

%%axis(ax)
hold on
%plot eel grass on great bay
        if yid > 0,
                [did] = plot_eelgrass_map_func(yid,llid);
        end;
if (infn)
	lonr=ncread(fn,'lon_rho');
	latr=ncread(fn,'lat_rho');
	DT = delaunayTriangulation([lonr(:) latr(:)]);
	[Nr, Er] = lltoutm(latr,lonr,23);
	angle = get_angle;
	R=6378.137e3;
	otu = ncreadatt(fn, 'ocean_time', 'units');
	dnuvs = datenum(str2num(otu(15:18)), str2num(otu(20:21)), str2num(otu(23:24)), str2num(otu(26:27)), str2num(otu(29:30)), str2num(otu(32:33)));
	t=ncread(fn,'ocean_time',1,2,1);
	dtsecs = t(2)-t(1);
end 
%pause %zoom in and press enter
pl = [];
while (1),

    if numel(pl), delete(pl); end
    for i=1:nloop
	       tind1 = tend + tstride;
               tend = tind1 + tstride*(wsize-1);
               [T, X, Y,A] = read_netcdf_func(fname, tind1, tcount, tstride, pind1, pcount, pstride);
               X = X';
               Y = Y';
               T = T(:); 
               if i == 1 & plid
		      [year1,month1,day1,~,~,~] = datevec(T(1));
                     sy = num2str(year1);
                     sm = num2str(month1, '%2.2d');
                     sd = num2str(day1, '%2.2d');
                     path = ['/home/achilds/drifters/data/' sy sm sd '/'];
                     prid = 1;
                     tavg = 0;
                     dnstart = [];
                     flag_lim = 1;
                     [dn1, latd1, lond1, Nd1, Ed1, submerged1, gps_flags1, sst1, batt1] = get_reef_func(1, prid, sy, sm, sd, tavg, dnstart, flag_lim);
                     [dn2, latd2, lond2, Nd2, Ed2, submerged2, gps_flags2, sst2, batt2] = get_reef_func(2, prid, sy, sm, sd, tavg, dnstart, flag_lim);
                     [dn3, latd3, lond3, Nd3, Ed3, submerged3, gps_flags3, sst3, batt3] = get_reef_func(3, prid, sy, sm, sd, tavg, dnstart, flag_lim);
                   
                     pl1 = plot(Ed1/1000, Nd1/1000, 'm');
                     pl1g = plot(Ed1(1)/1000, Nd1(1)/1000, 'go', 'markerfacecolor', 'g');
                     pl1r = plot(Ed1(end)/1000, Nd1(end)/1000, 'ro', 'markerfacecolor', 'r');
                     pl2 = plot(Ed2/1000, Nd2/1000, 'm');
                     pl2g = plot(Ed2(1)/1000, Nd2(1)/1000, 'go', 'markerfacecolor', 'g');
                     pl2r = plot(Ed2(end)/1000, Nd2(end)/1000, 'ro', 'markerfacecolor', 'r');
                     pl3 = plot(Ed3/1000, Nd3/1000, 'm');
                     pl3g = plot(Ed3(1)/1000, Nd3(1)/1000, 'go', 'markerfacecolor', 'g');
                     pl3r = plot(Ed3(end)/1000, Nd3(end)/1000, 'ro', 'markerfacecolor', 'r');
               end	
               if (infn),
	             p = floor((T(1)-dnuvs)*24*3600/dtsecs);
		     p2 = ceil((T(end)-dnuvs)*24*3600/dtsecs);
                     ntn = p2-p+2;
	             if sid == 1,
                        umr = ncread(fn, 'ubar_eastward', [1 1 p], [inf inf ntn], [1 1 1]);
                        vmr = ncread(fn, 'vbar_northward', [1 1 p], [inf inf ntn], [1 1 1]);
	             else
                        umr = ncread(fn, 'u_eastward', [1 1 1 p], [inf inf inf ntn], [1 1 1 1]);
                        vmr = ncread(fn, 'v_northward', [1 1 1 p], [inf inf inf ntn], [1 1 1 1]);
	             end
                     ot=ncread(fn,'ocean_time',p,ntn,1);
                     if dateid 
                          ot = ot./(24*3600)+dnuvs;
                     end
                     dn = find(isnan(umr));
                     if numel(dn),
                             umr(dn) = 0;
                     end
                     dn = find(isnan(vmr));
                     if numel(dn),
                             vmr(dn) = 0;
                     end
	             um2deg=1/(R*pi/180);
	             vm2deg=1./(R*pi/180*cosd(latr));
             
	             u=umr*cosd(angle)+vmr*sind(angle);
	             v=vmr*cosd(angle)-umr*sind(angle);
             
 	             u=u.*um2deg;
 	             v=v.*vm2deg;
               end
               for kk = 1:tcount
		    if (infn)
                    	tind1 = max(find(ot<= T(kk)));
                    	tind2 = min(find(ot>= T(kk)));
                    	u1 =  u(:,:,tind1);
                    	u2 =  u(:,:,tind2);
                    	um = ((ot(tind2)-T(kk))*u1+(T(kk)-ot(tind1))*u2)/(ot(tind2)-ot(tind1));
                    	v1 = v(:,:,tind1);
                    	v2 = v(:,:,tind2);
                    	vm = ((ot(tind2)-T(kk))*v1+(T(kk)-ot(tind1))*v2)/(ot(tind2)-ot(tind1));
                   end  
                    if ~(kk == 1 & i ==1)
		            if (infn), delete(qu);  end;
		            delete(pl); 
	            end
                    %%if (infn)
                     %%   [~,tid]=min(abs(t-T(k)));
	            %%end
                    %hold on
                    [N,E]=lltoutm(Y(kk,:),X(kk,:),23);
                    if llid == 0
	                if (infn)
            	            qu = quiver(Er/1000,Nr/1000, um, vm,scale,'w');
	                end;
                        pl=plot(E/1000,N/1000,'bo','MarkerFaceColor','b','MarkerSize',3);
                    else
	                if (infn)
            	            qu = quiver(lonr,latr, um, vm,scale,'w');
	                end;
                        pl=plot(X(k,:),Y(k,:),'bo','MarkerFaceColor','b','MarkerSize',3);
                    end;

			elapsedtime = (T(kk)-T1(1))*24;  %% in hours
	                title([datestr(T(kk)) '  Elapsed time: ' num2str(elapsedtime, '%6.2f') ' hrs'], 'fontsize', 24);

%pause(.25);
                  	drawnow
             	if movid == 1,
                  	frame = getframe(gcf);
                  	writeVideo(v, frame);
             	end;
             	if movid == 2,
                  	frame = getframe(figid);
                  	im = frame2im(frame);
                  	[imind,cm] = rgb2ind(im,256);
                  	if i == 1 && kk == 1
	                  	imwrite(imind,cm,vout,'gif', 'Loopcount', inf);
                  	else
	                  	imwrite(imind,cm,vout, 'gif', 'WriteMode', 'append', 'DelayTime', 0.001);
				%%fprintf(1, 'i = %d  kk = %d/n',i,kk);
                  	end
              	end
               end
     end
     %pause;
     inp= input('loop again? 0=no 1=yes  ');
     if (inp == 0 | numel(inp)==0)
          break
     else 
          tend = 1 - tstride;
     end
end

%pl=plot(X(:,10),Y(:,10),'r');

