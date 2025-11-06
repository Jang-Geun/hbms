function [did] = interp_Npts_func_tcl(fn, fnz, sid, A, dtmin, basename, basenameout, fileindx)
%% A = 1;  %%%1.06e-9;   %this is the dispersion
%% t_change = input('How much time between time steps? 1] 15 mins 10] 1.5 mins    ');
%% sid = 1] depth avg 2] surface 
if (0)
      	%sid = input('use 1] depth avg. 2] surface   ');
      	sid = 1;
      	if sid == 1,
       		fn='~/achilds/great_bay_modeling/data/gbe_0001_his_vars_test_30days_uvbar_untilt.nc'; % file name for model output (velocities)
      	else 
		fn = '~/achilds/great_bay_modeling/data/gbe_0001_his_vars_30days_u_v_untilt.nc';
      	end 
	fnz='~/achilds/great_bay_modeling/data/gbe_0001_his_vars_test_30days_zeta.nc'; % file name for model output (velocities)
    	A = 1.0;  %% dispersion m^2/s
    	%exponent = 0;
    	%A = 1.06*10^exponent;
    	%dtmin = input('How much time between time steps? (minutes)   ');
    	dtmin = 0.5;
	dtsecs = floor(dtmin*60);
    	%basename = [ input('Enter Location (ex: dover)   ', 's')    ];
    	basename = 'gb_meadnum_169_acres_53_npts_3079';
    	%%basename = ['N10kcirc_' basename];
    	fprintf(1, 'using:  %s\n', basename);
    	%basenameout = [ input('Enter basenameout for outfile:    ', 's')    ];
    	%basenameout = 'N1kcirc_20190809_Aminus0';
    	%basenameout = [basename '_' num2str(dtmin*60, '%3.3d') 's_Aminus' num2str(exponent)];
	basenameout = [basename '_' num2str(dtsecs) 'sec_A_' num2str(A, '%5.3f') ];

    	%fileindx = input('How many files do you want? (1 through 8 or [1:8])?  ');
    	fileindx = 1;
end

outfile = [ basenameout '_' num2str(fileindx)  '_traj.nc'];
outpath = '../image_code/';

did = 0;
rid = 0;

t=ncread(fn,'ocean_time');
otu = ncreadatt(fn, 'ocean_time', 'units');
dnuvs = datenum(str2num(otu(15:18)), str2num(otu(20:21)), str2num(otu(23:24)), str2num(otu(26:27)), str2num(otu(29:30)), str2num(otu(32:33)));
lonr=ncread(fn,'lon_rho');
latr=ncread(fn,'lat_rho');
pm = ncread(fn,'pm');
pn = ncread(fn,'pn');
h_sea = 24-ncread(fn, 'h');
%%zeta = ncread(fnz, 'zeta');
%%[Nr, Er] = lltoutm(latr,lonr,23);
%%h_sea = 24-h_bot;

xr = lonr*0;
yr = latr*0;
for ii = 1:size(xr,1)-1
	xr(ii+1,:) = xr(ii,:) + 2 ./ (pm(ii+1,:)+pm(ii,:));
end
for jj = 1:size(yr,2)-1
	yr(:,jj+1) = yr(:,jj) + 2 ./ (pn(:,jj+1)+pn(:,jj));
end
xu = convn(xr,[.5;.5],'valid');
yu = convn(yr,[.5;.5],'valid');
xv = convn(xr,[.5 .5],'valid');
yv = convn(yr,[.5 .5],'valid');


lonlim = [lonr(1,1) lonr(1,end) lonr(end,end) lonr(end,1) lonr(1,1)];
latlim = [latr(1,1) latr(1,end) latr(end,end) latr(end,1) latr(1,1)];
lonlim2 = [lonr(2,2) lonr(2,end-1) lonr(end-1,end-1) lonr(end-1,2) lonr(2,2)];
latlim2 = [latr(2,2) latr(2,end-1) latr(end-1,end-1) latr(end-1,2) latr(2,2)];

xulim = [xu(2,2) xu(2,end-1) xu(end-1,end-1) xu(end-1,2) xu(2,2)];
yulim = [yu(2,2) yu(2,end-1) yu(end-1,end-1) yu(end-1,2) yu(2,2)];
xrlim = [xr(2,2) xr(2,end-1) xr(end-1,end-1) xr(end-1,2) xr(2,2)];
yrlim = [yr(2,2) yr(2,end-1) yr(end-1,end-1) yr(end-1,2) yr(2,2)];

if exist('X'), clear X; end
if exist('Y'), clear Y; end
%%[Xa, Ya,dnS,dnE] = read_points_func(basename,fileindx);
[X, Y,dnS,dnE] = read_points_func(basename,fileindx);
dnS = datenum(2019,8,31,6,0,0);
dnE = datenum(2019,9,2,8,0,0);
ntimes = floor((dnE-dnS)*24*60/(dtmin))+1;
fprintf(1, 'ntimes = %d\n', ntimes);
%%dlim = inpolygon(Xa(:),Ya(:),lonlim2,latlim2);
dlim = inpolygon(X(:),Y(:),lonlim2,latlim2);

difft  = mean(diff(t));
dt = dtmin *60;
istrt = floor((dnS-dnuvs)/difft*24*3600)+1;
%%angle = -53;
angle = get_angle;
ind2 = 0;

%%%dg = find(dlim == 1);
%%%X(1,:) = Xa(dg)';
%%%Y(1,:) = Ya(dg)';
%X(1,:) = Xa';
%Y(1,:) = Ya';
X = X(:)';
Y = Y(:)';
%%%N_pts = length(dg);
%N_pts = length(Xa);
N_pts = length(X);

wsize = 100;
%%wsize = 10;

TT = zeros(wsize+1,1);
DTrr = delaunayTriangulation([lonr(:) latr(:)]);
Xi = my_interp_func(DTrr, xr, X(1,:), Y(1,:), 0)';
Yi = my_interp_func(DTrr, yr, X(1,:), Y(1,:), 0)';
%%DTr = delaunayTriangulation([lonr(:) latr(:)]);
%%Xi = my_interp_func(DTr, xr, X(1,:), Y(1,:), 0)';
%%Yi = my_interp_func(DTr, yr, X(1,:), Y(1,:), 0)';
X = zeros(wsize+1,N_pts);
Y = zeros(wsize+1,N_pts);
X(1,:) = Xi;
Y(1,:) = Yi;
%%T(1)=t(istrt);        % initial time for particle releasing
%%T(1) = dnS;
TT(1) = floor((dnS-dnuvs)*24*3600); %number of seconds from begining 
TTs = TT(1);
%%dnt = dnuvs+t/3600/24;

%%following code gives delt_x and delt_y with the right scale
vmax = 10;
dmax = vmax*dt;
lon_o = -70.8191;
lat_o = 43.1133;
[y_o,x_o,zone] = lltoutm(lat_o,lon_o,23);
[lat_o2, lon_o2] = utmtoll(y_o+dmax,x_o+dmax,19,23);
%%delt_x = abs(lon_o-lon_o2);
%%delt_y = abs(lat_o-lat_o2);
delt_x = dmax;
delt_y = dmax;

%%Xmin = min(xu(:));
%%Xmax = max(xu(:));
%%Ymin = min(yu(:));
%%Ymax = max(yu(:));
Xmin = min(xr(:));
Xmax = max(xr(:));
Ymin = min(yr(:));
Ymax = max(yr(:));
nX = 10;
dXY = (Xmax-Xmin)/nX;
nY = floor((Ymax-Ymin)/dXY)+1;
%%DTh = delaunayTriangulation([xu(:) yu(:)]);
if (exist('dfu')), clear dfu; end;
if (exist('dfv')), clear dfv; end;
if (exist('DTu')), clear DTu; end;
if (exist('DTv')), clear DTv; end;
DTrr = delaunayTriangulation([xr(:) yr(:)]);
pp = 0;
for n = 1:nY
        miny = Ymin + (n-1)*dXY;
        maxy = miny+ dXY;
        for j = 1:nX
                pp = pp +1;
                minx = Xmin + (j-1)*dXY;
                maxx = minx + dXY;
                df = find(yu>= miny-delt_y & yu<= maxy+delt_y & xu>= minx-delt_x & xu <= maxx+delt_x);
                dfu(pp).df = df;
                DT = delaunayTriangulation([xu(df) yu(df)]);
                DTu(pp).DT = DT;

                %%df = find(yv>= miny-delt_y & yv<= maxy+delt_y & xv>= minx-delt_x & xv <= maxx+delt_x);
                df = find(yr>= miny-delt_y & yr<= maxy+delt_y & xr>= minx-delt_x & xr <= maxx+delt_x);
                DT = delaunayTriangulation([xr(df) yr(df)]);
                DTz(pp).DT = DT;
		dfz(pp).df = df;

                df = find(yv>= miny-delt_y & yv<= maxy+delt_y & xv>= minx-delt_x & xv <= maxx+delt_x);
                DT = delaunayTriangulation([xv(df) yv(df)]);
                DTv(pp).DT = DT;
                dfv(pp).df = df;
        end;
end;

k = 0; %counter for reading in velocities
kk = 0; %counter for writing out trajectories
secmax = 1*3600; % 1 day at a time
readpts = floor(secmax/difft);
tic
pid0 = 0;

%%DTr = delaunayTriangulation([lonr(:) latr(:)]);


for i = 1:ntimes
	kk = kk + 1;
	%if mod(i,100) == 1, 
        %    fprintf(1, 'i=%d \n', i);
	%end
	m = mod(i*dt, secmax);
        if (m == dt)
            k = k + 1;
            fprintf(1, 'k=%d reading data...', k);
            p = istrt + (k-1)*readpts;
	    addpts = 2;
            if sid == 1,
                um = ncread(fn, 'ubar', [1 1 p], [inf inf readpts+addpts], [1 1 1]);
                vm = ncread(fn, 'vbar', [1 1 p], [inf inf readpts+addpts], [1 1 1]);
	    else
                um = ncread(fn, 'u', [1 1 1 p], [inf inf inf readpts+addpts], [1 1 1 1]);
                vm = ncread(fn, 'v', [1 1 1 p], [inf inf inf readpts+addpts], [1 1 1 1]);
	    end
            zm = ncread(fnz,'zeta',[1 1 p], [inf inf readpts+addpts], [1 1 1]);
            %%ot = ncread(fn, 'ocean_time', p, readpts + addpts, 1);
	    fprintf(1, ' done reading data\n');
            zm = 24 +  zm;
	    um(find(isnan(um))) = 0;
	    vm(find(isnan(vm))) = 0;
            %dn = find(isnan(um));
            %if numel(dn),
                    %um(dn) = 0;
            %end
            %dn = find(isnan(vm));
            %if numel(dn),
                %vm(dn) = 0;
            %end
    	    
        end

        idt0=max(find((t-TT(kk))<=0));
        %%idt0=max(find((dnt-T(kk))<=0));
        idt = idt0 - p + 1;
        t0=t(idt0);
        u0=um(:,:,idt);
        v0=vm(:,:,idt);
        z0 =double(zm(:,:,idt));
        t1=t(idt0+1);
        u1=um(:,:,idt+1);
        v1=vm(:,:,idt+1);
        pp = 0;
        for n = 1:nY
                miny = Ymin + (n-1)*dXY;
                maxy = miny+ dXY;
                for j = 1:nX
                        pp = pp +1;
                        minx = Xmin + (j-1)*dXY;
                        maxx = minx + dXY;
                        dg0 = find(Y(kk,:)>= miny & Y(kk,:)< maxy & X(kk,:)>= minx & X(kk,:) < maxx);
			if numel(dg0),
				%% check to see if pts are within xulim, yulim
                        	%%dglim = inpolygon(X(kk,dg0),Y(kk,dg0),xulim,yulim);
                        	dglim = inpolygon(X(kk,dg0),Y(kk,dg0),xrlim,yrlim);
				dg1 = find(dglim == 1);
				if numel(dg1), 
					 %%fprintf(1, 'size dg = %d\n', size(dg1));
					dg = dg0(dg1);
				else
					dg = [];
				end;
			else
				dg = [];
			end;

                        if numel(dg),
                                %dfuu = dfu(pp).df;
                                %dfvv = dfv(pp).df;
				%dfzz = dfz(pp).df;
                                %DTzz = DTz(pp).DT;
                                %DTuu = DTu(pp).DT;
                                %DTvv = DTv(pp).DT;

                                u0i = my_interp_func(DTu(pp).DT,u0(dfu(pp).df),X(kk,dg),Y(kk,dg),pid0);
                                u1i = my_interp_func(DTu(pp).DT,u1(dfu(pp).df),X(kk,dg),Y(kk,dg),pid0);
                                ku1=((t1-TT(kk))*u0i+(TT(kk)-t0)*u1i)/(t1-t0);
                                v0i = my_interp_func(DTv(pp).DT,v0(dfv(pp).df),X(kk,dg),Y(kk,dg),pid0);
                                v1i = my_interp_func(DTv(pp).DT,v1(dfv(pp).df),X(kk,dg),Y(kk,dg),pid0);
                                kv1=((t1-TT(kk))*v0i+(TT(kk)-t0)*v1i)/(t1-t0);

                                x1=X(kk,dg)+ku1*dt/2;
                                y1=Y(kk,dg)+kv1*dt/2;

                                u0i = my_interp_func(DTu(pp).DT,u0(dfu(pp).df),x1,y1,pid0);
                                u1i = my_interp_func(DTu(pp).DT,u1(dfu(pp).df),x1,y1,pid0);
                                ku2=((t1-(TT(kk)+dt/2))*u0i+((TT(kk)+dt/2)-t0)*u1i)/(t1-t0);

                                v0i = my_interp_func(DTv(pp).DT,v0(dfv(pp).df),x1,y1,pid0);
                                v1i = my_interp_func(DTv(pp).DT,v1(dfv(pp).df),x1,y1,pid0);
                                kv2=((t1-(TT(kk)+dt/2))*v0i+((TT(kk)+dt/2)-t0)*v1i)/(t1-t0);
                                x2=X(kk,dg)+ku2*dt/2;
                                y2=Y(kk,dg)+kv2*dt/2;

                                u0i = my_interp_func(DTu(pp).DT,u0(dfu(pp).df),x2,y2,pid0);
                                u1i = my_interp_func(DTu(pp).DT,u1(dfu(pp).df),x2,y2,pid0);
                                ku3=((t1-(TT(kk)+dt/2))*u0i+((TT(kk)+dt/2)-t0)*u1i)/(t1-t0);

                                v0i = my_interp_func(DTv(pp).DT,v0(dfv(pp).df),x2,y2,pid0);
                                v1i = my_interp_func(DTv(pp).DT,v1(dfv(pp).df),x2,y2,pid0);
                                kv3=((t1-(TT(kk)+dt/2))*v0i+((TT(kk)+dt/2)-t0)*v1i)/(t1-t0);

                                x3=X(kk,dg)+ku3*dt/2;
                                y3=Y(kk,dg)+kv3*dt/2;

                                u0i = my_interp_func(DTu(pp).DT,u0(dfu(pp).df),x3,y3,pid0);
                                u1i = my_interp_func(DTu(pp).DT,u1(dfu(pp).df),x3,y3,pid0);
	                        ku4=((t1-(TT(kk)+dt))*u0i+((TT(kk)+dt)-t0)*u1i)/(t1-t0);

                                v0i = my_interp_func(DTv(pp).DT,v0(dfv(pp).df),x3,y3,pid0);
                                v1i = my_interp_func(DTv(pp).DT,v1(dfv(pp).df),x3,y3,pid0);
                                kv4=((t1-(TT(kk)+dt))*v0i+((TT(kk)+dt)-t0)*v1i)/(t1-t0);

                                C_mag = sqrt(u0i.^2+v0i.^2);
                                C_max = nanmax(C_mag);
                                if C_max > 0,
                                        A_array = A.*C_mag./C_max;
                                else
                                        A_array = 0*C_mag;
                                end


                                X(kk+1,dg)=X(kk,dg) + dt/6*(ku1+2*ku2+2*ku3+ku4) + normrnd(0, sqrt(2*dt/6.*A_array));
                                Y(kk+1,dg)=Y(kk,dg) + dt/6*(kv1+2*kv2+2*kv3+kv4) + normrnd(0, sqrt(2*dt/6.*A_array));
                                TT(kk+1)=TT(kk)+dt;
                                %%dlim = inpolygon(X(kk+1,dg),Y(kk+1,dg),xulim,yulim);
                                dlim = inpolygon(X(kk+1,dg),Y(kk+1,dg),xrlim,yrlim);
                                outp = find(dlim < 1);
                                if numel(outp),
		     		 	 Y(kk+1,dg(outp)) = Y(kk,dg(outp));
                                        X(kk+1,dg(outp)) = X(kk,dg(outp));
                                      if (rid),  fprintf(1, 'resetting %d points outside box \n', length(outp)); end
				      %if numel(outp) >= 10 
					% return
				 	%figure 
					%plot(X(kk,dg(outp)),Y(kk,dg(outp)),'ro','Markersize',18)
					%pl = plot(xrlim,yrlim, 'k-')
				      %end
					
                                end

				if (0),
                                    z0i = my_interp_func(DTz(pp).DT,z0(dfz(pp).df),X(kk,dg),Y(kk,dg),pid0);
                                    %%h0i = my_interp_func(DTh,h_sea(:),X(kk+1,dg),Y(kk+1,dg),pid0);
                                    h0i = my_interp_func(DTrr,h_sea(:),X(kk+1,dg),Y(kk+1,dg),pid0);

                                    dz = find(h0i > z0i);
                                    if (numel(dz)),
                                        X(kk+1,dg(dz)) = X(kk,dg(dz));
                                        Y(kk+1,dg(dz)) = Y(kk,dg(dz));
                                        if (rid),    fprintf(1, 'resetting %d points \n',length(dz)); end
                                    end
				end
                        end
                end
        end
	if mod(i,wsize) == 0,
           fprintf(1, 'i=%d (of %d) \n', i, ntimes);
	   ind1 = i - wsize + 1;
  	   ind2 = i;
	   Xll = my_interp_func(DTrr,lonr(:),X(1:wsize,:),Y(1:wsize,:),pid0);
	   Yll = my_interp_func(DTrr,latr(:),X(1:wsize,:),Y(1:wsize,:),pid0);
	   Xll = reshape(Xll,wsize,N_pts);
	   Yll = reshape(Yll,wsize,N_pts);
	   start = ind1;
	   T = dnS+(TT(1:wsize)-TTs)/(3600*24);
	   gdone = write_netcdf_append_func(outfile, outpath, T(1:wsize), Xll, Yll, X(1:wsize,:), Y(1:wsize,:), start, 0);
	   if (~gdone), fprintf(1, 'did not write file properly. gdone = %d\n', gdone); end;
	   kk = 0;
	   TT(1) = TT(end);
	   X(1,:) = X(end,:);
	   Y(1,:) = Y(end,:);
	   if i == wsize

		nccreate([outpath outfile], 'A', 'Dimensions', {'nA', 1}, 'Datatype', 'single');
	        ncwrite([outpath outfile], 'A', single(A));
		ncwriteatt([outpath outfile], 'A', 'long_name', 'Dispersion Coeff');
		ncwriteatt([outpath outfile], 'A', 'units', 'm^2/s');
		ncwriteatt([outpath outfile],'/', 'vel_file',fn);	
	   end
	end
end

if (kk ~= 0 )
   ind1 = ind2 + 1;
   ind2 = ind2 + kk;
   start = ind2 + 1;
   Xll = my_interp_func(DTrr,lonr(:),X(1:kk,:),Y(1:kk,:),pid0);
   Yll = my_interp_func(DTrr,latr(:),X(1:kk,:),Y(1:kk,:),pid0);
   Xll = reshape(Xll,kk,N_pts);
   Yll = reshape(Yll,kk,N_pts);
   T = dnS+(TT(1:kk)-TTs)/(3600*24);
   gdone = write_netcdf_append_func(outfile, outpath, T(1:kk), Xll, Yll, X(1:kk,:), Y(1:kk,:), start, 0);
   if (~gdone), fprintf(1, 'did not write file properly. gdone = %d\n', gdone); end;
end

time_elapsed = toc


did = 1;
return


