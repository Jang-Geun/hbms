%%function [did] = interp_Npts_drifter_func_tcl(fn,fnz,sid,A,ntimes,svd,dtmin,basename,basenameout,fileindx,fcheck, lat_o, lon_o)
%% A = 1;  %%%1.06e-9;   %this is the dispersion
%% ntimes = input('How many time steps do you want? (e.g. 50)   ');
%% t_change = input('How much time between time steps? 1] 15 mins 10] 1.5 mins    ');
%%svd = input('save data to file?  1=yes, 0=no   ');
%%if (svd),
%%   outfile = [ input('enter file name (e.g. gb_traj_N_times_A ):    ', 's')  '.mat'];
%%end
%% sid = 1] depth avg 2] surface 
if (1)
      %sid = input('use 1] depth avg. 2] surface   ');
      sid = 1;
      if sid == 1,
       fn='../../../data/gbe_0001_his_vars_test_30days_uvbar_untilt.nc'; % file name for model output (velocities)
       else 
	fn = '../../../data/gbe_0001_his_vars_30days_u_v_untilt.nc';
       end 
	fnz='../../../data/gbe_0001_his_vars_test_30days_zeta.nc'; % file name for model output (velocities)
    %%A = 1.06e-9;   %this is the dispersion
    exponent = 0;
    A = 1.06*10^exponent;
    %ntimes = input('How many time steps do you want? (e.g. 50)   ');
    %dtmin = input('How much time between time steps? (minutes)   ');
    dtmin = 0.5;
    basename = [ input('Enter Location (ex: dover)   ', 's')    ];
    %%basename = ['N1kcirc_' basename];
    basename = ['N10kcirc_' basename];
    fprintf(1, 'using:  %=s\n', basename);
    %basenameout = [ input('Enter basenameout for outfile:    ', 's')    ];
    %basenameout = 'N1kcirc_20190809_Aminus0';
    basenameout = [basename '_' num2str(dtmin*60, '%3.3d') 's_Aminus' num2str(exponent)];
    %fileindx = input('How many files do you want? (1 through 8 or [1:8])?  ');
    fileindx = 1;
    %svd = input('save data to file?  1=yes, 0=no   ');
     svd = 0;
    %%fcheck = input('output file:  0] dont overwrite 1] overwrite  ');
    fcheck = 1;
    lat_o = 43.133652;
    lon_o = -70.902901;
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
h_bot = ncread(fn, 'h');
zeta = ncread(fnz, 'zeta');
[Nr, Er] = lltoutm(latr,lonr,23);
h_sea = 24-h_bot;
[No, Eo] = lltoutm(lat_o, lon_o, 23);
Ro = sqrt((Nr-No).^2+(Er-Eo).^2);
[Rmin, Rmini] = min(Ro);
riu = unique(Rmini);
nriu = length(riu);
rju = zeros(size(riu));
rm2 = zeros(size(riu));
for kk = 1:nriu 
	[rm2(kk),rju(kk)] = min(Ro(riu(kk),:));
end;
[rm2i, rjui] = min(rm2);
ri = riu(rjui);
rj = rju(rjui);
zeta_o = zeta(ri,rj);
%fig1 = plot(ot,reshape(zeta(ri,rj,:),length(ot),1));
%[x,y] = ginput(1)
ot_start = 7.0830e+15;
zeta_start = -22.7636;

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

R=6378.137e3;
um2deg=1/(R*pi/180);
vm2deg=1./(R*pi/180*cosd(latr));


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
[Xa, Ya,dnS,dnE] = read_points_func(basename,fileindx);
dnS = datenum(2019,8,31,6,0,0);
dnE = datenum(2019,9,2,8,0,0);
ntimes = floor((dnE-dnS)*24*60/(dtmin))+1;
dlim = inpolygon(Xa(:),Ya(:),lonlim2,latlim2);

difft  = mean(diff(t));
dt = dtmin *60;
istrt = floor((dnS-dnuvs)/difft*24*3600)+1;
%%angle = -53;
angle = get_angle;
ind2 = 0;

%dg = find(dlim == 1);
%X(1,:) = Xa(dg)';
%Y(1,:) = Ya(dg)';
X(1,:) = Xa';
Y(1,:) = Ya';
%N_pts = length(dg);
N_pts = length(Xa);

%M = matfile([outpath outfile], 'Writable', true);
wsize = 100;

TT = zeros(wsize+1,1);
DTr = delaunayTriangulation([lonr(:) latr(:)]);
Xi = my_interp_func(DTr, xr, X(1,:), Y(1,:), 0)';
Yi = my_interp_func(DTr, yr, X(1,:), Y(1,:), 0)';
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
DTh = delaunayTriangulation([xu(:) yu(:)]);
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

                df = find(yv>= miny-delt_y & yv<= maxy+delt_y & xv>= minx-delt_x & xv <= maxx+delt_x);
                DT = delaunayTriangulation([xr(df) yr(df)]);
                DTz(pp).DT = DT;
		dfz(pp).df = df;
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

DTr = delaunayTriangulation([lonr(:) latr(:)]);


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
            dn = find(isnan(um));
            if numel(dn),
                    um(dn) = 0;
            end
            dn = find(isnan(vm));
            if numel(dn),
                vm(dn) = 0;
            end
    	    
            %um=umr*cosd(angle)+vmr*sind(angle);
            %vm=vmr*cosd(angle)-umr*sind(angle);

            %um=um.*um2deg;
            %vm=vm.*vm2deg;
	    
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

				if (1),
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
           fprintf(1, 'i=%d \n', i);
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
	   %%M.Xll(ind1:ind2,1:N_pts) = Xll;
	   %%M.Yll(ind1:ind2, 1:N_pts) = Yll;
	   %%M.Xcart(ind1:ind2, 1:N_pts) = X(1:wsize,:);
	   %%M.Ycart(ind1:ind2, 1:N_pts) = Y(1:wsize,:);
	   %%M.T(ind1:ind2, 1) = T(1:wsize);
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
   %ind1 = ind2 + 1;
   %ind2 = ind2 + kk;
   start = ind2 + 1;
   Xll = my_interp_func(DTrr,lonr(:),X(1:kk,:),Y(1:kk,:),pid0);
   Yll = my_interp_func(DTrr,latr(:),X(1:kk,:),Y(1:kk,:),pid0);
   Xll = reshape(Xll,kk,N_pts);
   Yll = reshape(Yll,kk,N_pts);
   T = dnS+(TT(1:kk)-TTs)/(3600*24);
   gdone = write_netcdf_append_func(outfile, outpath, T(1:kk), Xll, Yll, X(1:kk,:), Y(1:kk,:), start, 0);
   if (~gdone), fprintf(1, 'did not write file properly. gdone = %d\n', gdone); end;
   %%M.Xll(ind1:ind2, 1:N_pts) = Xll;
   %%M.Yll(ind1:ind2, 1:N_pts) = Yll;
   %%M.Xcart(ind1:ind2, 1:N_pts) = X(1:kk,:);
   %%M.Ycart(ind1:ind2, 1:N_pts) = Y(1:kk,:);
   %%M.T(ind1:ind2, 1) = T(1:kk);
end

time_elapsed = toc


did = 1;
return


