ymd = [input('Enter YearMonthDay of file:   ','s')];
%%file1 = [input('file name?  ','s')]; %dispersion 1
file1 = ['N1kcirc_' ymd '_Aminus0_1_traj'];
%%file2 = [input('file name?  ','s')]; %dispersion 2
file2 = file1;
%%pid = input('Show particle motion? 0] no 1] yes  ');
pid = 0;
%%plid = input('Plot drifter trajectories? 0] no 1] yes   ');
plid = 1;
path = ['../image_code/'];
file1 = [path file1 '.nc']; 
file2 = [path file2 '.nc'];
finfo = ncinfo(file1);
[nt, np, ss] = finfo.Dimensions.Length;
wsize = 100;
tcount = wsize; 
tstride = 1; %every 50 time steps  
tend = 1 - tstride; 
pind1 = 1; 
pcount = inf; 
pstride = 50; % every 20 points 
nloop = floor(nt/(wsize*tstride));
nlo = nt - nloop*wsize*tstride;

[did] = plot_greatbay_image_func(11,1);
ax = [340 365 4760 4787];
title(ymd);
for i=1:nloop
       tind1 = tend + tstride;
       tend = tind1 + tstride*(wsize-1);
       [T1, X1, Y1,A1] = read_netcdf_func(file1, tind1, tcount, tstride, pind1, pcount, pstride);
       X1 = X1';
       Y1 = Y1';
       if i == 1
          [year1,month1,day1,~,~,~] = datevec(T1(1));
       end
       [T2, X2, Y2,A2] = read_netcdf_func(file2, tind1, tcount, tstride, pind1, pcount, pstride);
       X2 = X2';
       Y2 = Y2';
       [N1, E1, zone] = lltoutm(Y1, X1, 23);
       [N2, E2, zone] = lltoutm(Y2, X2, 23);
       pl1 = plot(E1/1000, N1/1000,'b-');
       pl2 = plot(E2/1000, N2/1000,'b--');
end
if nlo > 0 
    tind1 = tend + tstride;
    [T1, X1, Y1,A1] = read_netcdf_func(file1, tind1, nlo, tstride, pind1, pcount, pstride);
    [T2, X2, Y2,A2] = read_netcdf_func(file2, tind1, nlo, tstride, pind1, pcount, pstride);
     X1 = X1';
     Y1 = Y1';
     X2 = X2';
     Y2 = Y2';
     [N1, E1, zone] = lltoutm(Y1, X1, 23);
     [N2, E2, zone] = lltoutm(Y2, X2, 23);
     pl1 = plot(E1/1000, N1/1000,'b-');
     pl2 = plot(E2/1000, N2/1000,'b--');
end
if plid 
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

   pl1 = plot(Ed1/1000, Nd1/1000, 'w');
   pl1g = plot(Ed1(1)/1000, Nd1(1)/1000, 'go', 'markerfacecolor', 'g');
   pl1r = plot(Ed1(end)/1000, Nd1(end)/1000, 'ro', 'markerfacecolor', 'r'); 
   pl2 = plot(Ed2/1000, Nd2/1000, 'w');
   pl2g = plot(Ed2(1)/1000, Nd2(1)/1000, 'go', 'markerfacecolor', 'g');
   pl2r = plot(Ed2(end)/1000, Nd2(end)/1000, 'ro', 'markerfacecolor', 'r'); 
   pl3 = plot(Ed3/1000, Nd3/1000, 'w');
   pl3g = plot(Ed3(1)/1000, Nd3(1)/1000, 'go', 'markerfacecolor', 'g');
   pl3r = plot(Ed3(end)/1000, Nd3(end)/1000, 'ro', 'markerfacecolor', 'r'); 
end 

   
