N_input = input('How many points do you want?  ');
%%N_input = 1000;
nfiles = input('how many output files:   ');
%%nfiles = 1;
%%basename = [input('Enter basename for outfile:   ', 's')];
prid = 1;
R_o = input('Radius Size?  ');

lid = input('Location? 1] Near Site  2] West of Bridge  3] Center harbor  4] dye pump 20231031  ');
%%lid = 1;
if lid == 1
	basename = ['N_' num2str(N_input) '_circ_site'];
	lat_o = 42.888324366;
	lon_o = -70.825482967;
end
if lid == 2
	basename = ['N_' num2str(N_input) '_circ_bridge'];
	lat_o = 42.8970592218;
	lon_o = -70.8167521264;
end
if lid == 3
	basename = ['N_' num2str(N_input) '_circ_center'];
	lat_o = 42.8988680105;
	lon_o = -70.8210248385;
end
if lid == 4
	basename = ['N_' num2str(N_input) '_circ_dye'];
	lat_o = 42 + (53 + 11.9/60.)/60.
	lon_o = -(70 + (49 + 27.6/60.)/60.)
end

[No, Eo]= lltoutm(lat_o, lon_o, 23);

tavg = 0;
dnS= [];
dnE = [];
flag_lim = 1;

angle = get_angle;
theta = pi*2*(rand(1,N_input)-.5); %%100 is "n" number of particles in cloud
r = R_o*rand(1,N_input);
E= r.*cos(theta)+Eo;
N=r.*sin(theta)+No;
path = ['../pts_code/'];
pts = floor(N_input/nfiles);
N_pts = N_input;
for n = 1:nfiles
      ind1 = 1 + (n-1)*pts;
      ind2 = n*pts;
      %%X = E(ind1:ind2);
      %%Y = N(ind1:ind2);
      [Y,X] = utmtoll(N(ind1:ind2),E(ind1:ind2),19,23);

      fname = [path basename '_' num2str(n) '.mat'];
	fprintf(1, 'saving file to:  %s\n', fname);
      save(fname, 'X', 'Y','N_pts','dnS', 'dnE');
end



