
fn='~/tlippmann/hampton/ncks/hse_0001_his_vars_2days_spring_tide.nc'; % file name for model output (velocities)
%%fn='~/achilds/great_bay_modeling/data/gbe_0001_his_vars_test_30days_uvbar_untilt.nc'; % file name for model output (velocities)
%fn = '../data/gbe_0001_his_vars_30days_u_v_surface.nc'

fnz = fn;
%%fnz='~/achilds/great_bay_modeling/data/gbe_0001_his_vars_test_30days_zeta.nc'; % file name for model output (velocities)

%%A = 1.0*10^(0);
A = 1.0;

dtmin = 0.5; %% in minutes
dtsecs = floor(dtmin*60);
basename = 'N_1000_circ_site';
basenameout = [basename '_' num2str(dtsecs) 'sec_A_' num2str(A, '%5.3f') ];
%%fileindx = 1;
fprintf(1, 'fileindx = %d\n', fileindx);
sid = 1; %% 1] depth averaged  else] surface currents

Tf = 6;   %% time in hours for continuous feed, 0=all at t=0 (no feed)

did =  interp_Npts_master_func(fn, fnz, sid, A, dtmin, basename, basenameout, fileindx, Tf);

fprintf(1, 'did = %d\n', did);

