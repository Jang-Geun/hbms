
%%fn='~/tlippmann/hampton/ncks/hse_0001_his_vars_2days_spring_tide.nc'; % file name for model output (velocities)
%%fn='~/achilds/great_bay_modeling/data/gbe_0001_his_vars_test_30days_uvbar_untilt.nc'; % file name for model output (velocities)
%fn = '../data/gbe_0001_his_vars_30days_u_v_surface.nc'
fn = '/home/field/tlippmann/roms_his_hb202310.nc';

fnz = fn;
%%fnz='~/achilds/great_bay_modeling/data/gbe_0001_his_vars_test_30days_zeta.nc'; % file name for model output (velocities)

%%A = 1.0*10^(0);
%%A = 1.0;
A = 0.5;
%%A = 0.1;

dtmin = 0.5; %% in minutes
dtsecs = floor(dtmin*60);
%%basename = 'gb_meadnum_169_acres_53_npts_3079';
%%basename = 'gb_meadnum_52_acres_69_npts_5842';
%%basename = 'N_10000_circ_dover';
%%basename = 'N_1000_circ_site';
%%basename = 'N_50000_circ_site';
%%basename = 'N_50000_circ_bridge';
%%basename = 'N_50000_circ_center';
basename = 'N_50000_circ_dye';
basenameout = [basename '_' num2str(dtsecs) 'sec_A_' num2str(A, '%5.3f') ];
%fileindx = 1;
fprintf(1, 'fileindx = %d\n', fileindx);
sid = 1; %% 1] depth averaged  else] surface currents

%%did =  interp_Npts_func_tcl(fn, fnz, sid, A, dtmin, basename, basenameout, fileindx);
did =  interp_Npts_continuous_feed_func_tcl(fn, fnz, sid, A, dtmin, basename, basenameout, fileindx);

fprintf(1, 'did = %d\n', did);

