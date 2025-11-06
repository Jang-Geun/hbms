#!/bin/bash


matlab -nodisplay -nodesktop -nojvm -r 'fileindx=1;' < interp_Npts_control_master.m >& ../run_output/test1.out&
matlab -nodisplay -nodesktop -nojvm -r 'fileindx=2;' < interp_Npts_control_master.m >& ../run_output/test2.out&
matlab -nodisplay -nodesktop -nojvm -r 'fileindx=3;' < interp_Npts_control_master.m >& ../run_output/test3.out&
matlab -nodisplay -nodesktop -nojvm -r 'fileindx=4;' < interp_Npts_control_master.m >& ../run_output/test4.out&
matlab -nodisplay -nodesktop -nojvm -r 'fileindx=5;' < interp_Npts_control_master.m >& ../run_output/test5.out&
matlab -nodisplay -nodesktop -nojvm -r 'fileindx=6;' < interp_Npts_control_master.m >& ../run_output/test6.out&
matlab -nodisplay -nodesktop -nojvm -r 'fileindx=7;' < interp_Npts_control_master.m >& ../run_output/test7.out&
matlab -nodisplay -nodesktop -nojvm -r 'fileindx=8;' < interp_Npts_control_master.m >& ../run_output/test8.out&
