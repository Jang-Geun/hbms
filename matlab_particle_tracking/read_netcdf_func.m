function [T, Xll, Yll, A] = read_netcdf_func(fname, tind1, tcount, tstride, pind1, pcount, pstride)

T = ncread(fname, 'T', tind1, tcount, tstride);
Xll = ncread(fname, 'Xll', [pind1 tind1], [pcount tcount], [pstride tstride]);
Yll = ncread(fname, 'Yll', [pind1 tind1], [pcount tcount], [pstride tstride]);
A = ncread(fname, 'A');

