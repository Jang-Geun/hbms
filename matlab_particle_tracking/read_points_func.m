
function [Xa,Ya,dnS,dnE] = read_points_func(basename,fileindx)

nf = length(fileindx);
Xa = [];
Ya = [];

for n = 1:nf
    nfile = fileindx(n);
    path = ['../pts_code/' ];
    fname = [path basename '_' num2str(nfile) '.mat'];
    load(fname);
    Xa = [Xa(:);X(:)];
    Ya = [Ya(:);Y(:)];
end
if ~exist('dnS'),
    dnS = [];
    dnE = [];
end

return


