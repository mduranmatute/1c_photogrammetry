% Reverse of [Measurements_C.m]: reconstructs a .mat file containing
% Pmid, xq, yq, vq (as produced by Examples/Measure_underwater.m and
% Surface_calculation_template.m) from an archived measurement NetCDF
% file. Intended for someone who only has the published NetCDF dataset
% and wants to use the reconstructed surface directly.

Base_Loc = fullfile('..', 'Examples', 'Data_Files', 'NETCDF', 'Measurements');
SaveB_Loc = fullfile('..', 'Examples', 'Data_Files', 'MATLAB', 'Measurements');

files = dir(fullfile(Base_Loc, '*.nc'));
names = sort({files.name});
nz = numel(names);

for iz = 1:nz
    Filename = fullfile(Base_Loc, names{iz});

    Pmid = ncread(Filename, 'Measurements');
    xvals = ncread(Filename, 'X_int');
    yvals = ncread(Filename, 'Y_int');
    vq = ncread(Filename, 'Surface');
    [xq, yq] = meshgrid(xvals, yvals);

    Name_sep = regexp(names{iz}, '\.', 'split');
    MatName = [Name_sep{1} '.mat'];
    save(fullfile(SaveB_Loc, MatName), 'Pmid', 'xq', 'yq', 'vq');
end
