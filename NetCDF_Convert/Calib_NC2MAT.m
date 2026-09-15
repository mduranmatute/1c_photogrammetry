% Reverse of [Calib_C.m]: reconstructs a working camera/projector
% calibration .mat file (as produced by Camera_Calibration_Template.m /
% Projector_Calibration_Template.m, and consumed by Surface_calculation_
% template.m and Examples/Measure_underwater.m) from an archived
% calibration NetCDF file. Intended for someone who only has the
% published NetCDF dataset and wants to run the measurement pipeline
% against it.
%
% The reconstructed cx0/cy0/cdx/cdy are plain structs compatible with
% [Common_Functions/eval_poly.m], matching the poly33 model documented
% in Calib_C.m:
% val(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 +
%            p30*x^3 + p21*x^2*y + p12*x*y^2 + p03*y^3;
%
% Note: a NetCDF file converted before Plens and zlevels were added to
% Calib_C.m does not contain them; Plens/zlevels are then reconstructed
% as empty ([]), since that information was never archived in the first
% place.

Base_Loc = fullfile('..', 'Examples', 'Data_Files', 'NETCDF', ...
                     'Calibration_Data');
SaveB_Loc = fullfile('..', 'Examples', 'Data_Files', 'MATLAB', ...
                      'Calibration_data');

files = dir(fullfile(Base_Loc, '*.nc'));
names = sort({files.name});
nz = numel(names);

coeff_names = {'p00','p10','p01','p20','p11','p02','p30','p21','p12','p03'};

for iz = 1:nz
    Filename = fullfile(Base_Loc, names{iz});
    ncinfo_data = ncinfo(Filename);
    var_names = {ncinfo_data.Variables.Name};

    cx0 = read_poly33(Filename, 'cx0', coeff_names);
    cy0 = read_poly33(Filename, 'cy0', coeff_names);
    cdx = read_poly33(Filename, 'cdx', coeff_names);
    cdy = read_poly33(Filename, 'cdy', coeff_names);

    if ismember('Plens', var_names)
        Plens = ncread(Filename, 'Plens')';
    else
        Plens = [];
        warning('Calib_NC2MAT:noPlens', ...
            '%s has no Plens variable (archived before this was added); Plens set to [].', ...
            names{iz});
    end

    if ismember('zlevels', var_names)
        zlevels = ncread(Filename, 'zlevels')';
    else
        zlevels = [];
        warning('Calib_NC2MAT:noZlevels', ...
            '%s has no zlevels variable (archived before this was added); zlevels set to [].', ...
            names{iz});
    end

    Name_sep = regexp(names{iz}, '\.', 'split');
    MatName = [Name_sep{1} '.mat'];
    save(fullfile(SaveB_Loc, MatName), 'cx0', 'cy0', 'cdx', 'cdy', ...
        'Plens', 'zlevels');
end

function pf = read_poly33(Filename, varname, coeff_names)
    coeffs = ncread(Filename, varname);
    pf.model = 'poly33';
    for i = 1:numel(coeff_names)
        pf.(coeff_names{i}) = coeffs(i);
    end
end
