% Reverse of [Patterns_C.m]: reconstructs a .mat file containing
% DotsToProject (the projected dot pattern images used as Pattern_Fl /
% Pattern_images_file input to the calibration and measurement scripts)
% from an archived pattern NetCDF file. Intended for someone who only
% has the published NetCDF dataset and wants to run the pipeline
% against it.

Base_Loc = fullfile('..', 'Examples', 'Data_Files', 'NETCDF', 'Patterns');
SaveB_Loc = fullfile('..', 'Examples', 'Data_Files', 'MATLAB', 'Patterns');

files = dir(fullfile(Base_Loc, '*.nc'));
names = sort({files.name});
nz = numel(names);

for iz = 1:nz
    Filename = fullfile(Base_Loc, names{iz});

    % Patterns_C.m stores DotsToProject (a logical array) as numeric
    % NC_Byte values; cast back to logical to match the original.
    DotsToProject = logical(ncread(Filename, 'Pattern_Images'));

    Name_sep = regexp(names{iz}, '\.', 'split');
    MatName = [Name_sep{1} '.mat'];
    save(fullfile(SaveB_Loc, MatName), 'DotsToProject');
end
