% Reverse of [StepsM_C.m]: reconstructs the intermediate dot-position
% .mat caches (DotsInCam_*, Projected_DotLocation_case*, SortedDots_*)
% from their archived NetCDF files, so a measurement can be recomputed
% (via Examples/Measure_underwater.m or Surface_calculation_template.m)
% using only these intermediate results and a finished calibration and
% pattern file, without re-running dot detection from the raw
% photographs.
%
% Each pattern's row in X_pos/Y_pos is padded with NaN out to the widest
% row in the file (Max_Ind). NaN can also be a genuine "this dot was not
% matched" marker within a row (in SortedDots_*.nc; see
% Common_Functions/sort_projected_dots.m), so a row cannot simply have
% all its NaN entries removed - only the padding (trailing NaN beyond
% the last real value in that row) is trimmed; any NaN before that
% point is a real, meaningful "unmatched" marker and is kept in place,
% since its position aligns a photographed dot with its corresponding
% projected dot elsewhere in the pipeline.
%
% Note: DotsInCam_*.nc does not archive the original photograph file
% paths (image_names) - only the detected dot positions (camdots) -
% since StepsM_C.m never wrote them to NetCDF. camdots alone is
% sufficient to recompute a measurement once the corresponding
% SortedDots_*.nc file has also been reconstructed; image_names is only
% needed if the interactive sorting step itself must be redone.

Base_Loc = fullfile('..', 'Examples', 'Data_Files', 'NETCDF', ...
                     'Data_process_steps', 'Measurements');
SaveB_Loc = fullfile('..', 'Examples', 'Data_Files', 'MATLAB', ...
                      'Data_process_steps', 'Measurements');

files = dir(fullfile(Base_Loc, '*.nc'));
names = sort({files.name});
nz = numel(names);

for iz = 1:nz
    Filename = fullfile(Base_Loc, names{iz});
    Name_sep = regexp(names{iz}, '\.', 'split');
    Step_Name = regexp(Name_sep{1}, '\_', 'split');
    MatName = [Name_sep{1} '.mat'];

    Xpos = ncread(Filename, 'X_pos');
    Ypos = ncread(Filename, 'Y_pos');
    NumF = size(Xpos, 1);

    dots = struct('x', {}, 'y', {});
    for il = 1:NumF
        last_valid = find(~isnan(Xpos(il,:)), 1, 'last');
        if isempty(last_valid)
            last_valid = 0;
        end
        dots(il).x = Xpos(il,1:last_valid)';
        dots(il).y = Ypos(il,1:last_valid)';
    end

    if strcmp(Step_Name{1}, 'SortedDots')
        image_dots = dots; %#ok<NASGU>
        save(fullfile(SaveB_Loc, MatName), 'image_dots');
    elseif strcmp(Step_Name{1}, 'Projected')
        patdots = dots; %#ok<NASGU>
        save(fullfile(SaveB_Loc, MatName), 'patdots');
    elseif strcmp(Step_Name{1}, 'DotsInCam')
        camdots = dots; %#ok<NASGU>
        save(fullfile(SaveB_Loc, MatName), 'camdots');
    end
end
