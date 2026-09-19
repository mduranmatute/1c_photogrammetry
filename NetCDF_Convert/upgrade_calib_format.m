% Losslessly converts calibration .mat files saved with the old Curve
% Fitting Toolbox cfit-based representation (cx0/cy0/cdx/cdy as cfit
% objects) into the new plain-struct representation used by
% Common_Functions/fit_poly.m and eval_poly.m. Extracts the
% already-computed poly33 coefficients (no re-fitting) and rewraps
% them; Plens and zlevels are carried over unchanged (or set to [] with
% a warning if a particular file never had them, e.g. an older
% calibration saved before zlevels was routinely included).
%
% Usage:
%   upgrade_calib_format(in_dir, out_dir)
%
% Run this once against a Calibration_data/ folder produced by the old
% (Curve Fitting Toolbox-based) Camera_Calibration_Template.m /
% Projector_Calibration_Template.m, before using that data with the
% current pipeline. Pass a different out_dir than in_dir to inspect the
% result first; pass the same directory to upgrade in place.

function upgrade_calib_format(in_dir, out_dir)
    coeff_names = {'p00','p10','p01','p20','p11','p02','p30','p21','p12','p03'};
    files = dir(fullfile(in_dir, '*.mat'));

    for i = 1:numel(files)
        D = load(fullfile(in_dir, files(i).name));
        out = struct();
        for fn = {'cx0','cy0','cdx','cdy'}
            old = D.(fn{1});
            pf = struct();
            pf.model = 'poly33';
            for c = 1:numel(coeff_names)
                pf.(coeff_names{c}) = old.(coeff_names{c});
            end
            out.(fn{1}) = pf;
        end

        if isfield(D, 'Plens')
            out.Plens = D.Plens;
        else
            out.Plens = [];
            warning('upgrade_calib_format:noPlens', ...
                '%s has no Plens field; set to [].', files(i).name);
        end

        if isfield(D, 'zlevels')
            out.zlevels = D.zlevels;
        else
            out.zlevels = [];
            warning('upgrade_calib_format:noZlevels', ...
                '%s has no zlevels field; set to [].', files(i).name);
        end

        save(fullfile(out_dir, files(i).name), '-struct', 'out');
    end
end
