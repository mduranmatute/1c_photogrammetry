%% #######################################################################%
%----------- Function used to fit low-order polynomial models ------------%
%#########################################################################%
% Ordinary least-squares replacement for the subset of the MATLAB Curve
% Fitting Toolbox's fit(X, Z, model) used in this codebase, so that
% calibration no longer requires that toolbox to be licensed. Supported
% models:
%   fit_poly(x, z, 'poly1')   z = p1*x + p2                (x is Nx1)
%   fit_poly(X, z, 'poly11')  z = p00 + p10*x + p01*y       (X is Nx2)
%   fit_poly(X, z, 'poly33')  z = p00 + p10*x + p01*y + p20*x^2 + p11*x*y
%                               + p02*y^2 + p30*x^3 + p21*x^2*y
%                               + p12*x*y^2 + p03*y^3       (X is Nx2)
%
% Returns pf, a struct with a 'model' field and the named coefficients
% (matching the field names Curve Fitting Toolbox cfit objects expose,
% e.g. p00, p10, ..., as already relied upon in
% NetCDF_Convert/Calib_C.m), so a saved fit can be read the same way as
% before. Evaluate a fit with [Common_Functions/eval_poly.m].

function pf = fit_poly(X, Z, model)
    Z = Z(:);
    switch model
        case 'poly1'
            x = X(:);
            A = [x, ones(size(x))];
            c = A \ Z;
            pf.model = 'poly1';
            pf.p1 = c(1);
            pf.p2 = c(2);
        case 'poly11'
            x = X(:,1); y = X(:,2);
            A = [ones(size(x)), x, y];
            c = A \ Z;
            pf.model = 'poly11';
            pf.p00 = c(1); pf.p10 = c(2); pf.p01 = c(3);
        case 'poly33'
            x = X(:,1); y = X(:,2);
            A = [ones(size(x)), x, y, x.^2, x.*y, y.^2, ...
                 x.^3, x.^2.*y, x.*y.^2, y.^3];
            c = A \ Z;
            pf.model = 'poly33';
            names = {'p00','p10','p01','p20','p11','p02','p30','p21','p12','p03'};
            for i = 1:numel(names)
                pf.(names{i}) = c(i);
            end
        otherwise
            error('fit_poly:unknownModel', 'Unsupported model "%s".', model);
    end
end
