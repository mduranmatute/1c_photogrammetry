%% #######################################################################%
%--------- Function used to evaluate a fit_poly() polynomial model -------%
%#########################################################################%
% z = eval_poly(pf, x, y)
% Evaluate a polynomial fit produced by [Common_Functions/fit_poly.m] at
% the given points. For a 'poly1' fit pass only x (y is ignored); for
% 'poly11' and 'poly33' fits pass both x and y.

function z = eval_poly(pf, x, y)
    switch pf.model
        case 'poly1'
            z = pf.p1 .* x + pf.p2;
        case 'poly11'
            z = pf.p00 + pf.p10.*x + pf.p01.*y;
        case 'poly33'
            z = pf.p00 + pf.p10.*x + pf.p01.*y + pf.p20.*x.^2 + ...
                pf.p11.*x.*y + pf.p02.*y.^2 + pf.p30.*x.^3 + ...
                pf.p21.*x.^2.*y + pf.p12.*x.*y.^2 + pf.p03.*y.^3;
        otherwise
            error('eval_poly:unknownModel', 'Unsupported model "%s".', pf.model);
    end
end
