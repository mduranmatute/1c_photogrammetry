%% #######################################################################%
%------------- Function used to determine the location of the ------------%
%-------------- camera lens from the magnification of images -------------% 
%#########################################################################%

% The camera lens position is determined by assuming that the inverse 
% image magnification is linearly related to the distance to the pinhole 
% lens (idealized camera), and will be 0 at the pinhole. Uses the
% subfunction [derivatives]. As the name indicates, it is used to calculate
% derivatives in the x and y directions.

function [Mx,My] = magnification_A(fitp)
    
    % A displacement dX on the grid leads to a displacement of du, dv on 
    % the ccd. Total displacement on the ccd : 
    % ds = sqrt(du^2 + dv^2),
    % where du = (du/dX)*dX and dv = (dv/dX)*dX
    % Mx = ds/dX, ds = sqrt(((du/dX)*dX)^2 + ((dv/dX)*dX)^2)
    % Mx = sqrt((du/dX)^2 + (dv/dX)^2)
    
    CCDmid = [800,600]; %mid location in the ccd (pixels)

    P0 = [eval_poly(fitp.cX, CCDmid(1), CCDmid(2)), ...
          eval_poly(fitp.cY, CCDmid(1), CCDmid(2))];

    [dxdX,dxdY] = derivatives(fitp.cx, P0(1),P0(2));
    [dydX,dydY] = derivatives(fitp.cy, P0(1),P0(2));
    Mx = sqrt(dxdX^2 + dydX^2);
    My = sqrt(dxdY^2 + dydY^2);
end

function [dudx, dudy] = derivatives(pf,x,y)
    % Analytic partial derivatives of the poly33 model
    % z = p00 + p10 x + p01 y + p20 x^2 + p11 xy + p02 y^2
    %        + p30 x^3 + p21 x^2 y + p12 x y^2 + p03 y^3
    dudx = pf.p10 + 2*pf.p20*x + pf.p11*y ...
         + 3*pf.p30*x.^2 + 2*pf.p21*x.*y + pf.p12*y.^2;
    dudy = pf.p01 + pf.p11*x + 2*pf.p02*y ...
         + pf.p21*x.^2 + 2*pf.p12*x.*y + 3*pf.p03*y.^2;
end
