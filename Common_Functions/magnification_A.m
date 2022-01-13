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
    
    P0 = [fitp.cX(CCDmid), fitp.cY(CCDmid)];
    
    [dxdX,dxdY] = derivatives(fitp.cx, P0(1),P0(2));
    [dydX,dydY] = derivatives(fitp.cy, P0(1),P0(2));
    Mx = sqrt(dxdX^2 + dydX^2);
    My = sqrt(dxdY^2 + dydY^2);
end

function [dudx, dudy] = derivatives(c,x,y)
    names  = coeffnames(c);
    parval = coeffvalues(c);
    npar   = size(names,1);
    
    dudx = 0;
    dudy = 0;
    for i = 1:npar
        s = names{i};
        nx = s(2) - '0'; 
        ny = s(3) - '0';
        
        a = parval(i);
        if nx >= 1
            dudx = dudx + a * nx* x.^(nx-1) .* y.^ny;
        end
        if ny >= 1
            dudy = dudy + a * ny* x.^nx .* y.^(ny-1);
        end
    end
end
