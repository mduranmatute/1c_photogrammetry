%% #######################################################################%
%---------- Function used to determine the intersection of the -----------%
%----------- constructed virtual lines with a defined surface ------------%
%#########################################################################%
%[Ps, normal] = find_intersection_with_surface(S, r0, dv)
% Ps: position of the intersection
% normal: vector normal to the surface.

% S: Surface to be intersected. For the cases treated here, only a flat and
%    and parabolic surface are implemented.
% r0: Vector origin (x0,y0,z0)
% dv: Directional vector (dx,dy,dz)

% Ps = r0 + S * dv;
% where:
% Ps.x = x0 + S * dx; Ps.y = y0 + S * dy; Ps.z = z0 + S * dz;

function [Ps, normal] = find_intersection_with_surface(S, r0, dv)

    if isnan(S.a)       % no rotation (omega = 0), i.e. flat water surface
        Ps.x = r0(:,1) + dv(:,1) * S.z0;
        Ps.y = r0(:,2) + dv(:,2) * S.z0;
        Ps.z = r0(:,3) + dv(:,3) * S.z0;

        normal = repmat([0,0,1], numel(Ps.x), 1);
        
    else % for a parabolic surface (Omega ~= 0)

        x0 = r0(:,1);
        dx = dv(:,1);
        y0 = r0(:,2);
        dy = dv(:,2);

        z0 = S.z0;
        A  = S.a;  % the parabolic surface is now z = z0 + A*r

        % --- Determine intersection point ----------------------------
        % Equation of a line:
        % x = x0 + z*dx;
        % y = y0 + z*dy;
        % z = z;

        % Parabolic surface (axis == z-axis)
        % z = z0 + A*r^2
        % r = sqrt(x^2 + y^2)

        % z = z0 + A*( (x0 + z*dx)^2 + (y0 + z*dy)^2 );
        % Solving for z:
        %   z = (-b +/i sqrt(b^2 - 4*a*c)) /(2*a),
        % where
        %   a = A*(dx^2 + dy^2);
        %   b = A*(2*x0*dx + 2*y0*dy) - 1;
        %   c = z0 + A*(x0^2 + y0^2); 
        % ----------------------------------------------------

        a = A*(dx.^2 + dy.^2);
        b = A*(2*x0.*dx + 2*y0.*dy) - 1;
        c = z0 + A*(x0.^2 + y0.^2); 

        D = b.^2 - 4.*a.*c;
        % There are two solutions. It is necesary to select the best z. 
        % A criteria is that it must be within a given range.
        z = (-b + sqrt(D)) ./ (2*a);
        if z(1) < 0 || z(1) > 500
            z = (-b - sqrt(D)) ./ (2*a);
        end
        %r = sqrt((x0 + dx).^2 + (y0 + dy).^2);
        %z = S.z0 + S.a .* r.^2;
        Ps.x = x0 + z.*dx;
        Ps.y = y0 + z.*dy;
        Ps.z = z;

        % the normal vector in the point x,y,z
        % r = sqrt(x^2 + y^2);
        % z = z0 + A*r^2;
        % direction of a line tangent in x,y : dz/dr = 2*A*r
        % so the direction of the normal in x,y : -1/(dz/dr) = -1/(2*A*r)
        % The normal intersects the z-axis at deltaZ above z
        % deltaZ = -r * -1/(2*A*r) = 1/(2*A)
        % vector n = [-r, deltaZ] or n = [-x, -y, deltaZ]

        deltaZ = 1/(2*A);
        n = [-Ps.x, -Ps.y, ones(size(Ps.x))*deltaZ];
        %n = [-2*Ps.x*A,-2*A*Ps.y, ones(size(Ps.x))];
        
        a = sqrt( sum(n'.^2) );
        
        % n is then normalized
        %normal = n./a'; % this also works in Matlab >2018 versions
        for jk =1:3
            normal(:,jk) = n(:,jk)./a'; 
        end
        
    end
end
