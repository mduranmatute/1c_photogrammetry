function [pvo,isvalid] = snell(pvi, normal, refr1, refr2)
    % isvalid it set to true if the refracted ray exist.
    isvalid = false;
    pvo = 0;

    pvi = pvi/norm(pvi);
    n   = normal/norm(normal);	% must be a unit vector
    u   = pvi;                  % direction incedent ray

    % the pointing vector after refraction is in the plane n,u
    % v     = p*n + q*u
    % v x n / u x n = sin(theta2)/sin(theta1) = refr1/refr2
    % ie: (p*n + q*u) x n = refr1/refr2 * u x n -->
    q = refr1/refr2;

    % |v| = 1 --> (p*n + q*u).(p*n + q*u) = 1
    % p^2 + + 2pq(n.u) + q^2 - 1 = 0
    b = 2*q * dot(n,u);
    c = q*q - 1;
    D = b*b - 4*c;
    if D < 0                    % no solution
        return;
    end
    p1 = (-b + sqrt(D))/2;
    p2 = (-b - sqrt(D))/2;
    if (abs(p1) < abs(p2))
        pvo = n*p1 + u*q;
    else
        pvo = n*p2 + u*q;
    end
    isvalid = true;
end
