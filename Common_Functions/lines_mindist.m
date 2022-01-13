%% #######################################################################%
%------------ Function used to determine the minimum distance ------------%
%------------- (or intersection) between two "virtual" lines -------------%
%#########################################################################%
% line.p0 = point on the line
% line.dv = direction vector of the line
% d = nearest distance between the lines
% P = mid-point off the connecting line at the nearest distance 

function [Pmid d] = lines_mindist(line1, line2)

    % Starting point of line P and Q, respectively.
    P0 = line1.p0;
    Q0 = line2.p0;
    
    % Also include and normailize the unit vector of line P and Q
    u  = line1.dv./norm(line1.dv);
    v  = line2.dv./norm(line2.dv);
    
    a = dot(u,u);
    b = dot(u,v);
    c = dot(v,v);
    w0 = P0 - Q0;
    d = dot(u,w0);
    e = dot(v,w0);
    
    % Calculate the parametric length of line P and Q, which are
    % s and t, respectively:
    s = (b*e - c*d) / (a*c - b^2);
    t = (a*e - b*d) / (a*c - b^2);
    
    % Location of Wc = Pc - Qc
    Pc = P0 + s*u; 
    Qc = Q0 + t*v;
    
    % Find the Mid-point of Wc (closest point between the two lines)
    Pmid = (Pc + Qc)/2;
    d = norm(Pc-Qc);
   
end
