%% #######################################################################%
%--------------------- Function used to identify dots --------------------%
%#########################################################################%
% The presets used for the identification of dots were defined for the 
% example cases treated here. Since the projected dots retain a circular
% shape, the MATLAB function [imfindcircles] is used for the location of 
% dots. However, in the case of high dot deformation, the function 
% [regionprops] can be substituted to identify shapes such as ellipses.

function centres = find_circles(I, range)
    
    pol = 'bright'; %looks for bright spots 
    [Pcen,r] = imfindcircles(I, range, 'ObjectPolarity', pol, ...
                             'Method','TwoStage', 'Sensitivity', 0.7, ...
                             'EdgeThreshold', 0.01);
    if ~isempty(Pcen)
        centres.x = Pcen(:,1);
        centres.y = Pcen(:,2);
    else
        centres.x = [];
        centres.y = [];
    end
end
