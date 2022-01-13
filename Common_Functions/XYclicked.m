%% #######################################################################%
%-------------- Function used to get clicked pixel coordinates------------%
%#########################################################################%
function [x,y] = XYclicked(fig)
    H = datacursormode(fig);
    set(H,'SnapToDataVertex','off','Enable','on')
    waitforbuttonpress
    c = getCursorInfo(H);
    if ~isempty(c)
        x = c.Position(1);
        y = c.Position(2);
    else
        fprintf('Clicked outside the figure !!!\n');
        x = [];
        y = [];
    end
    
    set(H,'SnapToDataVertex','off','Enable','off')
end
