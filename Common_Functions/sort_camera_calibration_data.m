%% #######################################################################%
%----------- Function used to relate the calibration plate dots-----------%
%------------------- to the photogrpahs of the plate  --------------------%
%#########################################################################%
% This function requires the following values, 
% Ipat : Images (photographs) of the calibration plate
% grid : Real world position (x,y) of the dots contructed as a virtual grid
% camdots : Pixel position of the detected dots from uploaded photographs

% Dots from the photographed plate are sorted by relating virtual grid to
% the dots found through interpolation. This requires user input to relate 
% four "orientation" dots to their corresponding photographed dots.

function [griddots, camdots] = sort_camera_calibration_data(Icam, grid, ...
                                                                   camdots)
    clf();
    %% ========================================= Ask user for a little help
    % ---------------------------------- Get four corresponding orientation  
    % ------------------------------- points in both grid and camera images
    s = 'Click on 4 dots, first in the left figure, then in the right one'; 

    % The relation of the four dots will cycle until the user is satisfied
    % with the selection.
    satisfied = false;
    %-------------------------------------------------- First for the grid
    while ~satisfied
        subplot(1,2,1);
        plot(grid.X, grid.Y, 'b*', grid.X(1), grid.Y(1), 'r*');
        axis ij;  % flip Y-axis
        title(s);

        subplot(1,2,2);
        imagesc(Icam);
        % colormap for 8 bit grayvalue images
        colormap(([1:255; 1:255; 1:255])'/256);     
        fprintf('\n%s\n', s);

        Rsearch = 40;
        for i = 1:4
            subplot(1,2,1);
            [x,y] = XYclicked(gcf);

            [r,ir] = min( sqrt((grid.X - x).^2 + (grid.Y - y).^2) );
            if r < Rsearch
                xpat(i) = grid.X(ir);
                ypat(i) = grid.Y(ir);
            else
                xpat(i) = x;
                ypat(i) = y;
            end
            hold on; plot(xpat(i), ypat(i), 'rs', 'markersize', 15);
        end
        
        %----------------------- Next, the photographed calibration plate
        for i = 1:4
            subplot(1,2,2);
            [x,y] = XYclicked(gcf);

            [r,ir] = min( sqrt((camdots.x - x).^2 + (camdots.y - y).^2) );
            if r < Rsearch
                xcam(i) = camdots.x(ir);
                ycam(i) = camdots.y(ir);
            else
                xcam(i) = x;
                ycam(i) = y;
            end
            hold on; plot(xcam(i), ycam(i), 'ys', 'markersize', 15);
        end

        %% The relation between the clicked points is used to find the rest
        cx = fit([xpat',ypat'], xcam', 'poly11');
        cy = fit([xpat',ypat'], ycam', 'poly11');
        % estimated position of all the dots in the camera image
        x = cx(grid.X, grid.Y);   
        y = cy(grid.X, grid.Y);
        % Plot estimated position
        hold on; plot(x,y,'gx'); hold off
        
        % Relate the photographed dots in the plate to the estimated
        % position if they are within a given -small- radius.
        Rsearch = 50;
        % clear k
        k = nan*ones(numel(x),1);
        for i = 1:numel(x)
            [r,ir] = min( sqrt((camdots.x - x(i)).^2 + ...
                               (camdots.y - y(i)).^2) );
            if r < Rsearch
                k(i) = ir;
            end
        end
        k1 = find(~isnan(k));   % patdots indices of valid points
        k2 = k(~isnan(k));      % camdots indices of valid points
        
        hold on; 
        plot(camdots.x(k2), camdots.y(k2), 'rs', 'markersize',14);
        hold off;

        cx = fit([grid.X(k1), grid.Y(k1)], camdots.x(k2), 'poly33');
        cy = fit([grid.X(k1), grid.Y(k1)], camdots.y(k2), 'poly33');

        x = cx(grid.X, grid.Y);
        y = cy(grid.X, grid.Y);
        hold on; plot(x,y,'ys','markersize',10); hold off
        title('yellow squares are the predicted positions'); 
        
        % Once satisfied with the result exit. Either y and enter or just
        % Enter works.
        a = input('Are you satisfied with the result? <y> : ', 's');
        satisfied = isempty(a) || a(1) == 'y';
    end
    
    % place related values as result
    griddots.x = grid.X(k1);
    griddots.y = grid.Y(k1);
    camdots.x  = camdots.x(k2);
    camdots.y  = camdots.y(k2);
end
