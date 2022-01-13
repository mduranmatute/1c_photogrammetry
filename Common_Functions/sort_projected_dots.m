%% #######################################################################%
%----------- Function used to relate detected dots from photos -----------%
%------------------ to the projected dot pattern file  -------------------%
%#########################################################################%
% This function requires the following values, 
% Ipat : Dot Patterns to be projected 
% Icam : Images of the projected patterns captured by the camera
% patdots : Pixel position of the dots in the projected pattern file
% camdots : Pixel position of the detected dots from uploaded photographs

% Dots from photographs are sorted by relating the projected patterns to
% the camera dots found through interpolation. This requires user input to
% relate four "orientation" projector dots to their corresponding 
% photographed dots.

function camdots = sort_projected_dots(Ipat, Icam, patdots, camdots)
    %% ========================================= Ask user for a little help
    % --------------------------- Get four corresponding orientation points 
    % ----------------------------------- in both pattern and camera images
    s = 'Click on 4 dots, first in the top figure, then in the bottom one'; 
    
    % The relation of the four dots will cycle until the user is happy with
    % the selection.
    satisfied = false;
    clf();
    %------------------------------------ First for the projected patterns
    while ~satisfied
        subplot(2,1,1); imagesc(Ipat{1});
        title('pattern dots');
        subplot(2,1,2); imagesc(Icam{1}); 
        colormap(([1:255; 1:255; 1:255])'/256);
        title('image dots');
        fprintf('\n%s\n', s);

        patd = patdots(1);
        camd = camdots(1);
        Rsearch = 40;
        for i = 1:4
            subplot(2,1,1);
            found = 0;
            while ~found
                [x,y] = XYclicked(gcf);
                [r,ir] = min(sqrt((patd.x - x).^2 + (patd.y - y).^2));
                if r < Rsearch
                    xpat(i) = patd.x(ir);
                    ypat(i) = patd.y(ir);
                    found = 1;
                else
                    beep();
                    fprintf('Bad click, try again\n');
                end
            end
            hold on; plot(xpat(i), ypat(i), 'rs', 'markersize', 15);
        end

        %----------------------------------- Next, the photorpahed pattern 
        for i = 1:4
            subplot(2,1,2);
            found = 0;
            while ~found
                [x,y] = XYclicked(gcf);
                [r,ir] = min( sqrt((camd.x - x).^2 + (camd.y - y).^2) );
                if r < Rsearch
                    xcam(i) = camd.x(ir);
                    ycam(i) = camd.y(ir);
                    found = 1;
                else
                    beep();
                    fprintf('Bad click, try again\n');
                end
            end
            hold on; plot(xcam(i), ycam(i), 'ys', 'markersize', 15);
        end

        %% The relation between the clicked points is used to find the rest
        cx = fit([xpat',ypat'], xcam', 'poly11');
        cy = fit([xpat',ypat'], ycam', 'poly11');
        % estimated position of all the dots in the camera image
        x = cx(patd.x, patd.y);   
        y = cy(patd.x, patd.y);
        % Plot estimated positions to the previous image
        hold on; plot(x,y,'gx'); hold off
        
        % Relate the camera dots to the estimated position if they are 
        % within a given -small- radius.
        Rsearch = 50;
        % clear k
        k = nan*ones(numel(x),1);
        for i = 1:numel(x)
            [r,ir] = min( sqrt((camd.x - x(i)).^2 + (camd.y - y(i)).^2) );
            if r < Rsearch
                k(i) = ir;
            end
        end
        k1 = find(~isnan(k));   % patdots indices of valid points
        k2 = k(~isnan(k));      % camdots indices of valid points
        
        hold on; 
        plot(camd.x(k2), camd.y(k2), 'rs', 'markersize',14); 
        hold off;

        cx = fit([patd.x(k1), patd.y(k1)], camd.x(k2), 'poly33');
        cy = fit([patd.x(k1), patd.y(k1)], camd.y(k2), 'poly33');

        x = cx(patd.x, patd.y);
        y = cy(patd.x, patd.y);
        hold on; plot(x,y,'ys','markersize',10); hold off
        title('yellow squares are the predicted positions'); 
        
        % Once satisfied with the result exit. Either y and enter or just
        % Enter works.
        a = input('Are you satisfied with the result? <y> : ', 's');
        satisfied = isempty(a) || a(1) == 'y';
    end
    
    %% Sort all the patterns using the fit found
    clf();
    for ii = 1:numel(Ipat)
        imagesc(Icam{ii}); colormap(([1:255; 1:255; 1:255])'/256);
        title(sprintf('ii = %d', ii));

        patd = patdots(ii);
        camd = camdots(ii);
        
        % estimated position of all the dots in the camera image
        x = cx(patd.x, patd.y);
        y = cy(patd.x, patd.y);
        hold on; plot(x,y,'gx'); hold off

        Rsearch = 50;
        camdots(ii).x = nan * ones(numel(x),1);
        camdots(ii).y = nan * ones(numel(x),1);
        for i = 1:numel(x)
            [r,ir] = min( sqrt((camd.x - x(i)).^2 + (camd.y - y(i)).^2) );
            if r < Rsearch
                camdots(ii).x(i) = camd.x(ir);
                camdots(ii).y(i) = camd.y(ir);
            end
        end
        hold on; 
        plot(camdots(ii).x, camdots(ii).y, 'rs', 'markersize', 14); 
        hold off;
        pause(0.1);
    end
end
