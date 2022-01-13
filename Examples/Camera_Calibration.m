%% =======================================================================%
% This is an example script used to chow the steps taken to calibrate the
% camera for measurements. Each step taken in this script is briefly 
% described. For a complete description consult [Gonzalez-Vera, A. S.,
% Wilting, T. J. S., Holten, A. P. C., van Heijst, G. J. F. & Duran-Matute,
% M. "High-resolution single-camera photogrammetry: incorporation of 
% refraction at a fluid interface". Experiments in Fluids 61(1),3].
% Processed files are located in the directory:
% (Photogrammetry\Examples\Data_Files\MATLAB\Data_process_steps\).
% The resulting calibration file is saved in:
% (Photogrammetry\Examples\Data_Files\MATLAB\Calibration_data).  
% The resulting file is used in the "High" resolution case presented in the 
% (Photogrammetry\Examples\Measurements\) directory.
%=========================================================================%

%% The variables given here are given as 0 or 1 (logical) and are used to:
Lens_pass = 0; % Force the use of an additional plane in the calibration.
               % This is not usually recommended.
PlotCheck = 0; % Use to plot some of the results. Used for diagnostics
               %  and to visualize some of the results. 

%% Give function root directories and give (or load) needed data 
addpath('..\..\Common_Functions\'); % Folder containing the functions used
                                    % in the script.
SavLoc = '..\Data_Files\MATLAB'; % Root location where data is saved. 
               
imagedir = 'Photographs\Camera\'; % Location of the photographs used for
                                  % the calibration.
zlevels = [14.35, 59.375, 74.475, 104.425, 134.475, 149.25, 194.25];                                   
% These are the heights from the bottom (in mm) at which the calibration 
% plate was photographed.

%Name and location where the camera calibration file will be saved.
Save_name = [SavLoc '\Calibration_data\Camera_LinesFit_Example.mat'];

nz = numel(zlevels); % Number of levels at which the calibration plate was
                     % photographed.

% Approximate horizontal position (in pixels) of the camera lens.
Plens_pixels = [811.7, 628.9];

% Approximate vertical position (in mm) of the camera lens.
Zlens = [2060, 2060]; % mm

%% This part of the script is used to find dots in an image
infofile = [SavLoc '\Data_process_steps\Calibration\' ...
              'PlateDotsRaw_CamCal.mat']; % data file where the unprocessed
                                          % dot location is stored.
create_new = true;
if exist(infofile, 'file') % Find if file already exists.
    a = input('using the saved dot positions <y> ? ','s');
    if isempty(a) || a(1) == 'y' % If file already exists and chosen, load
                                 % data.
        load(infofile, 'I', 'dotpos');
        create_new = false; 
   end
end

if create_new %If file does not exist find dots and create file.
    list = ls([imagedir 'z*.png']); % make a list with the images.
    for iz = 1:nz % cycle though each image in the list
        imgname = [imagedir, list(iz,:)];
        im = imread(imgname);
        I{iz} = im;
        dotpos(iz) = find_circles(im, [5,12]); %Use function to identify
        % and locate the position (in pixels) of dots with a size range of
        % [5,12] pixels. This size range is chosen for the calibration 
        % example cases shown here. This range was optimized by trial and
        % error for these particular calibration images.
        
        if(PlotCheck) % Plot photographs and dot positions found 
            image(im); 
            % colormap for 8 bit grayvalue images
            colormap(([1:255; 1:255; 1:255])'/256);     
            title(sprintf('iz = %d\n', iz));
            hold on; plot(dotpos(iz).x, dotpos(iz).y, 'r+'); hold off
            pause(0.2);
        end
    end
    save(infofile, 'I', 'dotpos', 'zlevels'); % Save data to file
end


%% grid coordinates in world coordinates of the calibration plate
nx = 23;    % number of dots in x and y
ny = 23;
dx = 25;    % mesh size in mm
dy = 25;
xrange = 300 + (-11:11)*dx;       % take the centre of the axes 0,0
yrange = 300 + (-11:11)*dy;       % y-axis to down
[X,Y] = meshgrid(xrange, yrange);

clear grid
for iz = 1:nz
    grid(iz).X = X(:);
    grid(iz).Y = Y(:); % reshape the matrices to 1 dim. arrays
end

%% Sorting the dots found in the images and eliminating the dots not 
% belonging to points on the calibration plate. Images are necessary since
% human input is required to relate the dots in the photographs to the ones
% found by the [find_circles] function.

infofile = [SavLoc '\Data_process_steps\Calibration\' ...
            'Sorted_CamDots_CamCal.mat'];

create_new = true; % As before, determine if file already exists and load
                   % if chosen.
if exist(infofile, 'file')
    a = input('using the saved clicked positions <y> ? ','s');
    if isempty(a) || a(1) == 'y'
        load(infofile, 'grid_dots', 'img_dots');
        create_new = false;
   end
end

% If file does not exist or chosen, sort which dots found in the images 
% correspond to actual dots in the calibration grid. 
if create_new
    for iz = 1:nz
        [grid_dots(iz), img_dots(iz)] = sort_camera_calibration_data(...
            I{iz}, grid(iz), dotpos(iz));
    end
    save(infofile, 'grid_dots', 'img_dots'); % Save data to file.
end

%% fitting Z-planes 
%  loop 1 : fit without plate shift corrections
%  loop 2 : fit with plate shift corrections
fittype = 'poly33';
for loop = 1:2
    for iz = 1:nz
        % only fit existing coordinates
        X  = grid_dots(iz).x; Y  = grid_dots(iz).y; 
        xp = img_dots(iz).x;  yp = img_dots(iz).y; 
       
        % fitting grid coordinates -> pixel coordinates
        fitp(iz).cx = fit([X,Y], xp, fittype);
        fitp(iz).cy = fit([X,Y], yp, fittype);
        
        xfit = fitp(iz).cx(X,Y);
        yfit = fitp(iz).cy(X,Y);
        
        if(PlotCheck)
        % showing the fit result
        image(I{iz});
        colormap(([1:255; 1:255; 1:255])'/256);
        title(sprintf('pixel coordinate grid points, z = %f',zlevels(iz)));
        hold on; plot(xfit, yfit, 'b+'); hold off
        end

        % fitting pixel coordinates -> real-world coordinates
        fitp(iz).cX = fit([xp,yp], X, fittype);
        fitp(iz).cY = fit([xp,yp], Y, fittype);
        
        xfit = fitp(iz).cX(xp,yp); 
        yfit = fitp(iz).cY(xp,yp);
        
        if(PlotCheck)
        plot(X,Y,'o');
        hold on; plot(xfit, yfit, 'b+'); hold off
        end
    end

    if loop == 1
        %% Calculations of the plate shifts
        % A point straight under the lens must have the same grid  
        % coordinates independent from the height. 
        % The position of the plate at the bottom is taken as reference.
        clear X Y Xshift Yshift
        for iz = 1:nz
            X(iz) =  fitp(iz).cX(Plens_pixels(1), Plens_pixels(2)); 
            Y(iz) =  fitp(iz).cY(Plens_pixels(1), Plens_pixels(2)); 
        end
        Xshift = X - X(1);
        Yshift = Y - Y(1);
        Plens = [X(1), Y(1)];

        % correct the grid positions for the shift and redo the calibration
        for iz = 1:nz
            grid_dots(iz).x = grid_dots(iz).x - Xshift(iz);
            grid_dots(iz).y = grid_dots(iz).y - Yshift(iz);
        end
     end
end

%% Z location of the lens from the magnification
% assumption: the inversed magnification goes linear with the distance to
% the pinhole lens and will be 0 at the pinhole. Uses function 
% [magnification_A].

for iz = 1:nz
    [Mx(iz),My(iz)] = magnification_A(fitp(iz));
end
meanMinv = 1./mean([Mx;My]);

if(PlotCheck)
    plot(zlevels, 1./Mx, '*-', zlevels, 1./My, '*-');
end

c = fit(meanMinv', zlevels', 'poly1');
fprintf('by magn. : Z-lens = %.3f\n', c(0));
Plens(3) = c(0);

%% --- Finding formulas to convert pixel coordinates to lines in space ---
% For a more in depth description of the conversion via third order
% polynomial see (Soloff, S. M., Adrian, R. J., and Liu, Z.-C. (1997). 
% "Distortion compensation for generalized stereoscopic particle image 
% velocimetry". Measurement Science and Technology, 8(12):1441.

%------ create an artificial grid within the pixel area populated with dots
x = img_dots(nz).x;     % highest plane
y = img_dots(nz).y;
x1 = min(x); x2 = max(x); dx = (x2-x1)/22;
y1 = min(y); y2 = max(y); dy = (y2-y1)/22;
[x,y] = meshgrid(x1:dx:x2, y1:dy:y2);
art_xp = x(:);
art_yp = y(:);

xp = art_xp;
yp = art_yp;
Nz = nz;

% The calibration can be finetuned by using the pinhole of the lens as an 
% extra plane. In other words, force the calibration lines to go through
% the lens centre, by adding an extra zlevel. This is not recomended unless
% the position of the camera is well determined. For the cases treated
% here, this option was not chosen.
if(Lens_pass)
    Nz = nz+1;
end

ic = 1;
for iz = 1:Nz
    % calculate the real-world coordinates for xp,yp in each iz plane.
    if iz <= nz
        X = fitp(iz).cX(xp,yp);
        Y = fitp(iz).cY(xp,yp);
        Z = zlevels(iz) * ones(size(X));
    else
        X = Plens(ic,1) * ones(size(X));
        Y = Plens(ic,2) * ones(size(X));
        Z = Plens(ic,3) * ones(size(X));
    end
    if iz == 1
        Sx  = X;    Sy  = Y;    Sz = Z;
        Sxz = X.*Z; Syz = Y.*Z; Szz = Z.^2;
    else
        Sx  = Sx+X;     Sy  = Sy+Y;     Sz  = Sz+Z;
        Sxz = Sxz+X.*Z; Syz = Syz+Y.*Z; Szz = Szz+Z.^2;
    end
    
    if(PlotCheck)
        U{iz} = X;      % only used for plotting
        V{iz} = Y;
        W{iz} = Z;
    end 
end
xrc = (Nz*Sxz - Sx.*Sz)./(Nz.*Szz - Sz.^2);
xof = (Sx - xrc.*Sz)./Nz;		

yrc = (Nz*Syz - Sy.*Sz)./(Nz*Szz - Sz.^2);
yof = (Sy - yrc.*Sz)./Nz;		

cx0 = fit([xp,yp], xof, fittype);
cdx = fit([xp,yp], xrc, fittype);
cy0 = fit([xp,yp], yof, fittype);
cdy = fit([xp,yp], yrc, fittype);

% % x0,y0 is the intersection point with the z=0 plane
% % dx,dy is the pointing vector, with dz = 1

%% Save calibration fit to file. 
% Save vector origin A0 and its components a, the camera Lens Location and
% the heights at which the calibration plate was located.
save(Save_name, 'cx0', 'cy0', 'cdx', 'cdy', 'Plens', 'zlevels');

%% Show results: 
% construct lines in 3D from pixel coordinates and plot also the
% corresponding 3D points. A good calibration has all lines passing through
% (or very close to) a small region that closely matches the real-world
% position of the camera.

if(PlotCheck)
    xp = art_xp;
    yp = art_yp;

    x0 = cx0(xp,yp);
    dx = cdx(xp,yp);
    y0 = cy0(xp,yp);
    dy = cdy(xp,yp);

    figure()
    h1 = 0;
    h2 = 2500;
    plot3([x0 + dx*h1, x0 + dx*h2]', [y0 + dy*h1, y0 + dy*h2]', ...
            ones(size(x0))*[h1,h2], '-');
    hold on
    % grid planes
    for iz = 1:nz
        plot3(U{iz},V{iz},W{iz},'k.');
        hold on
    end
    hold off
end
