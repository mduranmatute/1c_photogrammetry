%% =======================================================================%
% This is an example script used for the digital Projector calibration. The 
% case presented is related to the camera calibration files:
% Camera_LinesFit_Example. Each step taken in this script is briefly 
% described. For a complete description consult [Gonzalez-Vera,A. S.,
% Wilting,T. J. S., Holten,A. P. C., van Heijst,G. J. F. & Duran-Matute,M. 
% "High-resolution single-camera photogrammetry: incorporation of 
% refraction at a fluid interface". Experiments in Fluids 61(1),3]. 
% Processed files are located in the directory: (Photogrammetry\Examples\
% Data_Files\MATLAB\Data_process_steps). The resulting calibration file
% is saved in (Photogrammetry\Examples\Data_Files\MATLAB\Calibration_data)
% The resulting calibration is used in the "high" resolution measurement
% of the ondulated plate example case presented in the directory:
% (Photogrammetry\Examples\Measurements\).  
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
ImageDir    = 'Photographs\Projector\'; % Location of the photographs used
                                                    %  for the calibration.
% zlevels are the heights in mm (from the bottom) at which the patterns 
% were photographed.                                                    
zlevels = [2.150,15.450,60.475,75.575,105.525,135.575,150.350,195.350];
% Name of the file that contains the dot pattern images        
PatternName = 'CalibrationPattern_Example.mat';        
npos = 120; % Total number of patterns in the file. 
            % The value of npos should match the number of images  
            % photographed and be ordered accordingly for a correct
            % calibration.
% The following values are used to finds circles with radii with the range
% specified. Dependent on the projected and photographed dot size.
Im_Rrange = [7,16]; % For dot detection in the images sent to the projector 
Pr_Rrange = [5,9];  % For the dot detection of the photographed patterns

% Name and location where the Projector calibration file will be saved and
% name of the camera calibration file (.mat) to load.
Save_name = [SavLoc '\Calibration_data\Projector_LinesFit_Example.mat'];
Cam_File = [SavLoc '\Calibration_data\Camera_LinesFit_Example.mat'];
         
% Load the polynomials from the camera calibration
cm = load(Cam_File, 'cx0', 'cy0', 'cdx', 'cdy', 'Plens');

% Number of levels at which patterns where projected.
nz = numel(zlevels);

%Polynomial fitting type.
fittype = 'poly33';

% Approximate position (in mm) of the projector lens.
Plens = [4, 560, 1100]; % Unlike the camera calibration, the position of 
                        % the projector is given by the user.

%% Determine the positions of the dots in each of the patterns sent to the 
%  projector.

% Data file where the (pixel) position of the dot patterns sent to the
% projector are saved.
infofile = [SavLoc '\Data_process_steps\Calibration\'...
                        'DotPosInPattern_ProjCal.mat']; 
        
create_new = true;
if exist(infofile, 'file')% Find if file with dot positions already exists.
    a = input('using the saved pattern dots positions <y> ? ','s');
    if isempty(a) || a(1) == 'y'
        % If file already exists (and chosen), load the data
        load(infofile, 'PatternName', 'patdots');
        % If no file is found (or not chosen), a new file will be created 
        create_new = false;
   end
end

%
if create_new
    % Load the file with the dot pattern images.
    Pattern_images_file = [SavLoc '\Patterns\' PatternName];
    PatternPrj = load(Pattern_images_file);
    
    % Determine the amount of pattern images
    npos = size(PatternPrj.DotsToProject,3);

    for ii = 1:npos
        Ipat = PatternPrj.DotsToProject(:,:,ii);
        % Determine the location of each of the dots in the images sent to
        % the projector.
        patdots(ii) = find_circles(Ipat, Pr_Rrange); 
        
        if(PlotCheck) % visualize location process
            image(Ipat * 255);
            % colormap for 8 bit grayvalue images
            colormap(([1:255; 1:255; 1:255])'/256);     
            title(sprintf('ii = %d\n', ii));
            hold on
            plot(patdots(ii).x, patdots(ii).y, 'r+');
            hold off
            pause(0.1);
        end
    end
    % Save data 
    save(infofile, 'PatternName', 'patdots');
end

%% Determine the positions of dots of the photographed patterns

% Name the data file where the (pixel) position of the photographed dot 
% patterns will be saved.
infofile = [SavLoc '\Data_process_steps\Calibration\'...
                    'DotPosInPhoto_ProjCal.mat']; 

% Load values or create file with the position of the dots
create_new = true;
if exist(infofile, 'file')
    a = input('using the saved image dots positions <y> ? ','s');
    if isempty(a) || a(1) == 'y'
        load(infofile, 'image_names', 'camdots');
        create_new = false;
   end
end

if create_new
    for iz = 1:nz % Do for each level 
        % Add all the folders in the directory with z_1, z_2, ... z_n
        % in a list.
        imgfold = sprintf('z_%d/', iz);
        list = ls([ImageDir imgfold, '*.png']);
        fprintf('processing iz = %d\n', iz);
        for ii = 1:npos % for all photographed patterns
            imgloc = [ImageDir imgfold, list(ii,:)];
            im = imread(imgloc);
            imgname = [imgfold list(ii,:)];
            if(PlotCheck) % Show pattern
                image(im); 
                % colormap for 8 bit grayvalue images
                colormap(([1:255; 1:255; 1:255])'/256);     
                title(sprintf('iz = %d, ii = %d\n', iz, ii));
            end

            im = im2bw(im,.3); % Threshold the brightness of the image. 
                               % helps in the detection.

            image_names{iz,ii} = imgname;
            camdots(iz,ii) = find_circles(im, Im_Rrange); % locate dots and
                                                          % save position
            if(PlotCheck)
                hold on
                plot(camdots(iz,ii).x, camdots(iz,ii).y, 'r+');
                hold off
                pause(0.1);
            end 
        end
    end
    % Save the data.
    save(infofile, 'image_names', 'camdots', 'zlevels');
end

%% Sort the dots found in the photographs by relating the patterns in the 
% photogrpahs with the patterns sent to the projector. Dots not beloging 
% are eliminated.

% Name the data file where the position of the sorted dots will be saved
infofile = [SavLoc '\Data_process_steps\Calibration\' ...
                    'SortedDotPos_ProjCal.mat']; 

% Load or create new file.        
create_new = true;
if exist(infofile, 'file')
    a = input('using the saved clicked positions <y> ? ','s');
    if isempty(a) || a(1) == 'y'
        load(infofile, 'image_dots');
        create_new = false;
   end
end
if create_new 
    % If new file is created, clear previous variables containing dot 
    % positions and detections and read from previously saved files.
    clear Ipat Icam
    if ~exist('DotsToProject', 'var')
        % if pattern file does not exist load.
        PatternPrj = load([SavLoc '\Patterns\' PatternName]);
    end
    % Compare each dot pattern projected with its photograph for all levels
    for iz = 1:nz
        fprintf('processing iz = %d\n', iz);
        for ii = 1:npos
            % Pattern sent to projector.
            Ipat{ii} = PatternPrj.DotsToProject(:,:,ii); 
            % Pattern captured by the camera.
            Icam{ii} = imread([ImageDir image_names{iz,ii}]);
        end
        % The comparisson is done by the function 
        % [sort_projector_calibration_data]. It requires user input to 
        % relate the projected dots to the ones captured by the camera. 
        % Graphical interface is required for the input.
        image_dots(iz,:) = sort_projected_dots(Ipat, Icam, ...
                                                patdots, camdots(iz,:));
    end
    % Save file
    save(infofile, 'image_dots');
end

%% Combine all the dots from each level and convert the projected and 
% photographed dots coordinates to 3D

if(PlotCheck) %clear figures if present
    clf()
end

% Initialize a variable where the positions of the dots sent to the
% projector will be combined
clear pattern
pattern.x = [];
pattern.y = [];

% Fill varaible with viable dot positions
for ii = 1:npos
    pattern.x = [pattern.x; patdots(ii).x];
    pattern.y = [pattern.y; patdots(ii).y];
end

% Now a variable where the positions of the dots found in the photographs
% will be combined is initialized.
for iz = 1:nz
    cam.x = [];
    cam.y = [];
    for ii = 1:npos
        cam.x = [cam.x; image_dots(iz,ii).x]; 
        cam.y = [cam.y; image_dots(iz,ii).y]; 
    end
    
    % Here, the pixel coordinates of the dots in the photographs that were 
    % found to belong to the projected dots are converted to 3D-planes
    x = cam.x;
    y = cam.y;
    observed(iz).x = cm.cx0(x,y) + zlevels(iz) * cm.cdx(x,y);
    observed(iz).y = cm.cy0(x,y) + zlevels(iz) * cm.cdy(x,y);
    observed(iz).z = zlevels(iz) * ones(size(x));
    
    if(PlotCheck) % plot the viewed dots
        plot(observed(iz).x, observed(iz).y, 'b*');
    end
    
    % Since the calibration plate did not cover the complete projection 
    % area, the dots that fall over the plate are removed. The position 
    % values at which the dots are eliminated are given by the plates
    % dimensions. 
    k = find(observed(iz).x < 25 | observed(iz).x > 575 | ...
             observed(iz).y < 25 | observed(iz).y > 575);
    observed(iz).x(k) = nan;
    observed(iz).y(k) = nan;
    observed(iz).z(k) = nan;
    
    if(PlotCheck)
        plot3(observed(iz).x, observed(iz).y, observed(iz).z, '.');
        hold on
    end
end
%hold off

%% Find formulas to convert projector pixel coordinates to lines in space:

% 1. Remove pattern dots which were not detected on all images
np = numel(observed(iz).x);
k = zeros(nz,np);
for iz = 1:nz
    xp = observed(iz).x;
    k(iz,:) = ~isnan(xp);   % i.e. 1 = detected, 0 = not detected
end
kk = find(sum(k) == nz)';   % indices of the points visible in all planes

% 2. Fit the checked data set 
nZ = nz;
if(Lens_pass)
    % lines have to go throught the lens centre, so an extra zlevel
    % is added to force this
    nZ = nz+1;
end

% Pixel coordinates of the dot patterns sent to the projector
xp = pattern.x(kk);          
yp = pattern.y(kk);

for iz = 1:nZ
    % create the data set to fit
    if iz <= nz
        % 3d-coordinates in the (iz) plane 
        X = observed(iz).x(kk);            
        Y = observed(iz).y(kk);
        Z = observed(iz).z(kk);
    else
        % This is the added plane given by the projector lens 
        X = Plens(1) * ones(size(kk));
        Y = Plens(2) * ones(size(kk));
        Z = Plens(3) * ones(size(kk));
    end
    if iz == 1
        Sx  = X;    Sy  = Y;    Sz = Z;
        Sxz = X.*Z; Syz = Y.*Z; Szz = Z.^2;
    else
        Sx  = Sx+X;     Sy  = Sy+Y;     Sz  = Sz+Z;
        Sxz = Sxz+X.*Z; Syz = Syz+Y.*Z; Szz = Szz+Z.^2;
    end
end
xrc = (nZ*Sxz - Sx.*Sz)./(nZ.*Szz - Sz.^2);
xof = (Sx - xrc.*Sz)./nZ;		

yrc = (nZ*Syz - Sy.*Sz)./(nZ*Szz - Sz.^2);
yof = (Sy - yrc.*Sz)./nZ;		

cx0 = fit([xp,yp], xof, fittype);
cdx = fit([xp,yp], xrc, fittype);
cy0 = fit([xp,yp], yof, fittype);
cdy = fit([xp,yp], yrc, fittype);

%% Save projector calibration fit to file. 
% Save vector origin a0 and its components a, the projector Lens Location 
% and the heights at which the projections photographed were located.
save(Save_name, 'cx0', 'cy0', 'cdx', 'cdy', 'zlevels', 'Plens');

%% Show calibration: 
if(PlotCheck)
%  1- Construct 3D lines from projected dots pixel coordinates
%  2- Plot the 3D pattern dots detected by the camera 

% convert the dot-positions to lines in space
    xp = patdots(1).x;
    yp = patdots(1).y;

    x0 = cx0(xp,yp);
    dx = cdx(xp,yp);
    y0 = cy0(xp,yp);
    dy = cdy(xp,yp);

    h1 = 0;             % draw the lines between the planes h1 and h2
    h2 = 1500;

    X = [x0,x0] + dx*[h1,h2]; %X = x0 + dx*[h1,h2]; %Matlab 2017 and above
    Y = [y0,y0] + dy*[h1,h2]; %Y = y0 + dy*[h1,h2];
    Z = ones(size(xp)) * [h1,h2];
    plot3(X',Y',Z','-');
    hold on

% add the observed dots in the planes to the plot
    for iz = 1:nz
        x = image_dots(iz,1).x;
        y = image_dots(iz,1).y; 
        cm_x0 = cm.cx0(x,y);
        cm_y0 = cm.cy0(x,y);
        cm_dx = cm.cdx(x,y);
        cm_dy = cm.cdy(x,y);

        X = cm_x0 + cm_dx*zlevels(iz);
        Y = cm_y0 + cm_dy*zlevels(iz);
        Z = zlevels(iz) * ones(size(X));

        plot3(X, Y, Z, 'kx');
    end
end