%% #######################################################################%
%------------- MATLAB Script used to calibrate the projector -------------%
%#########################################################################%
% This is an example script used for the projector calibration. For a
% complete description, consult [Gonzalez-Vera, A. S., Wilting, T. J. S., 
% Holten, A. P. C., van Heijst, G. J. F. & Duran-Matute,M. "High-resolution
% single-camera photogrammetry: incorporation of refraction at a fluid 
% interface". Experiments in Fluids 61(1),3]. This is a simple template 
% that can be used by specifying directories, values and files necessary.
% See directory (Examples\Calibration) for a projector calibration 
% example. Files, image directories, and values necessary are included in 
% for each example.
%
%% The variables given here are given as 0 or 1 (logical) and are used to:
Lens_pass = 0; % Force the use of an additional plane in the calibration.
               % This is not usually recommended.

%% Data and function root directories
addpath('Common_Functions\'); % Folder containing the functions used
                                    % in the script.
SavLoc = ''; % Root location where data is saved. 
ImageDir = ''; % Location of the photographs used for the calibration.
zlevels = []; % these are the heights in mm (from the bottom) at which 
              % the patterns were photographed.
% Name of the file that contains the dot pattern images
PatternName = ''; 
npos = 120; % Total number of patterns in the file. 
            % The value of npos should match the number of images  
            % photographed and be ordered accordingly for a correct 
            % calibration. However, a lower value can be selected if 
            % pattern and images are consistent.
        
% The following values are used to finds circles with radii with 
% the range specified. Dependent on dot size.
Im_Rrange = []; % For dot detection in the images sent to the projector. 
                % Range values depend on the size of the dots in the
                % images. size [2x1]. recomended above 4px for minimum
                % value.
Pr_Rrange = []; % For the dot detection of the photographed patterns
% Name and location where the Projector calibration file will be saved and
% name of the camera calibration file (.mat) to load.
Save_name = '';
Cam_File = '';
         
% Load the polynomials from the camera calibration
cm = load(Cam_File, 'cx0', 'cy0', 'cdx', 'cdy', 'Plens');

% Number of levels at which patterns where projected.
nz = numel(zlevels);

%Polynomial fitting type.
fittype = 'poly33';

% Approximate position (in mm) of the projector lens.
Plens = []; % 3x1 : (x,y,z)

%% Determine the positions of the dots in each of the patterns sent to the 
%  projector.

% Data file where the (pixel) position of the dot patterns sent to the
% projector are saved.
infofile = ''; 
        
create_new = true;
if exist(infofile, 'file')% Find if file with dot positions already exists.
    a = input('using the saved pattern dots positions <y> ? ','s');
    if isempty(a) || a(1) == 'y'
        % If file already exists (and chosen), load the data
        load(infofile, 'Pattern_images_file', 'patdots');
        % If no file is found (or not chosen), a new file will be created 
        create_new = false;
   end
end

%
if create_new
    % Load the file with the dot pattern images.
    Pattern_images_file = '';
    PatternPrj = load(Pattern_images_file);
    
    % Determine the amount of pattern images
    npos = size(PatternPrj.DotsToProject,3);

    for ii = 1:npos
        Ipat = PatternPrj.DotsToProject(:,:,ii);
        % Determine the location of each of the dots in the images sent to
        % the projector.
        patdots(ii) = find_circles(Ipat, Pr_Rrange); 
    end
    % Save data 
    save(infofile, 'Pattern_images_file', 'patdots');
end

%% Determine the positions of dots of the photographed patterns

% Name the data file where the (pixel) position of the photographed dot 
% patterns will be saved.
infofile = ''; 

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
        for ii = 1:npos % for all photographed patterns
            %figure
            imgname = [ImageDir imgfold, list(ii,:)];
            im = imread(imgname);
            im = im2bw(im,.3); % Threshold the brightness of the image. 
                               % helps in the detection.
            image_names{iz,ii} = imgname;
            camdots(iz,ii) = find_circles(im, Im_Rrange); % locate dots and
                                                          % save position
        end
    end
    % Save the data.
    save(infofile, 'image_names', 'camdots', 'zlevels');
end

%% Sort the dots found in the photographs by relating the patterns in the 
% photogrpahs with the patterns sent to the projector. Dots not beloging 
% are eliminated.

% Name the data file where the position of the sorted dots will be saved
infofile = ''; 

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
        PatternPrj = load(Pattern_images_file);
    end
    % Compare each dot pattern projected with its photograph for all levels
    for iz = 1:nz
        fprintf('processing iz = %d\n', iz);
        for ii = 1:npos
            % Pattern sent to projector.
            Ipat{ii} = PatternPrj.DotsToProject(:,:,ii); 
            % Pattern captured by the camera.
            Icam{ii} = imread(image_names{iz,ii});
        end
        % The comparisson is done by the function 
        % [sort_projector_calibration_data]. It requires user input to 
        % relate the projected dots to the ones captured by the camera. 
        % Graphical interface is required for the input.
        image_dots(iz,:) = sort_projector_calibration_data(Ipat, Icam, ...
                                                patdots, camdots(iz,:));
    end
    % Save file
    save(infofile, 'image_dots');
end

%% Combine all the dots from each level and convert the projected and 
% photographed dots coordinates to 3D

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
end

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