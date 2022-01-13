%% #######################################################################%
%------------ MATLAB Script used to measure submerged surfaces -----------%
%#########################################################################%
% This is an example template to measure different surfaces, for a  
% description see [ReadMe.txt]. It requires the input of the the 
% calibration files given by the Projector and Camera calibration scripts, 
% the projected patterns file, and the photographs of the patterns on the
% surface to be measured. Each step taken here is briefly described. For a
% complete description consult [Gonzalez-Vera,A. S., Wilting,T. J. S.,
% Holten,A. P. C., van Heijst,G. J. F. & Duran-Matute,M.  "High-resolution
% single-camera photogrammetry: incorporation of  refraction at a fluid
% interface". Experiments in Fluids 61(1),3]. Examples of camera and
% projector calibrations and a few test measurements can be found it the
% (Examples\) directory. The template includes the added calculations 
% necessary to for the presence of a flat or parabolic water interface.

addpath('Common_Functions\'); % Folder containing the functions used
                                    % in the script.
zwater = 0; % Level, or depth, of the water layer (in metres).
theta = 0.0; % A small roation of the pattern is implemented.
Omega = 0.0; % Rotation rate (rad/s) of the set-up. Examples shown here did
             % not rotate except for a single case (measurement = 3).
Normal = [0,0,1]; % Define the Normal vector (points upward)
load('refraction', 'Nwater', 'Nair'); % Load air and water index of 
                                      % refraction values

% Directory where pattern photographs are located
imagedir    = '';
% Name of files where part of the image processing steps are saved
camdotsfile = '';
sortedfile  = '';
% Name of file where the result (reconstruction) is saved
Reconstruct = '';
% Pattern used for measurements
Pattern_Fl  = '';
        
thrsh = 0; % Thereshold value that will be ignored when processing 
                   % the photographs of the projected dots.
Rrange = []; % Range of radii [2x1] used for the dot detection in the 
                         % images sent to the projector
DotR = []; % Range of radii [2x1] used for the dot detection in the
                       % photographed patterns
Interp_Div = 200; %Used for grid interpolation.

%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
%=================== Determine positions of the dots sent to the projector%
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

%-------------------------- Data file where the (pixel) position of the dot 
%-------------------------- patterns sent to the projector are saved.
PatternFile  = Pattern_Fl;
infofile = ''; % Rename

% Find if file with dot positions already exists. If file already exists 
% (and chosen), load the data. If no file is found (or not chosen), a  
% variable marks that a new file will be created. 
create_new = true;
if exist(infofile, 'file')
    a = input('Use the saved dot pattern positions <y> ? ','s');
    if isempty(a) || a(1) == 'y'
        load(infofile, 'patdots');
        npos = numel(patdots);
        create_new = false;
   end
end

% Create a file containing the positions (x,y) of the dots to be projected. 
if create_new
    % Load the file with the dot patterns (those sent to the projector).
    PatternPrj = load(PatternFile);
    DotPatterns = PatternPrj.DotsToProject;
    
    % Determine the amount of pattern images
    npos = size(DotPatterns,3);
    for ii = 1:npos
        Ipat = DotPatterns(:,:,ii);
        image(Ipat * 255);
        % colormap for 8 bit grayvalue images
        colormap(([1:255; 1:255; 1:255])'/256);     
        title(sprintf('ii = %d\n', ii));
        
        % Determine the location of each of the dots in the images sent to
        % the projector.
        patdots(ii) = find_circles(Ipat, Rrange);
        hold on; plot(patdots(ii).x, patdots(ii).y, 'r+'); hold off
        pause(0.5);
    end
    save(infofile, 'patdots');
end

%----------------- combine all the projected patterns into single data set.
pattern.x = [];
pattern.y = [];
for ii = 1:npos
    pattern.x = [pattern.x; patdots(ii).x];
    pattern.y = [pattern.y; patdots(ii).y];
end

%-------------------------------- Due to the movement of the set-up a small 
%----------------------- rotation of the projected dot patterns is applied.
cenx = mean(pattern.x); % Find mean values along x and y directions
ceny = mean(pattern.y);

pattern.x = cenx + cos(theta)*(pattern.x - cenx) ...
                 - sin(theta)*(pattern.y - ceny);
pattern.y = ceny + sin(theta)*(pattern.x - cenx) ...
                 + cos(theta)*(pattern.y - ceny);

%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
%========================= Determine the position of the dots photographed%
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

% Name the data file where the (pixel) position of the photographed dot 
% patterns will be saved.
infofile = camdotsfile;

% Load values or indicate that a file with the position of the dots
% requires to be created.
create_new = true;
if exist(infofile, 'file')
    a = input('using the saved image dots positions <y> ? ','s');
    if isempty(a) || a(1) == 'y'
        load(infofile, 'camdots', 'image_names');
        create_new = false;
   end
end

% Create file containing the positions (x,y) of the dots photographed. 
if create_new    
    % Make a list of all saved photographed patterns
    list = ls([imagedir, '*.png']);
    clear camdots  
    for ii = 1:npos % cycle through each pattern projected
        imgname = [imagedir, list(ii,:)]; 
        im = imread(imgname);
        %--------------------------- Display photograph, find dots and plot
        imagesc(im);
        colormap(([1:255; 1:255; 1:255])'/256); 
        title(sprintf('Dot pattern at start. Image ii = %d\n', ii));
        im(im<thrsh) = 0; % Remove image noise to help in the dot detection
        camdots(ii) = find_circles(im,DotR); %#ok<SAGROW>
        hold on; plot(camdots(ii).x, camdots(ii).y, 'r+',...
            'MarkerSize',7,'LineWidth',1.5); hold off
        %------------------------------------------------------------------
        image_names{ii} = imgname;
        pause(0.0);
    end
    save(infofile, 'camdots', 'image_names');
end

%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
%======================= Sort dots found in the photographs and eliminiate%
%==========================  measurements not belonging to the dot pattern%
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

% Name the data file where the position of the sorted dots will be saved
infofile = sortedfile;
create_new = true;

% Load values or indicate that a file with the sorted dot posiotion needs
% to be created.
if exist(infofile, 'file')
    a = input('Use the saved clicked positions <y> ? ','s');
    if isempty(a) || a(1) == 'y'
        load(infofile, 'image_dots');
        create_new = false;
   end
end

if create_new
    % If new file is created, clear previous variables containing dot 
    % positions and detections and read from previously saved files.
    PatternPrj = load(PatternFile);
    clear Ipat Icam image_dots
    for ii = 1:npos
        % Pattern sent to projector.
        Ipat{ii} = PatternPrj.DotsToProject(:,:,ii); %#ok<SAGROW>
        % Pattern captured by the camera.
        Icam{ii} = imread(image_names{ii}); %#ok<SAGROW>
    end
    % The comparisson is done by the function 
    % [sort_projector_calibration_data]. It requires user input to 
    % relate the projected dots to the ones captured by the camera. 
    % Graphical interface is required for the input.
    image_dots = sort_projector_calibration_data(Ipat, Icam, ...
                    patdots, camdots);
    
    save(infofile, 'image_dots');
end

%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
%================================ Determine the shape of the water surface
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

h0 = zwater;     % water layer height in metres (initial)
%--------------------------- Calculate the shape of the surface in the tank
% If Omega = 0 (no rotation) the surface is flat.
if Omega == 0 
    % Values are saved to variable S (surface). S.z0 indicates the lowest
    % level of the water surface and S.a how the surface varies with the 
    % distance from the center of the tank.
    S.z0 = h0;
    S.a  = nan;
else
    Lx = 0;     % Tank Length
    Ly = 0;     % Tank Width
    g = 9.81;       % Gravity
    % When the tank is in rotation the water surface will deform as a
    % paraboloid. Following the description of "Newtons Bucket", the water
    % surface will follow the following height:
    % z = h_min + (Omega r)^2 / (2g),
    % where r = sqrt(x^2 +y^2) is the distance from the center of rotation
    % and h_min is the height of the water layer at the center. In this 
    % case, h_min is obtained by volume conservation, which yields for a 
    % rectangular tank:
    % h_min = h0 - (Omega^2 / 24g) (Lx^2 + Ly^2).
    S.z0 = h0 - (1/(24.*g)).*(Omega.^2).*(Lx.^2 + Lx.^2);
    S.a  = (1/(2.*g)) * Omega^2;

    % Convert values to mm
    S.z0 = S.z0 * 1e3;
    S.a  = S.a * 1e-3;
end

%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
%====================================== Construct virtual "lines" in space 
%================================== from the uninterrumpted projected dots 
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

% Set the pixel position of each of the projected dot patterns
xp = pattern.x;
yp = pattern.y; 

% Obtain the direction vector and origin from the projector calibration and
% dot pixel position.
pr_x0 = pr.cx0(xp, yp);
pr_y0 = pr.cy0(xp, yp);
pr_z0 = zeros(size(xp));
pr_dx = pr.cdx(xp, yp);
pr_dy = pr.cdy(xp, yp);
pr_dz = ones(size(xp));

% When the setup is rotating, the distance from the center of the tank is
% included in the calculations since the water surface deforms as a 
% paraboloid.
deltaX = 200;
deltaY = 500;

%=========================== Find where the "lines" in space intersect the
%======================== water surface. If the tank rotates, use function
%======================== [find_intersection_with_surface] to determine the 
%========================== intersection points with the parabolic surface.
%============================ Otherwise, extend lines until the flat water
%========================================================== surface height.

if Omega ~= 0 % Intersection points with the parabolic water surface 
    r0 = [pr_x0 - deltaX, pr_y0 - deltaY, pr_z0];
    dv = [pr_dx, pr_dy, pr_dz];
    [Ps, Normal] = find_intersection_with_surface(S, r0, dv);
    pr_xw = Ps.x + deltaX;
    pr_yw = Ps.y + deltaY;
    pr_zw = Ps.z;   
else % Intersection points with a flat water surface.
    z = zwater;
    pr_xw = pr_x0 + pr_dx*z;  
    pr_yw = pr_y0 + pr_dy*z; 
    pr_zw = pr_z0 + pr_dz*z;
end

%================= Calculate the direction vectors of the incoming "lines"
dvi = [pr_dx, pr_dy, pr_dz];
for i = 1:numel(xp) % Do line for line
    if Omega ~= 0
        nrm = Normal(i,:);
    else
        nrm = Normal; %[0,0,1];
    end
    pr_dvo(i,:) = snells_law(dvi(i,:), nrm, Nair, Nwater); %#ok<SAGROW>
end

%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
%====================================== Construct virtual "lines" in space 
%====================================== from the photographed dot patterns 
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

% First combine all the detected dots from images
x = []; 
y = []; 
for ii = 1:npos
    x = [x; image_dots(ii).x]; %#ok<AGROW>
    y = [y; image_dots(ii).y]; %#ok<AGROW>
end

% Obtain the direction vector and origin from the camera calibration and
% dot pixel position.
cm_x0 = cm.cx0(x, y);
cm_y0 = cm.cy0(x, y);
cm_z0 = zeros(size(x));
cm_dx = cm.cdx(x, y);
cm_dy = cm.cdy(x, y);
cm_dz = ones(size(x));

%=========================== Find where the "lines" in space intersect the
%======================== water surface. As before, depends on the water
%======================== surface shape.
if Omega ~= 0
    r0 = [cm_x0 - deltaX, cm_y0 - deltaY, cm_z0];
    dv = [cm_dx, cm_dy, cm_dz];
    [Ps, normal] = find_intersection_with_surface(S, r0, dv);

    cm_xw = Ps.x + deltaX;
    cm_yw = Ps.y + deltaY;
    cm_zw = Ps.z;
else
    z = zwater;
    cm_xw = cm_x0 + cm_dx*z; 
    cm_yw = cm_y0 + cm_dy*z; 
    cm_zw = cm_z0 + cm_dz*z;
end

%================= Calculate the direction vectors of the incoming "lines"
dvi = [cm_dx, cm_dy, cm_dz];
for i = 1:numel(x)
    
    if Omega ~= 0
        nrm = Normal(i,:);
    else
        nrm = Normal; %[0,0,1];
    end
    [cm_dvo(i,:),~] = snells_law(dvi(i,:), nrm, Nair, Nwater); %#ok<SAGROW>
    
end

%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
% =============================== Find the intersection between the virtual
% =============================== lines drawn by the dots in the projection
%================================ file and the lines obtained from the dots
%================================ in the photographs
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

%--- Use the function [lines_mindist] to determine the point in space where
%------- the minimum distance (or the intersection) between the constructed 
%------------------------------------------------ virtual lines is located. 

clear Pmid
for i = 1:numel(pr_x0)
    lin1.p0 = [pr_xw(i), pr_yw(i), pr_zw(i)];
    lin1.dv = pr_dvo(i,:);
    lin2.p0 = [cm_xw(i), cm_yw(i), cm_zw(i)];
    lin2.dv = cm_dvo(i,:);
    [Pmid(i,:), ~] = lines_mindist(lin1, lin2); %#ok<SAGROW>
end

%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
% =========================================== Mesh data obtained and save %
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%\

% Since the measurements obtained are not distributed uniformly, the data
% is interpolated into a equally spaced grid. The value of [Interp_Div] is 
% used to set the resolution of the interpolated grid.

%------------------------------------------------------ Morphology at start
output = Pmid;
% first remove not existing points
k = find(~isnan(output(:,1)));
x = output(k,1);
y = output(k,2);
v = output(k,3);

xrange = linspace(min(output(:,1)),max(output(:,1)),Interp_Div);
yrange = linspace(min(output(:,2)),max(output(:,2)),Interp_Div);
[xqS,yqS] = meshgrid(xrange, yrange);
vqS = griddata(x,y,v,xqS,yqS, 'cubic');

% save both the original measurements and the re-grided interpolation
save(Reconstruct,'Pmid','xqS','yqS','vqS');
