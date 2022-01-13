%% #######################################################################%
%------------ MATLAB Script used to measure submerged surfaces -----------%
%#########################################################################%
% This is an example script used to measure different surfaces. Seven cases
% are presented, three of an undulated plate (Base and High resolution and 
% under rotation), two of standing cylinders (mess and arranged) and two of
% a sediment bed. It requires the calibration files given by base camera
% and projector calibration and those obtained from the example calibration
% script. Each step taken here is briefly described. For a complete
% description consult [Gonzalez-Vera, A. S., Wilting, T. J. S., Holten,
% A. P. C., van Heijst, G. J. F. & Duran-Matute, M. "High-resolution
% single-camera photogrammetry: incorporation of refraction at a fluid 
% interface". Experiments in Fluids 61(1),3]. 
% Processed files are located in the directory: (Photogrammetry\Examples\
% Data_Files\MATLAB\Data_process_steps). The necesary photographed patterns
% for each case are located in the (Case_Photos\) directory. The resulting 
% measurement files are saved in 
% (Photogrammetry\Examples\Data_Files\MATLAB\Measurements).

addpath('..\..\Common_Functions\'); % Folder containing the functions used
                                    % in the script.
BaseLoc = '..\Data_Files\MATLAB\'; % Base directory to save data
StepSave = 'Data_process_steps\Measurements\'; % Extension to the directory 
                                               % where the data is saved to
                                               % place each partial step in
                                               % the data process.

zwater = 216; % Level, or depth, of the water layer (in mm).
theta = 0.0015; % A small roation of the pattern is implemented.
measurement = 7; % State the case that will be processed.
Omega = 0; % Rotation rate (rad/s) of the set-up. Examples shown here did
             % not rotate except for a single case (measurement = 3).
Normal = [0,0,1]; % Define the Normal vector (points upward)
load('refraction', 'Nwater', 'Nair'); % Load air and water index of 
                                      % refraction values
                                      
cm = load([BaseLoc 'Calibration_data\Camera_LinesFit_Base'],...
             'cx0', 'cy0', 'cdx', 'cdy');
pr = load([BaseLoc 'Calibration_data\Projector_LinesFit_Base'],...
             'cx0', 'cy0', 'cdx', 'cdy');
%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
%=================================================================== Cases:
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
switch measurement
%=========================================================================%
%--------------------------------- Underwater Undulated Plate Measurements%
%=========================================================================%
    case 1 %--------------------------------------------- "Base resolution"
        % Directory where pattern photographs are located
        imagedir    = 'Case_Photos\Undulated_Plate\Low_Resolution\';
        % Name of files where part of the image processing steps are saved
        camdotsfile = 'DotsInCam_Undulated_BaseR.mat';
        sortedfile  = 'SortedDots_Undulated_BaseR.mat';
        % Name of file where the result (reconstruction) is saved
        Reconstruct = 'Measurements\Reconstruct_Undulated_BaseRes.mat';
        % Pattern used for measurements
        Pattern_Fl  = 'Pattern_Base.mat';
        
        thrsh = 0; % Thereshold value that will be ignored when processing 
                   % the photographs of the projected dots.
        Rrange = [6,10]; % Range of radii used for the dot detection in the 
                         % images sent to the projector 
        DotR = [6 10]; % Range of radii used for the dot detection in the
                       % photographed patterns
        Interp_Div = 200; %Used for grid interpolation.
                     
    case 2 %--------------------------------------------- "High Resolution"
        imagedir    = 'Case_Photos\Undulated_Plate\High_Resolution\';
        camdotsfile = 'DotsInCam_Undulated_HighR.mat';
        sortedfile  = 'SortedDots_Undulated_HighR.mat';
        Reconstruct = 'Measurements\Reconstruct_Undulated_HighRes.mat';
        Pattern_Fl  = 'Pattern_High.mat';
        thrsh = 50;
        Rrange = [5,10];
        DotR = [3 14];
        Interp_Div = 550;
        cm = load([BaseLoc 'Calibration_data\Camera_LinesFit_Example'],...
             'cx0', 'cy0', 'cdx', 'cdy');
        pr=load([BaseLoc 'Calibration_data\Projector_LinesFit_Example'],...
             'cx0', 'cy0', 'cdx', 'cdy');
         

    case 3 %------------------------------------- In a rotating square tank
        imagedir    = 'Case_Photos\Undulated_Plate\With_Rotation\';
        camdotsfile = 'DotsInCam_Undulated_Rotating.mat';
        sortedfile  = 'SortedDots_Undulated_Rotating.mat';
        Reconstruct = 'Measurements\Reconstruct_Undulated_Rotating.mat';
        Pattern_Fl  = 'Pattern_Base.mat';
        thrsh = 0;
        Rrange = [5,10];
        DotR = [6 10];
        Interp_Div = 200;
        Omega = 1.496; %(rad/s)
        
%=========================================================================%
%------------------------------------------ Measurement of Three Cylinders%
%=========================================================================%
    case 4 %------------------------------------ Standing cylinders 1st try
        imagedir    = 'Case_Photos\Cylinders\cylinders_1\';
        camdotsfile = 'DotsInCam_Cylinders_1.mat';
        sortedfile  = 'SortedDots_Cylinders_1.mat';
        Reconstruct = 'Measurements\Reconstruct_Cylinders1_mess.mat';
        Pattern_Fl  = 'Pattern_Base.mat';
        thrsh = 0;
        Rrange = [5,10];
        DotR = [5 35];
        Interp_Div = 200;
        
    case 5 %------------------ Standing cylinders. Artificial dot selection
        imagedir    = 'Case_Photos\Cylinders\cylinders_1\';
        camdotsfile = 'DotsInCam_Cylinders_1.mat';
        sortedfile  = 'SortedDots_Cylinders_1_Artificial.mat';
        Reconstruct = 'Measurements\Reconstruct_Cylinders1_Artificial.mat';
        Pattern_Fl  = 'Pattern_Base.mat';
        thrsh = 0;
        Rrange = [5,10];
        DotR = [5 35];
        Interp_Div = 200;
        
%=========================================================================%
%----------------------------------------------- Sediment Bed Measurements%
%=========================================================================%
    case 6 %------------------------------------------- "Flat" Sediment Bed
        imagedir    = 'Case_Photos\Sediment\Flat\';
        camdotsfile = 'DotsInCam_Sediment_Flat.mat';
        sortedfile  = 'SortedDots_Sediment_Flat.mat';
        Reconstruct = 'Measurements\Reconstruct_Sediment_Flat.mat';
        Pattern_Fl  = 'Pattern_Base.mat';
        thrsh = 0;
        Rrange = [4,9];
        DotR = [6 15];
        Nwater = 1.33; % For the sediment measurments, a slight change in
                       % the refraction value is implemented.  
        Interp_Div = 200;                
        
    case 7 %---------------------- Sediment bed Distubed by a strong vortex
        imagedir    = 'Case_Photos\Sediment\Disturbed\';
        camdotsfile = 'DotsInCam_Sediment_VortexDisturbed.mat';
        sortedfile  = 'SortedDots_Sediment_VortexDisturbed.mat';
        Reconstruct = 'Measurements\Reconstruct_Sediment_Disturbed.mat';
        Pattern_Fl  = 'Pattern_Base.mat';
        thrsh = 0;  
        Rrange = [4,9];
        DotR = [6 15];
        Nwater = 1.33;
        Interp_Div = 200;
end

%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
%=================== Determine positions of the dots sent to the projector%
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

%-------------------------- Data file where the (pixel) position of the dot 
%-------------------------- patterns sent to the projector are saved.
PatternFile  = [BaseLoc 'Patterns\' Pattern_Fl];
infofile = [BaseLoc StepSave 'Projected_DotLocation_case' ...
    num2str(measurement) '.mat'];

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
infofile = [BaseLoc StepSave camdotsfile];

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
        camdots(ii) = find_circles(im,DotR);
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
infofile = [BaseLoc StepSave sortedfile];
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
        Ipat{ii} = PatternPrj.DotsToProject(:,:,ii);
        % Pattern captured by the camera.
        Icam{ii} = imread(image_names{ii});
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

h0 = zwater/1000;     % water layer height in metres (initial)
%--------------------------- Calculate the shape of the surface in the tank
% If Omega = 0 (no rotation) the surface is flat.
if Omega == 0 
    % Values are saved to variable S (surface). S.z0 indicates the lowest
    % level of the water surface and S.a how the surface varies with the 
    % distance from the center of the tank.
    S.z0 = h0;
    S.a  = nan;
else
    Lx = 1.000;     % Tank Length
    Ly = 1.000;     % Tank Width
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
    S.z0 = h0 - (1/(24.*g)).*(Omega.^2).*(Lx.^2 + Ly.^2);
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
deltaY = 500+5;

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
    pr_dvo(i,:) = snells_law(dvi(i,:), nrm, Nair, Nwater);
end

%% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%
%====================================== Construct virtual "lines" in space 
%====================================== from the photographed dot patterns 
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!%

% First combine all the detected dots from images
x = []; 
y = []; 
if measurement ~= 5
   for ii = 1:npos
    x = [x; image_dots(ii).x]; %#ok<AGROW>
    y = [y; image_dots(ii).y]; %#ok<AGROW>
    end 
else % The artificial selection already combined all dots 
    x = image_dots.x;
    y = image_dots.y;
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
    [Ps, Normal] = find_intersection_with_surface(S, r0, dv);

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
    [cm_dvo(i,:),~] = snells_law(dvi(i,:), nrm, Nair, Nwater);
    
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
[xq,yq] = meshgrid(xrange, yrange);
vq = griddata(x,y,v,xq,yq, 'cubic');

% save both the original measurements and the re-grided interpolation
save([BaseLoc Reconstruct],'Pmid','xq','yq','vq');

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%==================== Plot resulting reconstruction ======================%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------- Plot the obtained measurements(unevenly spaced)
figure(1)
scatter3(Pmid(:,1),Pmid(:,2),Pmid(:,3),15,'filled')
view(41,42);
daspect([1 1 .1])
set(gca, 'YDir','reverse')
axis([min(Pmid(:,1))-1 max(Pmid(:,1))+1 min(Pmid(:,2))-1 ...
    max(Pmid(:,2)+1)])
xlabel({'$x$ (mm)'},'FontUnits','points','interpreter','latex',...
        'FontSize',14,'FontName','Times')
ylabel({'$y$ (mm)'},'FontUnits','points','interpreter','latex',...
        'FontSize',14,'FontName','Times')

%------------------------ Plot the regrided measurements (uniformly spaced)
figure(2)
surf(xq, yq, vq,'EdgeColor','none');
daspect([1 1 .1])
grid on
set(gca, 'YDir','reverse')
axis([floor(min(min(xq))) ceil(max(max(xq))) floor(min(min(yq))) ...
    ceil(max(max(yq))) floor(min(min(vq))) ceil(max(max(vq)))])
xlabel({'$x$ (mm)'},'FontUnits','points','interpreter','latex',...
        'FontSize',14,'FontName','Times')
ylabel({'$y$ (mm)'},'FontUnits','points','interpreter','latex',...
        'FontSize',14,'FontName','Times')
zlabel({'$z$ (mm)'},'FontUnits','points','interpreter','latex',...
        'FontSize',14,'FontName','Times')
   