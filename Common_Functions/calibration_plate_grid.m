%% #######################################################################%
%------------ Function describing the calibration plate geometry ---------%
%#########################################################################%
% Single source of truth for the physical dot grid on the calibration
% plate used by both the camera and projector calibration steps. The
% camera calibration builds its world-coordinate grid from these values,
% and the projector calibration uses the same values to work out which
% plate area the projected pattern could physically have covered (used to
% discard projector dots detected outside the plate). Keeping both derived
% from this one function means changing the physical plate only requires
% an edit here.

function plate = calibration_plate_grid()
    plate.nx = 23;    % number of dots in x
    plate.ny = 23;    % number of dots in y
    plate.dx = 25;    % dot spacing in x (mm)
    plate.dy = 25;    % dot spacing in y (mm)
    plate.center_x = 300; % x position of the grid centre (mm)
    plate.center_y = 300; % y position of the grid centre (mm)
end
