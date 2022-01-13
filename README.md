# 1c_photogrammetry
This repository contains the code to perform photogrammetry using one camera and one digital projector as described by:

González-Vera, A. S., Wilting, T. J. S., Holten, A. P. C., van Heijst, G. J. F., & Duran-Matute, M. (2020). High-resolution single-camera photogrammetry: incorporation of refraction at a fluid interface. Experiments in Fluids, 61(1), 1-19.


The following steps are a simple description of how the code functions.
For a in-depth explanation see: [Gonzalez-Vera,A. S., Wilting,T. J. S.,
Holten,A. P. C., van Heijst,G. J. F. & Duran-Matute,M. "High-resolution
single-camera photogrammetry: incorporation of refraction at a fluid 
interface". Experiments in Fluids 61(1),3]. 

01.- Upload file containing the patterns that are projected.

02.- Find the center of the dots that will be projected in each of the
     the patterns. result is given in pixels. 
     
03.- Upload the directories that contain the photographs of the patterns
     that were projected on the surface to be measured.
     
04.- Find dots in the photographs and establish their pixel position.

05.- Sort which dots belong to the projected patterns. Eliminate 
     wrong or not found dots.
     
06.- Calculate the location of the interface if present.

07.- Transform the pixel position of the projected and photographed dots
     into vectors with an origin r0 and direction dv. These can be 
     considered as the virtual lines in space. This requires the
     interpolation obtained from the calibration of the camera and 
     projector. 
     
08.- Find the intersection of the each of the virtual lines with the
     interface if it is present. Extend the virtual lines until the 
     intersection positions.
     
09.- Calculate the directional vectors from the incoming lines.

10.- Find where the lines from the projected dots and their respective
     photographs intersect (or find the location of minimum distance).
     Each of the intersections is a measurement point.
     
11.- Remove non-existent intersections (NaNs).

12.- Due to the non-uniform distribution of the measurements, a cubic 
     interpolation is used to redefine the result into a grid.

Example scripts for a few test case measurements (and their respective 
camera and projector calibration procedure and files) can be found in 
the (Example\) directory. The resulting processed files and results 
were obtained with MATLAB 2015b.

The data was obtained from experiments performed by the authors and 
described in the article: 
Gonzalez-Vera, A. S., Wilting, T. J. S., Holten, A. P. C., van Heijst,
G. J. F. & Duran-Matute, M. "High-resolution single-camera photogrammetry:
 incorporation of refraction at a fluid interface". Experiments in Fluids 
61(1), 3.
