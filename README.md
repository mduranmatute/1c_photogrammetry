# 1c_photogrammetry

MATLAB code to perform single-camera, single-projector photogrammetry, as described in:

González-Vera, A. S., Wilting, T. J. S., Holten, A. P. C., van Heijst, G. J. F., & Duran-Matute, M. (2020). High-resolution single-camera photogrammetry: incorporation of refraction at a fluid interface. *Experiments in Fluids*, 61(1), 3.

A digital projector projects a pattern of dots onto the surface to be measured; a single camera photographs it. Comparing the pixel positions of the projected dots to their photographed positions, together with a calibration of both the camera and the projector, yields the 3D shape of the surface - including, optionally, a flat or parabolic fluid interface between the camera and the surface (e.g. a water layer) and the refraction it introduces.

For the full derivation and method, see the paper above.

## Requirements

- MATLAB. Developed against MATLAB 2015b; also verified working on MATLAB R2024b.
- **Image Processing Toolbox** (used for dot detection via `imfindcircles` and thresholding via `im2bw`).
- The Curve Fitting Toolbox is **not** required - the polynomial calibration fits are implemented directly in `Common_Functions/fit_poly.m` and `Common_Functions/eval_poly.m`.

## Repository structure

| Path | Contents |
|---|---|
| `Camera_Calibration_Template.m` | Template script to calibrate the camera. Start here to calibrate a new camera. |
| `Projector_Calibration_Template.m` | Template script to calibrate the projector. Requires a finished camera calibration. |
| `Surface_calculation_template.m` | Template script to measure a surface. Requires finished camera and projector calibrations. |
| `Common_Functions/` | Shared functions used by the templates and examples (dot detection, calibration-plate geometry, polynomial fitting, refraction, line intersection, interactive point-matching). |
| `Examples/` | Filled-in versions of the three templates above, showing complete worked calibrations and measurements (undulated plate, standing cylinders, sediment bed). These scripts expect a photographs/data directory tree that is **not** included in this repository (see below). |
| `NetCDF_Convert/` | Scripts to convert `.mat` files produced by the pipeline into NetCDF for long-term archival (`Calib_C.m`, `Patterns_C.m`, `Measurements_C.m`, `StepsM_C.m`), and the reverse - `Calib_NC2MAT.m`, `Patterns_NC2MAT.m`, `Measurements_NC2MAT.m` reconstruct working `.mat` calibration, pattern, and reconstruction files from an archived NetCDF dataset, for someone who only has the published NetCDF data and wants to run the pipeline against it. |
| `LICENSE` | Apache License 2.0. |

The `Examples/` scripts reference photographs and processed `.mat` files (calibration plate photos, projected pattern photos, dot positions, reconstructions) that belong to the experiments described in the paper. That data is not part of this repository; the scripts are included to show a complete, working pipeline and as a starting point for adapting to your own data.

If you only have a published NetCDF archive (no `.mat` files), running the `*_NC2MAT.m` scripts in `NetCDF_Convert/` reconstructs working calibration, pattern, and reconstruction `.mat` files you can use directly as `Cam_File`/`Proj_File`/`Pattern_Fl` inputs elsewhere in the pipeline. One caveat: a calibration NetCDF file archived before `Plens` and `zlevels` were added to `Calib_C.m` won't have them, so `Calib_NC2MAT.m` reconstructs those two as empty for such a file (with a warning) - the calibration coefficients themselves are unaffected either way.

## Usage

The pipeline has three stages, each with its own template script. Run them in order:

1. **Camera calibration** - `Camera_Calibration_Template.m`. Fill in the blank variables near the top (photograph directory, calibration plate heights, where to save results, etc.) before running; each required field is checked with an `assert` that explains what it's for. Produces a camera calibration `.mat` file.
2. **Projector calibration** - `Projector_Calibration_Template.m`. Requires the camera calibration file from step 1. Produces a projector calibration `.mat` file.
3. **Surface measurement** - `Surface_calculation_template.m`. Requires both calibration files from steps 1 and 2, plus photographs of the projected pattern on the surface to be measured. Produces the reconstructed surface.

A few things worth knowing before running any of these:

- Some steps require interactive input: you'll be asked to click four corresponding points in two images to relate a projected/plate pattern to what the camera photographed (`Common_Functions/sort_camera_calibration_data.m` and `Common_Functions/sort_projected_dots.m`). These steps need a real display and cannot be run headlessly.
- Intermediate results (detected dot positions, sorted dot positions) are cached to the `.mat` files you specify. On a re-run, you'll be asked whether to reuse the cached file or redo that step - useful when iterating on a later step without redetecting dots from scratch.
- The calibration plate's physical geometry (grid size, dot spacing, position) is defined once in `Common_Functions/calibration_plate_grid.m`. If you use a different calibration plate, update it there rather than in the calibration scripts.
- `Examples/` contains complete, filled-in versions of all three templates, worth reading alongside the templates even without the underlying photo data.

## Method overview

1. Upload the file containing the patterns that are projected.
2. Find the center of the dots that will be projected in each pattern. The result is given in pixels.
3. Upload the directories that contain the photographs of the patterns that were projected on the surface to be measured.
4. Find the dots in the photographs and establish their pixel position.
5. Sort which dots belong to the projected patterns; eliminate wrong or undetected dots.
6. Calculate the location of the interface, if present.
7. Transform the pixel positions of the projected and photographed dots into vectors with an origin `r0` and direction `dv`. These are the virtual lines in space, obtained via the camera and projector calibrations.
8. Find the intersection of each virtual line with the interface, if present, extending the virtual lines until the intersection position.
9. Calculate the directional vectors of the incoming (refracted) lines.
10. Find where the lines from the projected dots and their respective photographs intersect (or the point of minimum distance between them). Each intersection is a measurement point.
11. Remove non-existent intersections (`NaN`s).
12. Because the measurements are not uniformly distributed, cubic interpolation redefines the result onto a regular grid.

## Data

The photographs and processed data behind the paper's example cases (undulated plate, standing cylinders, sediment bed - see `Examples/`) were obtained from experiments performed by the authors and described in the article above. The results in `Examples/` were originally obtained with MATLAB 2015b.

## License

Apache License 2.0 - see [LICENSE](LICENSE).
