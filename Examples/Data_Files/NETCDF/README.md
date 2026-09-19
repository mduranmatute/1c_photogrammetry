# Dataset description: NetCDF example data for single-camera photogrammetry

This document describes the collection of NetCDF (`.nc`) files in this
directory, together with the companion `Case_Photos.zip` archive of raw
photographs (see Section 11.5), for the purpose of archiving them in a FAIR
data repository. It accompanies, but is independent of, the
`1c_photogrammetry` code repository.

## 1. Title

Example data for "High-resolution single-camera photogrammetry: incorporation
of refraction at a fluid interface"

## 2. Description / abstract

This dataset contains example input, intermediate, and output data files for
the single-camera photogrammetry method described in the publication below.
The method reconstructs the 3D shape of a (possibly submerged) surface from
photographs of a sequence of dot patterns projected onto it with a digital
projector, using a two-step calibration of the camera and the projector.

The data are organized into four groups that follow the processing pipeline
of the method, from the patterns that are projected, through the intermediate
dot-detection and dot-sorting steps, to the final reconstructed surfaces, for
several example measurement cases (rigid cylinders, a sediment bed, and an
undulated plate, with and without rotation and at two spatial resolutions).
Camera and projector calibration coefficients are included separately. All
files are provided in classic NetCDF format and are self-describing (each
variable carries `long_name`, `units`, and `description` attributes, and each
file carries global attributes with authorship and provenance information).

Corresponding data are also available in MATLAB `.mat` format in
`Examples/Data_Files/MATLAB/` of the source code repository; both formats
contain the same content.

The raw camera photographs from which the `DotsInCam_*.nc` files were
derived are included as a separate ZIP archive, `Case_Photos.zip` (see
Section 11.5), so that the full pipeline can be reproduced from the original
images rather than only from the already-processed dot positions.

## 3. Related publication

González-Vera, A. S., Wilting, T. J. S., Holten, A. P. C., van Heijst,
G. J. F., & Duran-Matute, M. (2020). High-resolution single-camera
photogrammetry: incorporation of refraction at a fluid interface.
*Experiments in Fluids*, 61(1), 3.

DOI: https://doi.org/10.1007/s00348-019-2826-y

## 4. Related software

The code used to generate and process this data is available at:
https://github.com/mduranmatute/1c_photogrammetry

## 5. Creators / authors

| Name | Affiliation |
|---|---|
| González-Vera, A. S. | Eindhoven University of Technology, Department of Applied Physics |
| Wilting, T. J. S. | Eindhoven University of Technology, Department of Applied Physics |
| Holten, A. P. C. | Eindhoven University of Technology, Department of Applied Physics |
| van Heijst, G. J. F. | Eindhoven University of Technology, Department of Applied Physics |
| Duran-Matute, M. | Eindhoven University of Technology, Department of Applied Physics |

Contact / correspondence: m.duran.matute@tue.nl

## 6. Source

Laboratory experiments performed by the authors, Eindhoven University of
Technology, Department of Applied Physics.

## 7. Keywords

photogrammetry, structured light projection, 3D surface reconstruction,
camera calibration, refraction, fluid interface, digital projector, NetCDF

## 8. License

CC BY 4.0 (Creative Commons Attribution 4.0 International). Reuse is
permitted provided the source (the publication above) is credited. This
applies to the *data*, including the NetCDF files and the `Case_Photos.zip`
raw-photograph archive; the accompanying processing code is released under
the Apache License 2.0 (see the `LICENSE` file of the code repository).

## 9. Format and conventions

- Format: NetCDF classic (CDF-1), readable with any NetCDF library
  (`netCDF4`/`xarray` in Python, `ncread`/`ncdisp` in MATLAB, `ncdump` from
  the Unix command line, `RNetCDF` in R).
- Every variable has `long_name`, `units`, and `description` attributes.
  `units = "[]"` denotes a dimensionless quantity (e.g. an index or a
  polynomial coefficient with mixed physical units already captured
  elsewhere).
- The following global attributes are present in every file:

  | Attribute | Meaning |
  |---|---|
  | `Title` | Dataset/publication title |
  | `Description` | File-specific description, including provenance |
  | `Authors` | Author list (see Section 5) |
  | `email for correspondence` | Contact e-mail |
  | `Institution` | Author institution |
  | `Source` | Data source (laboratory experiments) |
  | `Creation date` | Date/time the file was generated (MATLAB `datestr` format) |

- Missing/non-existent values (e.g. dots that were not detected) are encoded
  as `NaN` in floating-point variables.
- Total collection size: ~377 MB across 35 NetCDF files, plus a ~610 MB
  companion `Case_Photos.zip` archive of raw PNG photographs (Section 11.5).

## 10. Directory structure and file naming

```
Examples/
├── Data_Files/
│   └── NETCDF/                      (this directory)
│       ├── README.md                this file
│       ├── Calibration_Data/        camera & projector calibration coefficients
│       ├── Patterns/                projected dot-pattern images
│       ├── Data_process_steps/
│       │   └── Measurements/        intermediate dot-detection/sorting results
│       └── Measurements/            final reconstructed 3D surfaces
└── Measurements/
    └── Case_Photos.zip              raw photographs (companion archive, Section 11.5)
```

Filenames encode the measurement case:

| Case label | Meaning |
|---|---|
| `Cylinders`, `Cylinders_1` | Rigid submerged cylinders test case |
| `Sediment_Flat` | Flat sediment bed |
| `Sediment_VortexDisturbed` / `Sediment_Disturbed` | Sediment bed disturbed by a vortex |
| `Undulated_BaseR` / `Undulated_BaseRes` | Undulated plate, base (low) spatial resolution |
| `Undulated_HighR` / `Undulated_HighRes` | Undulated plate, high spatial resolution |
| `Undulated_Rotating` | Undulated plate under rotation |
| `_Artificial` / `_mess` suffix | Alternative/auxiliary reconstruction of the same case (e.g. artificially generated or unfiltered ("messy") point set) |
| `_Base` / `_Example` / `_High` suffix (Patterns, Calibration) | Variant of the pattern or calibration used (baseline vs. worked example vs. high-resolution pattern set) |

## 11. File and variable dictionary

Files within each group share an identical variable structure; only the
dimension *lengths* differ between cases (they depend on the number of dots
detected/measured). NaN-padding is used where the number of dots differs
between patterns within the same file.

### 11.1 `Calibration_Data/` — 4 files

`Camera_LinesFit_Base.nc`, `Camera_LinesFit_Example.nc`,
`Projector_LinesFit_Base.nc`, `Projector_LinesFit_Example.nc`

Coefficients of the 3rd-order polynomial (`Poly33`) fit obtained from
calibration, mapping a pixel position `(x, y)` to the origin and direction of
the corresponding line of sight in 3D space:
`val(x,y) = p00 + p10*x + p01*y + p20*x^2 + p11*x*y + p02*y^2 + p30*x^3 + p21*x^2*y + p12*x*y^2 + p03*y^3`

| Variable | Dimensions | Units | Description |
|---|---|---|---|
| `Num_Coeff` | `Num_Coeff` (=10) | – | Coefficient index `[p00, p10, p01, p20, p11, p02, p30, p21, p12, p03]` |
| `cdx` | `Num_Coeff` | – | Coefficients of the line-direction vector, X component |
| `cdy` | `Num_Coeff` | – | Coefficients of the line-direction vector, Y component |
| `cx0` | `Num_Coeff` | mm | Coefficients of the line-origin vector, X component |
| `cy0` | `Num_Coeff` | mm | Coefficients of the line-origin vector, Y component |
| `Plens` | `Lens_dim` (=3) | mm | Approximate position (x, y, z) of the camera or projector lens determined during calibration |
| `zlevels` | `Num_Levels` (case-dependent) | mm | Heights at which the calibration plate or projected pattern was photographed. Absent from `Camera_LinesFit_Base.nc`, whose source calibration predates this value being routinely recorded |

### 11.2 `Patterns/` — 3 files

`Pattern_Base.nc`, `Pattern_High.nc`, `CalibrationPattern_Example.nc`

Binary (logical) images of the dot patterns projected onto the calibration
target / measured surface.

| Variable | Dimensions | Units | Description |
|---|---|---|---|
| `Num_Pattern` | `Num_Pattern` (=120) | – | Pattern index |
| `X_p` | `X_p` (=1366) | pixels | Horizontal resolution of the projector |
| `Y_p` | `Y_p` (=768) | pixels | Vertical resolution of the projector |
| `Pattern_Images` | `Num_Pattern` × `X_p` × `Y_p` | pixels (0/1) | Binary images of the projected patterns |

### 11.3 `Data_process_steps/Measurements/` — 21 files

Intermediate results of the processing pipeline (see the code repository
`README.md`, steps 02–05), for each measurement case.

**`DotsInCam_*.nc`** (7 files) — pixel positions of dots detected in the
photographs of the projected patterns:

| Variable | Dimensions | Units | Description |
|---|---|---|---|
| `N_Patterns` | `Num_Patterns` (=120) | – | Pattern index |
| `Max_Dots` | `Max_Ind` (case-dependent) | – | Maximum number of dots detected across all photographs of this case |
| `X_pos` | `Max_Ind` × `Num_Patterns` | pixels | X pixel position of each detected dot (NaN = not found) |
| `Y_pos` | `Max_Ind` × `Num_Patterns` | pixels | Y pixel position of each detected dot (NaN = not found) |

**`Projected_DotLocation_case1..7.nc`** (7 files) — same variable structure
as above, but for the pixel positions of dots in the projected pattern
images themselves (rather than the photographs).

**`SortedDots_*.nc`** (7 files) — same variable structure as above, but for
dots after they have been matched (sorted) between the projected patterns
and the corresponding photographs, with unmatched/invalid dots removed.

### 11.4 `Measurements/` — 7 files

`Reconstruct_Cylinders1_Artificial.nc`, `Reconstruct_Cylinders1_mess.nc`,
`Reconstruct_Sediment_Disturbed.nc`, `Reconstruct_Sediment_Flat.nc`,
`Reconstruct_Undulated_BaseRes.nc`, `Reconstruct_Undulated_HighRes.nc`,
`Reconstruct_Undulated_Rotating.nc`

Final output of the pipeline: the reconstructed 3D measurement points and the
interpolated surface elevation grid.

| Variable | Dimensions | Units | Description |
|---|---|---|---|
| `Num_Points` | `Num_Points` (case-dependent) | – | Measured-point index |
| `X_int` | `X_int` (case-dependent) | mm | Horizontal coordinate of the interpolation grid |
| `Y_int` | `Y_int` (case-dependent) | mm | Vertical coordinate of the interpolation grid |
| `Measurements` | 3 × `Num_Points` | mm | 3D coordinates `(x, y, z)` of each measured point, stacked as `x:(1:Num_Points,1)`, `y:(1:Num_Points,2)`, `z:(1:Num_Points,3)` |
| `Surface` | `Y_int` × `X_int` | mm | Surface elevation interpolated onto the regular `(X_int, Y_int)` grid |

### 11.5 `Case_Photos.zip` — companion raw-photograph archive

Location: `Examples/Measurements/Case_Photos.zip` (a sibling of this
`Data_Files/NETCDF/` directory, not inside it). Format: ZIP, containing 720
uncompressed-inside PNG photographs (~610 MB total, ~611 MB unpacked)
organized by measurement case:

```
Case_Photos/
├── Cylinders/
│   └── cylinders_1/                 120 photographs
├── Sediment/
│   ├── Flat/                        120 photographs
│   └── Disturbed/                   120 photographs
└── Undulated_Plate/
    ├── Low_Resolution/              120 photographs
    ├── High_Resolution/             120 photographs
    └── With_Rotation/               120 photographs
```

These are the original camera photographs of each of the 120 projected dot
patterns, for each measurement case, taken with the surface (cylinders,
sediment bed, or undulated plate) present. They are the raw input to
pipeline step 04 ("find dots in the photographs"; see Section 12) and are
each case's rawest form of data: `DotsInCam_<case>.nc` (Section 11.3) is
the pixel-position output of running dot detection on these images, so the
photographs are not needed to reproduce the reconstruction from
`DotsInCam_*.nc` onward, but they are needed to reproduce or re-run the
dot-detection step itself, or to apply alternative image-processing methods.

The case subdirectory names correspond to the case labels in Section 10 as
follows:

| `Case_Photos/` subdirectory | Corresponding `DotsInCam_*.nc` file |
|---|---|
| `Cylinders/cylinders_1/` | `DotsInCam_Cylinders.nc`, `DotsInCam_Cylinders_1.nc` |
| `Sediment/Flat/` | `DotsInCam_Sediment_Flat.nc` |
| `Sediment/Disturbed/` | `DotsInCam_Sediment_VortexDisturbed.nc` |
| `Undulated_Plate/Low_Resolution/` | `DotsInCam_Undulated_BaseR.nc` |
| `Undulated_Plate/High_Resolution/` | `DotsInCam_Undulated_HighR.nc` |
| `Undulated_Plate/With_Rotation/` | `DotsInCam_Undulated_Rotating.nc` |

The example script `Examples/Measure_underwater.m` of the code repository
reads photographs from this same relative path structure
(`Case_Photos/<case>/...`).

## 12. Methodology

The data were produced following the processing pipeline described in the
code repository `README.md` (pattern upload → dot detection in patterns →
photograph upload → dot detection in photographs → dot sorting → interface
location → pixel-to-line conversion via calibration → line/interface
intersection → line/line intersection to obtain measurement points → removal
of non-existent intersections → cubic interpolation onto a regular grid),
using MATLAB 2015b. Full methodological detail is given in the related
publication (Section 3).

## 13. How to cite

González-Vera, A. S., Wilting, T. J. S., Holten, A. P. C., van Heijst,
G. J. F., & Duran-Matute, M. (2020). High-resolution single-camera
photogrammetry: incorporation of refraction at a fluid interface.
*Experiments in Fluids*, 61(1), 3.

## 14. Version history

| Version | Date | Change |
|---|---|---|
| 1.0 | 2020-11-10 | Initial data files generated (see per-file `Creation date` global attribute) |
| 1.1 | 2026-09-17 | `Calibration_Data/*.nc` regenerated to also include `Plens` and `zlevels` (Section 11.1), previously omitted from the NetCDF conversion though always present in the corresponding `.mat` files |
