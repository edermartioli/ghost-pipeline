# GHOST Pipeline

A tool-kit to reduce Gemini High-resolution Optical SpecTrograph (GHOST)
data and perform a cross-correlation (CCF) analysis to obtain precise radial
velocities.

The GHOST Pipeline wraps the [Gemini DRAGONS](https://dragons.readthedocs.io/)
reduction of raw GHOST exposures and provides a set of command-line tools to
build 1D spectrum (`s1d`) products, construct high signal-to-noise template
spectra, compute CCFs with weighted stellar line masks, and measure radial
velocities, bisector spans, FWHMs, and activity indicators. It also includes
tools for differential RV analysis with dual-IFU observations and for Doppler
tomography of the Rossiter–McLaughlin effect.

## Installation

Clone the repository and install locally with pip:

```bash
git clone https://github.com/edermartioli/ghost-pipeline.git
cd ghost-pipeline
pip install -U .
```

This installs the `ghostpipe` Python package, the standard CCF masks, and
the command-line tools listed below. The core dependencies (`numpy`,
`scipy`, `matplotlib`, `astropy`, `scikit-learn`, `uncertainties`) are
installed automatically.

Optional extras:

```bash
pip install -U ".[transit]"   # adds batman-package, used by ghost_dual_timeseries
```

**Note:** the raw-data reduction step (`ghost_pipeline`) requires
[DRAGONS](https://dragons.readthedocs.io/) with the GHOST data-reduction
package to be installed and available in the same environment (the
`reduce`, `dataselect`, `caldb`, and `showd` commands). See the
[GHOST DR tutorial](https://ghost-drtutorial.readthedocs.io/) for
instructions. All other tools work on reduced products and do not require
DRAGONS.

## Command-line tools

| Tool | Description |
|---|---|
| `ghost_pipeline` | Run the DRAGONS recipes to reduce raw GHOST exposures (bias, flat, arc, standard, science) |
| `ghost_log` | Print an observation log for a set of raw GHOST FITS files |
| `ghost_build_s1d` | Stitch echelle orders of `*_calibrated.fits` products into `s1d` 1D spectra, computing BJD and BERV |
| `ghost_template_s1d` | Build a template spectrum out of a time series of `s1d` spectra |
| `ghost_ccf_pipeline` | CCF radial-velocity analysis of a time series of `s1d` spectra (RV, bisector, and FWHM time series) |
| `ghost_get_simple_template` | Inspect blue+red template products and export a combined simple template |
| `ghost_dual_timeseries` | Photometric flux-ratio time series from dual-IFU observations |
| `ghost_diff_rv_timeseries` | Differential RV time series between a target and a simultaneous comparison |
| `ghost_rm_tomography` | Doppler tomography of the Rossiter–McLaughlin effect from a CCF time series |
| `ghost_plot_dragons` | Quick-look plot of GHOST spectra in DRAGONS format |
| `sun_ccf_pipeline` | CCF analysis of a set of IAG Solar Atlas spectra |

All tools accept `-h`/`--help` for the full list of options. They can also
be run as Python modules, e.g. `python -m ghostpipe.scripts.ghost_ccf_pipeline --help`.

## Typical workflow

### 1. Reduce the raw data with DRAGONS

First set up the DRAGONS calibration database. In `~/.geminidr/`, create or
edit the configuration file `rsys.cfg` with the following content:

```
[calibs]
standalone = True
database_dir = <path_to_my_data>
```

Then run the reduction. For example, to reduce a standard star and then a
science target:

```bash
ghost_pipeline --input=/path/to/rawdata/*.fits \
               --reduced_data_dir=/path/to/playground \
               --object="CD -32 9927" -s -bfao1r

ghost_pipeline --input=/path/to/rawdata/*.fits \
               --reduced_data_dir=/path/to/playground \
               --object="XX Oph" \
               --red_std_spec=S20230416S0073_red001_standard.fits \
               --blue_std_spec=S20230416S0073_blue001_standard.fits \
               -o1
```

The switches select the reduction steps: `-b` bias, `-f` flat, `-a` arc,
`-o` object, `-1` 1D spectra, `-s` object is a standard, `-r` reset the
calibration database. **Attention:** `-r` runs `caldb init --wipe`, which
erases previous information in the calibration database.

An equivalent step-by-step shell workflow using the DRAGONS command-line
utilities directly is provided in [`reduce.sh`](reduce.sh) for reference.

### 2. Build the s1d products

```bash
ghost_build_s1d --input="/path/to/*_calibrated.fits"
```

Add `-2` if the object was also observed with IFU2, `-s` if a sky fiber is
present, and `-b` to apply the barycentric correction to the wavelengths.

**Note:** the BJD and BERV are computed with the JPL DE430 ephemeris; the
first run will download the ephemeris file (~115 MB) through astropy, which
is then cached locally.

### 3. Run the CCF analysis to obtain radial velocities

```bash
ghost_ccf_pipeline --ccf_mask=G2_nm.mas --nknots=40 \
    --input="*red*_s1d.fits" \
    --output_ccfs_file=target_ghost_red_ccfs.fits \
    --output_rv_file=target_ghost_red_ccfrv.rdb \
    --output_bis_file=target_ghost_red_ccfbis.rdb \
    --output_fwhm_file=target_ghost_red_ccffwhm.rdb \
    --output_obslog_file=target_ghost_red_log.txt -pv
```

RV, bisector, and FWHM time series are saved in `rdb` format, and the CCFs
in a FITS product. The blue and red arms are analyzed independently
(e.g. run again with `--input="*blue*_s1d.fits"`).

**Note:** the CCF RV analysis is designed for a *time series* and requires
at least two exposures, since the statistical weights of the mask lines are
derived from the dispersion of the residuals around the template. To
measure the RV of a single spectrum, use the CCF machinery directly from
Python (see below).

### 4. Build a template spectrum (optional)

Once a first RV solution is available, a rest-frame template can be built
and used to iterate the analysis:

```bash
ghost_template_s1d --nknots=40 --input="*red*_s1d.fits" \
    --rv_file=target_ghost_red_ccfrv.rdb \
    --output=target_ghost_red_template.fits -v
```

### 5. Further analysis (optional)

* `ghost_diff_rv_timeseries` corrects the target RVs by the drift measured
  on a simultaneous comparison star (dual-IFU mode).
* `ghost_dual_timeseries` produces flux-ratio (differential photometry)
  time series from dual-IFU observations.
* `ghost_rm_tomography` maps the Doppler shadow of a transiting planet from
  the time series of CCF residuals.

## CCF masks

Standard stellar line masks are bundled with the package (in
`ghostpipe/masks/`) for spectral types G2, K0, K5, and M2, in Angstrom
(`G2.mas`, ...) and in nanometer (`G2_nm.mas`, ...) versions, along with a
few target-specific masks. From Python, their paths can be obtained with:

```python
import ghostpipe
mask = ghostpipe.get_mask_path("G2_nm")
```

## Test data

The [`data/`](data/) directory contains a minimal GHOST dataset of the RV
standard star HD 65907 (blue and red `*_calibrated.fits` products reduced
with DRAGONS, and the corresponding `s1d` products) that can be used to try
out the tools. For example, to rebuild the s1d products from the calibrated
files:

```bash
ghost_build_s1d --input="data/*_calibrated.fits"
```

Since this dataset contains a single exposure per arm, it is suited for
testing the s1d format and the single-spectrum CCF measurement (see the
Python example below); the full RV time-series analysis with
`ghost_ccf_pipeline` requires two or more exposures.

## Using the libraries from Python

The package can also be used as a library:

```python
from ghostpipe import ghostlib, reduc_lib, ccf_lib, spectrallib

spectrum = ghostlib.load_spectrum("S20260208S0016_red001_HD65907_s1d.fits")
```

For example, to measure the RV of a single spectrum with the CCF:

```python
import numpy as np
import ghostpipe
from ghostpipe import reduc_lib, ccf_lib

arr = reduc_lib.load_array_of_ghost_spectra(
    ["data/S20260208S0016_red001_HD65907_s1d.fits"], apply_berv=True)
sp = arr["spectra"][0]

params = ccf_lib.set_ccf_params(ghostpipe.get_mask_path("G2_nm"))
params["CCF_WIDTH"] = 150.0
weights = np.ones_like(sp["flux"])
mask = ccf_lib.apply_weights_to_ccf_mask(
    params, sp["wl"], sp["flux"], sp["fluxerr"], weights)
ccf = ccf_lib.run_ccf_eder(params, sp["wl"], sp["flux"],
                           sp["header"], mask, normalize_ccfs=True)
print("RV = {:.3f} km/s".format(ccf["header"]["RV_OBJ"]))
# -> RV = 14.877 km/s for HD 65907
```

The main modules are:

* `ghostpipe.ghostlib` — core utilities for handling GHOST spectra (order
  limits, s1d FITS I/O, BJD/BERV calculation);
* `ghostpipe.reduc_lib` — reduction of time series of s1d spectra
  (resampling, alignment, template construction, normalization);
* `ghostpipe.ccf_lib` — CCF machinery for radial-velocity measurements
  (adapted from the APERO/SPIRou CCF codes);
* `ghostpipe.spectrallib` — spectral quantities and stellar activity
  indicators (S-index, log R'HK, H-alpha, etc.).

## Contact

Eder Martioli — <emartioli@lna.br>

Laboratório Nacional de Astrofísica (LNA/MCTI), Brazil
Institut d'Astrophysique de Paris, France

## License

This project is distributed under the terms of the
[GNU General Public License v3.0](LICENSE).

The CCF routines in `ghostpipe/ccf_lib.py` were adapted from the CCF codes
of the APERO data reduction software for SPIRou.
