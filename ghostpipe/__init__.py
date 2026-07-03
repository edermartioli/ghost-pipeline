# -*- coding: utf-8 -*-
"""
GHOST Pipeline: reduction and radial-velocity analysis of Gemini/GHOST spectra.

This package provides a tool-kit to reduce Gemini High-resolution Optical
SpecTrograph (GHOST) data with DRAGONS, build 1D spectrum (s1d) and template
products, and perform a cross-correlation (CCF) analysis to obtain precise
radial velocities, bisector spans, FWHMs, and activity indicators.

Modules
-------
ghostlib
    Core utilities for handling GHOST spectra (order limits, s1d FITS I/O,
    barycentric corrections).
reduc_lib
    Spectral reduction utilities for time series of s1d spectra (resampling,
    alignment, template construction, normalization).
ccf_lib
    Cross-correlation function (CCF) toolkit for radial-velocity
    measurements (adapted from the APERO/SPIRou CCF codes).
spectrallib
    Spectral quantities and stellar activity indicators (S-index, log R'HK,
    H-alpha, etc.).

License
-------
Distributed under the terms of the GNU General Public License v3.0. See the
LICENSE file at the root of this project or
<https://www.gnu.org/licenses/gpl-3.0.html> for details.
"""

import os

__version__ = "1.0.0"

__all__ = ["ghostlib", "ccf_lib", "reduc_lib", "spectrallib",
           "masks_dir", "get_mask_path"]


def masks_dir():
    """Return the absolute path of the directory with the bundled CCF masks.

    Returns
    -------
    str
        Path of the ``masks`` directory installed with the package,
        containing the standard stellar line masks (G2, K0, K5, M2, in
        Angstrom and in nm versions) and a few target-specific masks.
    """
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "masks")


def get_mask_path(name):
    """Return the absolute path of a bundled CCF mask by name.

    Parameters
    ----------
    name : str
        Mask file name (e.g. ``"G2_nm.mas"``) or mask name without the
        extension (e.g. ``"G2_nm"``).

    Returns
    -------
    str
        Absolute path of the mask file.

    Raises
    ------
    FileNotFoundError
        If no mask with the given name is found.
    """
    if not name.endswith(".mas"):
        name = name + ".mas"
    path = os.path.join(masks_dir(), name)
    if not os.path.exists(path):
        available = sorted(f for f in os.listdir(masks_dir())
                           if f.endswith(".mas"))
        raise FileNotFoundError(
            "CCF mask '{}' not found. Available masks: {}".format(
                name, ", ".join(available)))
    return path
