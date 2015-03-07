"""
This package produces aberration-correction holograms for use with a digital
spatial light modulator (SLM), for the purpose of achieving uniform focal
conditions independent of focal depth during femtosecond laser irradiation
within transparent samples. Supports irradiation through an additional
window layer, as in the case of a sample within a high temperature stage.

Details of the derivation and experimental implementation are provided in

    Adam Stone et al. "Multi-layer aberration correction for depth-independent
    3D crystal growth in glass by femtosecond laser heating." JOSA B, Vol. 30,
    Issue 5, pp. 1234-1240 (2013). <http://dx.doi.org/10.1364/JOSAB.30.001234>

The interface module contains those items which might be imported in a typical
usage, including classes to store experimental configuration information and
the create_bitmaps function which is used to initiate the calculations. These
will be imported into the aberration namespace when 'import aberration' is
called from within the parent directory.
"""

from __future__ import division
import numpy as np
import os

from core import optimize_d_f, calculate_phasemap, save


class ObjectiveLens(object):
    """
    Convenience class for passing objective lens parameters to create_bitmaps.
    """
    def __init__(self, NA, f, wavelength, n_a=1):
        """
        Parameters
        ----------

        NA : float
            Numerical aperture of the objective lens.

        f : float
            Focal length of objective lens (in m).

        wavelength : float
            Wavelength of incident laser (in m).

        n_a : float
            Refractive index of the ambient.
        """
        self.NA = NA
        self.f = f
        self.wavelength = wavelength
        self.n_a = n_a
        self.alpha = np.arcsin(self.NA/self.n_a)


class SLM(object):
    """
    Convenience class for passing SLM parameters to create_bitmaps.
    """
    def __init__(self, px, nx, ny, m_tel=2, correction_file=None, scale=254):
        """
        Parameters
        ----------

        px : float
            SLM pixel size (in m).

        nx : int
            Number of SLM pixels in x direction.

        ny : int
            Number of SLM pixels in y direction.

        m_tel : float
            Magnification of the internal telescope within the SLM module,
            relative to the imaging direction (rays from the focus are
            magnified by m_tel before reaching the SLM screen; the laser
            traveling from the SLM screen to the focus is demagnified by
            1/m_tel).

        correction_file : str
            Filename of distortion correction hologram to be added to the
            aberration correction hologram to to account for phase distortion
            caused by the SLM surface not being perfectly flat.

        scale : float
            Rescales the hologram such that the range of phase shifts
            between 0 and 2pi can be tuned within the 256 channels of the
            bitmap. The value expected is the channel at which the phase
            shift should correspond to 2pi, which may vary
            (consult SLM documentation).
        """
        self.px = px
        self.nx = nx
        self.ny = ny
        self.m_tel = m_tel
        self.correction_file = correction_file
        self.scale = scale


class Sample(object):
    """
    Convenience class for passing sample parameters to create_bitmaps.
    """
    def __init__(self, n_s, focal_depths):
        """
        Parameters
        ----------

        n_s : float
            Refractive index of the sample.

        focal_depths : list
            List of intended focal depth values (in m) below the
            sample surface. For example, to irradiate at focal depths of
            250 micron intervals with aberration correction applied,
            use [250e-6, 500e-6, 750e-6, 1000e-6].
        """
        self.n_s = n_s
        self.focal_depths = focal_depths


class Window(object):
    """
    Convenience class for passing window parameters to create_bitmaps.
    """
    def __init__(self, n_w, thicknesses):
        """
        Parameters
        ----------

        n_w : float
            Refractive index of the window.

        thicknesses : list
            List of window thickness values (in m) to accommodate. For example,
            in order to calculate phase maps for both the cases of no window
            and a 1 mm window, use [0, 1e-3].
        """
        self.n_w = n_w
        self.thicknesses = thicknesses


class Slit(object):
    """
    Convenience class for passing slit beamshaping
    parameters to create_bitmaps.
    """
    def __init__(self, width, orientation='y'):
        """
        Parameters
        ----------

        width : float
            Slit width, as a fraction of total beam diameter.

        orientation : str
            Sets slit to horizontal or vertical orientation.
            Expects 'x' or 'y'.
        """

        self.width = width
        self.orientation = orientation


def create_bitmaps(lens, slm, sample, window=None, slit=None):
    """ The main function of the script, which initiates the calculation
    of hologram bitmaps based on the input conditions.

    Parameters
    ----------

    lens : `ObjectiveLens` instance
        Carries details about objective lens and laser parameters.

    slm : `SLM` instance
        Carries details about SLM parameters.

    sample : `Sample` instance
        Carries details about sample and focal depths of interest.

    window : `Window` instance
        Carries details about window, if applicable.

    slit : `Slit` instance
        Carries details about slit beamshaping, if applicable.
    """

    if not window:
        window = Window(1, [0])

    for d_w in window.thicknesses:
        for d_s in sample.focal_depths:

            if d_w != 0:
                print ('Processing window thickness {} micron and sample '
                       'depth {} micron...').format(d_w*1e6, d_s*1e6)
            else:
                print 'Processing sample depth {} micron...'.format(d_s*1e6)

            d_f = optimize_d_f(lens.wavelength, d_s, lens.f, lens.alpha,
                               sample.n_s, lens.n_a, d_w, window.n_w)
            d_shift = -optimize_d_f(lens.wavelength, 0, lens.f, lens.alpha,
                                    sample.n_s, lens.n_a, d_w, window.n_w)
            z = d_shift + d_f

            ab_corrected, ab_dist_corrected = calculate_phasemap(
                slm.px, slm.nx, slm.ny, lens.wavelength, d_s, lens.f,
                lens.alpha, sample.n_s, lens.n_a, d_w, window.n_w, slm.m_tel,
                slm.correction_file, slm.scale, slit)

            if d_w != 0:
                filename = 'window{}_depth{}_z{}'.format(
                    int(d_w*1e6), int(d_s*1e6), int(z*1e6))
            else:
                filename = 'depth{}_z{}'.format(int(d_s*1e6), int(z*1e6))

            try:
                save('bitmaps/aberration_corrected/'+filename, ab_corrected)
            except:
                os.makedirs('bitmaps/aberration_corrected')
                save('bitmaps/aberration_corrected/'+filename, ab_corrected)

            if ab_dist_corrected is not None:
                try:
                    save('bitmaps/aberration_distortion_corrected/'+filename,
                         ab_dist_corrected)
                except:
                    os.makedirs('bitmaps/aberration_distortion_corrected')
                    save('bitmaps/aberration_distortion_corrected/'+filename,
                         ab_dist_corrected)
