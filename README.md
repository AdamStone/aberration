aberration.py
=========

This script produces aberration-correction holograms for use with a digital 
spatial light modulator (SLM), for the purpose of achieving uniform focal 
conditions independent of focal depth during femtosecond laser irradiation
within transparent samples. Supports irradiation through an additional 
window layer, as in the case of a sample within a high temperature stage.


Dependencies
-------------

Python 2.x, matplotlib, numpy, and scipy.


Getting Started
-------------

Experimental parameters are specified in the inputs.py script, in the 
manner illustrated by the included example. Parameters are grouped 
into several classes and passed as instances to create_bitmaps().

A successful run will calculate the aberration correction patterns and
place them in a folder called 'bitmaps'. If a distortion correction 
file is specified, an additional set of patterns will be calculated which
include both aberration and distortion corrections. 

The distortion correction is used to compensate for phase distortion 
due to the SLM surface not being perfectly flat. An example 'distortion.bmp' 
pattern is included, but in practice each SLM module should have its own 
pattern provided by the manufacturer to be used here. The distortion-
corrected patterns should be used for actual experiments, while the pre-
distorted patterns are useful for aligning the SLM module. 

Output filenames are formatted as

    window{thickness}_depth{focal depth}_z{focal displacement}

Focal depth here refers to the target depth of the focus below the 
sample surface. Focal displacement refers to the stage displacement 
required after focusing on the surface in order to achieve the intended 
focal depth after the refraction and aberration correction is accounted 
for. Values are given in microns. For example:

    window1000_depth500_z263

This hologram would be used when there is a 1 mm window to focus through, 
and a 500 micron focal depth after correction is desired. The sample should 
be positioned by focusing on the sample surface, then raising the stage 
by 263 microns, in order to achieve a real focal depth of (approximately)
500 microns.


Implementation Details
------------

The algorithm presented by Itoh et al. (http://dx.doi.org/10.1364/OE.17.014367) 
was generalized to account for multiple refracting layers, in order to apply 
aberration correction while irradiating a sample inside a high-temperature 
stage through a silica window. Details of the derivation and experimental 
implementation are provided in 

    Adam Stone et al. "Multi-layer aberration correction for depth-independent 
    3D crystal growth in glass by femtosecond laser heating." JOSA B, Vol. 30, 
    Issue 5, pp. 1234-1240 (2013). <http://dx.doi.org/10.1364/JOSAB.30.001234>

Example Outputs
------------

An example aberration-correction pattern with slit beamshaping before applying distortion correction:

[Example aberration correction hologram](/example output/aberration_corrected/window1000_depth1200_z543.bmp)

After applying an additional distortion correction:

[Example aberration correction hologram with distortion correction](/example output/aberration_distortion_corrected/window1000_depth1200_z543.bmp)