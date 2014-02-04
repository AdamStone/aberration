""" Example implementation of aberration and distortion correction, including
irradiation through window and simulated beam shaping. """

from aberration import ObjectiveLens, SLM, Sample, Window, Slit, create_bitmaps

# lens setup
lens = ObjectiveLens(NA = 0.55,                 # numerical aperture         
                     f = 4e-3,                  # focal length, in meters  
                     wavelength = 800e-9)       # laser wavelength, in meters

# SLM setup
slm = SLM(px = 20e-6,                           # pixel size, in meters
          nx = 792,                             # number of pixels in x direction
          ny = 600,                             # number of pixels in y direction
          m_tel = 2.,                           # magnification of internal telescope
          correction_file = 'distortion.bmp')   # filename of distortion correction bmp (if available)

# sample setup
sample = Sample(n_s = 2.06,                     # refractive index of the sample
                focal_depths = [ 100e-6,        # focal depths to calculate, in meters
                                 200e-6, 
                                 500e-6, 
                                 1000e-6, 
                                 1100e-6, 
                                 1200e-6 ])

# window setup (if applicable)
window = Window(n_w = 1.399,                    # refractive index of window
                thicknesses = [0, 1000e-6])     # window thicknesses, in meters

# slit beamshaping setup (if applicable)                                     
slit = Slit(width = 0.5,                        # as a fraction of total pattern width 
            orientation = 'y')                  # 'x' or 'y'


# pass configuration objects to create_bitmaps
create_bitmaps(lens, slm, sample, window, slit)