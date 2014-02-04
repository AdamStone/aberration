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
    
The core module includes utility functions and functions which perform
the calculations. In most cases, a user will not need to import any of
these functions directly.
"""

from __future__ import division
import numpy as np
import scipy.misc as misc
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.figure import Figure
import random
import PIL
import os
    
    
def save(filename, np_array, png=False): 
    """ saves numpy array as bitmap and/or png image. A list of details
    may be included to save relevant info about the images in a CSV file
    of the same name. By default only a bitmap is produced. 
    
    Parameters
    ----------
    
    filename : string
        Filename (without extension). 
        
    np_array : numpy.array
        2D numpy array (dtype=uint8) representing grayscale image to save.
        
    png : bool
        By default, a bmp file is produced. If png=True, a png file will be 
        retained along with the bitmap.
    """

    value_min=0
    value_max=255
    color_map='gray'
    
    fig = Figure(figsize=np_array.shape[::-1], dpi=1, frameon=False)
    fig.canvas = FigureCanvas(fig)
    fig.figimage(np_array, cmap=color_map, vmin=value_min, vmax=value_max, origin=None)
    fig.savefig(filename, dpi=1, format=None)
    
    # save bitmap (convert from PNG)
    img = PIL.Image.open(filename + '.png')
    if len(img.split()) == 4:           # prevent IOError: cannot write mode RGBA as BMP
        r, g, b, a = img.split()
        img = PIL.Image.merge("RGB", (r, g, b))      
    greyscale = img.convert("L")        # L indicates PIL's greyscale mode
    greyscale.save(filename + '.bmp')
            
    # delete PNG
    if not png:
        os.remove(filename + '.png')


def diff_theta(theta_1, d_f, d_s, f, n_s, n_a=1, d_w=0, n_w=1): 
    """ The derivative of theta (incident angle of a ray in the absence of 
    refraction) with respect to theta_1 (the corresponding convergence angle 
    after correction is applied). 
    
    Theta_1 as a function of theta is desired, but this cannot be solved 
    analytically. An iterative method such as Newton's Method can 
    be applied to obtain an approximation of theta_1 as a function of theta, 
    and Newton's Method requires this derivative. The analytical expression 
    was obtained by differentiating the expression for theta using the SAGE 
    mathematics suite (http://www.sagemath.org).

    Parameters
    ----------
    
    theta_1 : float
        Incident angle (in rads) after correction corresponding to a given theta.
        
    d_f : float
        Depth (in m) of hypothetical focal point below sample surface 
        in the absence of refraction (i.e. focal length minus the distance 
        between the objective lens and the sample surface).
        
    d_s : float
        Intended focal depth (in m) below sample surface (after correction).
        
    f : float
        Focal length of objective lens (in m).
    
    n_s : float
        Refractive index of the sample.
        
    n_a : float
        Refractive index of the ambient.

    d_w : float
        Thickness of the window (in m), used when irradiating through a window
        (e.g. sample inside a heated stage with viewport). Use 0 for no window.
        
    n_w : float
        Refractive index of the window.
        
    Returns
    ----------
    
    ret : float
        Value of derivative of theta with respect to theta_1.
    """
    
    return ((d_w*n_a*np.sin(theta_1)/np.sqrt(-n_a**2*np.sin(theta_1)**2 + n_w**2) + 
            d_s*n_a*np.sin(theta_1)/np.sqrt(-n_a**2*np.sin(theta_1)**2 + n_s**2) - 
            (d_f + d_w)*np.tan(theta_1))*np.sin(theta_1)/f - 
            (d_w*n_a**3*np.sin(theta_1)**2*np.cos(theta_1)/(-n_a**2*np.sin(theta_1)**2 + n_w**2)**(3./2) + 
            d_s*n_a**3*np.sin(theta_1)**2*np.cos(theta_1)/(-n_a**2*np.sin(theta_1)**2 + n_s**2)**(3./2) + 
            d_w*n_a*np.cos(theta_1)/np.sqrt(-n_a**2*np.sin(theta_1)**2 + n_w**2) + 
            d_s*n_a*np.cos(theta_1)/np.sqrt(-n_a**2*np.sin(theta_1)**2 + n_s**2) - (np.tan(theta_1)**2 + 
            1)*(d_f + d_w))*np.cos(theta_1)/f)/np.sqrt(-(d_w*n_a*np.sin(theta_1)/np.sqrt(-n_a**2*np.sin(theta_1)**2 + 
            n_w**2) + d_s*n_a*np.sin(theta_1)/np.sqrt(-n_a**2*np.sin(theta_1)**2 + n_s**2) - 
            (d_f + d_w)*np.tan(theta_1))**2*np.cos(theta_1)**2/f**2 + 1) - 1


def newton_theta(theta, d_f, d_s, f, n_s, n_a=1, d_w=0, n_w=1):
    """ Newton's Method to find theta_1 for a given theta. 
    
    Theta_1 as a function of theta is desired, but cannot be solved 
    analytically. Newton's Method yields a numerical approximation.

    Parameters
    ----------
    
    theta : float
        Incident angle (in rads) of a ray in the absence of refraction; 
        effectively specifies a radial position on the objective lens.
        
    d_f : float
        Depth (in m) of hypothetical focal point below sample surface 
        in the absence of refraction (i.e. focal length minus the distance 
        between the objective lens and the sample surface).
        
    d_s : float
        Intended focal depth (in m) below sample surface (after correction).
        
    f : float
        Focal length of objective lens (in m).

    n_s : float
        Refractive index of the sample.
        
    n_a : float
        Refractive index of the ambient.

    d_w : float
        Thickness of the window (in m), used when irradiating through a window
        (e.g. sample inside a heated stage with viewport). Use 0 for no window.
        
    n_w : float
        Refractive index of the window.
    
    Returns
    ----------
    
    ret : float
        Numerical approximation of the value of theta_1 for a given theta.
    """
    
    def g_th1(theta_1, d_s, d_w): #separating out a large term from F
        return d_s*n_a*np.sin(theta_1)/np.sqrt(n_s**2 - n_a**2*np.sin(theta_1)**2) + \
                d_w*n_a*np.sin(theta_1)/np.sqrt(n_w**2 - n_a**2*np.sin(theta_1)**2)
    
    iterations = 0
    theta_1 = theta # initial 'guess'
    F = theta - (theta_1 + np.arcsin((np.cos(theta_1)/f)*(g_th1(theta_1, d_s, d_w) - (d_f+d_w)*np.tan(theta_1))))
    while abs(F) > 1e-10: # precision limit
        F = theta - (theta_1 + np.arcsin((np.cos(theta_1)/f)*(g_th1(theta_1, d_s, d_w) - (d_f+d_w)*np.tan(theta_1))))
        theta_1 -= F/diff_theta(theta_1, d_f, d_s, f, n_s, n_a, d_w, n_w)
        iterations += 1
    return theta_1
    
    
def PHI(theta, d_s, d_f, f, n_s, n_a=1, d_w=0, n_w=1):
    """ Calculates total optical path length (geometric distance of path multiplied 
    by refractive index of the medium) from a point on the lens (specified by theta) to 
    the intended focal point (after correction).
    
    Parameters
    ----------
    
    theta : float
        Incident angle (in rads) of a ray in the absence of refraction; 
        effectively specifies a radial position on the objective lens.
        
    d_s : float
        Intended focal depth (in m) below sample surface (after correction).
        
    d_f : float
        Depth (in m) of hypothetical focal point below sample surface 
        in the absence of refraction (i.e. focal length minus the distance 
        between the objective lens and the sample surface).
        
    f : float
        Focal length of objective lens (in m).
    
    n_s : float
        Refractive index of the sample.
        
    n_a : float
        Refractive index of the ambient.    

    d_w : float
        Thickness of the window (in m), used when irradiating through a window
        (e.g. sample inside a heated stage with viewport). Use 0 for no window.
        
    n_w : float
        Refractive index of the window.
        
    Returns
    ----------
    
    ret : float
        Value of optical path length.
    """
    
    def DE(d_s, theta_2):                   # length of ray DE
        return d_s/np.cos(theta_2)
    
    def BC(d_w, theta_3):                   # length of ray BC
        return d_w/np.cos(theta_3)
        
    def ABpCD(theta, theta_1, d_f, f, d_w):   # length of rays AB+CD
        return (f*np.cos(theta)-d_f-d_w)/np.cos(theta_1)
        
    if theta != 0:
        theta_1 = newton_theta(theta, d_f, d_s, f, n_s, n_a, d_w, n_w)
        theta_2 = np.arcsin(n_a*np.sin(theta_1)/n_s)
        theta_3 = np.arcsin(n_a*np.sin(theta_1)/n_w)
        ans = n_a*ABpCD(theta, theta_1, d_f, f, d_w) + n_s * DE(d_s, theta_2) + n_w*BC(d_w, theta_3)
    else:
        ans = n_a*(f-d_f-d_w) + n_s*d_s + n_w*d_w 
    return ans
    
    
def DeltaPHI(wavelength, theta, d_s, d_f, f, n_s, n_a=1, d_w=0, n_w=1, wrap=False): 
    """ Calculates the optical path length difference between PHI(theta) and 
    PHI(theta=0). This is equivalent to the phase shift required to bring the 
    ray at incident angle theta into phase with the ray along the beam axis
    at the intended focal point.
    
    Parameters
    ----------
    
    wavelength : float
        Wavelength of incident laser (in m).
    
    theta : float
        Incident angle (in rads) of a ray in the absence of refraction; 
        effectively specifies a radial position on the objective lens.
        
    d_s : float
        Intended focal depth (in m) below sample surface (after correction).
        
    d_f : float
        Depth (in m) of hypothetical focal point below sample surface 
        in the absence of refraction (i.e. focal length minus the distance 
        between the objective lens and the sample surface).
        
    f : float
        Focal length of objective lens (in m).
    
    n_s : float
        Refractive index of the sample.

    n_a : float
        Refractive index of the ambient.    

    d_w : float
        Thickness of the window (in m), used when irradiating through a window
        (e.g. sample inside a heated stage with viewport). Use 0 for no window.
        
    n_w : float
        Refractive index of the window.
        
    wrap : bool
        If True, phase shift values are constrained between 0 and 2pi.
        
    Returns
    ----------
    
    ret : float
        Value of the optical path length difference (in radians) for the given theta.
    """
    
    ans = -( PHI(theta, d_s, d_f, f, n_s, n_a, d_w, n_w) - PHI(0, d_s, d_f, f, n_s, n_a, d_w, n_w) )
    ans = ans/wavelength*2*np.pi 
    if wrap:
        ans %= 2*np.pi
    return ans
    
    
def optimize_d_f(wavelength, d_s, f, alpha, n_s, n_a=1, d_w=0, n_w=1):
    """ The value of d_f for an intended d_s could be chosen arbitrarily, 
    but a choice may be considered optimal which minimizes the peak-
    to-valley (PV) value of the phase pattern (i.e. the maximum phase shift 
    required to bring the most disparate rays into phase at the intended 
    focus). 
    
    This function finds the value of d_f for a given d_s which minimizes the 
    PV value by iteratively pinning DeltaPHI = 0 where theta = alpha. 

    Parameters
    ----------
    
    wavelength : float
        Wavelength of incident laser (in m).
    
    d_s : float
        Intended focal depth (in m) below sample surface (after correction).
        
    f : float
        Focal length of objective lens (in m).
        
    alpha : float
        Maximum convergence angle of the lens (in rad), as determined by 
        alpha = asin(NA/n_a).
    
    n_s : float
        Refractive index of the sample.

    n_a : float
        Refractive index of the ambient.    

    d_w : float
        Thickness of the window (in m), used when irradiating through a window
        (e.g. sample inside a heated stage with viewport). Use 0 for no window.
        
    n_w : float
        Refractive index of the window.
        
    Returns
    ----------
    
    ret : float
        The value of d_f optimized for a given d_s.
    """
    
    d_f = d_s # initial guess
    ans = DeltaPHI(wavelength, alpha, d_s, d_f, f, n_s, n_a, d_w, n_w)
    step = d_f/2. + 1e-6
    if ans <= 0:
        flip = 1
    else:
        flip = -1
    while abs(ans) > d_s*1e-5: # precision threshold
        while flip*ans <= 0:
            d_f -= flip*step
            ans = DeltaPHI(wavelength, alpha, d_s, d_f, f, n_s, n_a, d_w, n_w)
        d_f += flip*step
        ans = DeltaPHI(wavelength, alpha, d_s, d_f, f, n_s, n_a, d_w, n_w)
        step /= 2.
    return d_f


def calculate_phasemap(px, nx, ny, wavelength, d_s, f, alpha, n_s, n_a=1, d_w=0, n_w=1, 
                       m_tel=2, correction_file=None, scale=254, slit=None, wrap=True): 
    """ Calculates the aberration correction hologram for a given focal condition, 
    for use with a digital SLM. 
    
    If correction_file is specified (e.g. the filename of a distortion 
    correction map), an additional hologram will be produced that incorporates
    the additional correction.
    
    Slit beamshaping can be simulated by providing a `Slit` instance.    
    
    Parameters
    ----------
    
    px : float
        SLM pixel size (in m).
    
    nx : int
        Number of SLM pixels in x direction.
    
    ny : int
        Number of SLM pixels in y direction.    
    
    wavelength : float
        Wavelength of incident laser (in m).
    
    d_s : float
        Intended focal depth (in m) below sample surface (after correction).
        
    f : float
        Focal length of objective lens (in m).
        
    alpha : float
        Maximum convergence angle of the lens (in rad), as determined by 
        alpha = asin(NA/n_a).
    
    n_s : float
        Refractive index of the sample.
        
    n_a : float
        Refractive index of the ambient.

    d_w : float
        Thickness of the window (in m), used when irradiating through a window
        (e.g. sample inside a heated stage with viewport). Use 0 for no window.
        
    n_w : float
        Refractive index of the window.
    
    m_tel : float
        Magnification of the internal telescope within the SLM module, relative
        to the imaging direction (rays from the focus are magnified by m_tel
        before reaching the SLM screen; the laser traveling from the SLM screen 
        to the focus is demagnified by 1/m_tel).

    correction_file : str
        Filename of distortion correction hologram to be added to the aberration
        correction hologram to to account for phase distortion caused by the SLM 
        surface not being perfectly flat. 
        
    scale : float
        Rescales the hologram such that the range of phase shifts 
        between 0 and 2pi can be tuned within the 256 channels of the 
        bitmap. The value expected is the channel at which the phase shift 
        should correspond to 2pi. 
        
    slit : `Slit` instance
        If provided, some fraction the sides of the pattern will be cut off by 
        replacing those pixels with random phase values, in order to simulate 
        the effect of passing the beam through a slit.
        
    wrap : bool
        If True, phase shift values are constrained between 0 and 2pi.
        
    Returns
    ----------
    
    final : 'array'        
        2D numpy array of 8-bit integers containing the aberration correction 
        phase shifts, with the range of 0-2pi represented by the integer range 
        0-255. Intended for saving as a grayscale bitmap.
        
    corrected : 'array'
        2D numpy array of 8-bit integers containing the combined aberration and 
        distortion correction phase shifts, with the range of 0-2pi represented 
        by the integer range 0-255. Intended for saving as a grayscale bitmap.
    """
    
    d_f = optimize_d_f(wavelength, d_s, f, alpha, n_s, n_a, d_w, n_w)
    h_max = m_tel*f*np.sin(alpha)
    px_range = int(h_max/px)
    
    phase_NW = np.zeros((px_range, px_range))
    
    #diagonal:
    for i in range(px_range):
        d_axial = (px_range - i - 0.5)*px
        r_SLM = np.sqrt(2*(d_axial)**2)     
        if r_SLM <= h_max:
            theta_r = np.arcsin(r_SLM/m_tel/f)
            phase_NW[i,i] = DeltaPHI(wavelength, theta_r, d_s, d_f, f, n_s, 
                                     n_a, d_w, n_w, wrap=wrap)
        else:
            phase_NW[i,i] = random.random()*2*np.pi

        
    #corner:
    for j in range(px_range): #y
        d_y = (px_range - j - 0.5)*px
        for i in range(1+j, px_range): #x
            d_x = (px_range - i - 0.5)*px
            r_SLM = np.sqrt((d_x)**2 + (d_y)**2)     
            if r_SLM <= h_max:
                theta_r = np.arcsin(r_SLM/m_tel/f)
                phase_NW[j,i] = DeltaPHI(wavelength, theta_r, d_s, d_f, f, 
                                         n_s, n_a, d_w, n_w, wrap=wrap)
            else:
                phase_NW[j,i] = random.random()*2*np.pi
            phase_NW[i,j] = phase_NW[j,i]


    #copy to other quadrants:
    NE = phase_NW.copy()
    NE = np.fliplr(NE)
    NW_NE = np.hstack((phase_NW, NE))
    
    SW_SE = NW_NE.copy()
    SW_SE = np.flipud(SW_SE)
    
    NSWE = np.vstack((NW_NE, SW_SE))
    
    NSWE = np.array(np.round(NSWE*255/(2*np.pi)), dtype='uint8')
    img = PIL.Image.fromarray(NSWE)
    basewidth = int(img.size[0])
    hpercent = 46.5/48 # account for slight anisotropy from mirror tilt
    hsize = int(round((float(img.size[1]) * float(hpercent))))
    if hsize - int(hsize/2)*2 != 0:
        hsize += 1
    img = img.resize((basewidth, hsize), PIL.Image.NEAREST)
    
    NSWE = np.array(img, dtype='float')
    NSWE = np.array(NSWE/255*(2*np.pi))
    
    #insert into full-size grid for SLM
    if not slit:
        active_width = len(NSWE[0])
        active_height = len(NSWE)
        hmargin = (nx - active_width)/2
        vmargin = (ny - active_height)/2
            
        upper = np.random.rand(vmargin, nx)*2*np.pi
        left_mid = np.random.rand(active_height, hmargin)*2*np.pi
        right_mid = np.random.rand(active_height, hmargin)*2*np.pi
        middle = np.hstack((left_mid, NSWE, right_mid))
        lower = np.random.rand(vmargin, nx)*2*np.pi
    
    else:
        if slit.orientation == 'x':
            active_width = len(NSWE[0])
            active_height = int(len(NSWE)*slit.width)
            if np.mod(active_height, 2):           # if odd
                active_height += 1
            slit_margin = (len(NSWE) - active_height)/2
        if slit.orientation == 'y':
            active_width = int(len(NSWE[0])*slit.width)
            active_height = len(NSWE)        
            if np.mod(active_width, 2):            # if odd
                active_width += 1
            slit_margin = (len(NSWE[0]) - active_width)/2
                    
        hmargin = (nx - active_width)/2
        vmargin = (ny - active_height)/2
        
        upper = np.random.rand(vmargin, nx)*2*np.pi
        left_mid = np.random.rand(active_height, hmargin)*2*np.pi
        right_mid = np.random.rand(active_height, hmargin)*2*np.pi
        lower = np.random.rand(vmargin, nx)*2*np.pi
        if slit.orientation == 'x':
            middle = np.hstack((left_mid,NSWE[slit_margin:-slit_margin],right_mid))
        if slit.orientation == 'y':
            middle = np.hstack((left_mid,NSWE[:, slit_margin:-slit_margin],right_mid))
    
    final = np.vstack((upper, middle, lower))
    final = np.array(np.round(final*255/(2*np.pi)), dtype='uint8')
    
    if correction_file: 
        correction = misc.imread(correction_file)
        correction = np.array(np.round(correction[:,:,0]), dtype='uint8')
        corrected = np.array(np.round((scale/256.)*(final+correction)), dtype='uint8')
    else:
        corrected = None
        
    final = np.array(np.round((scale/256)*(final)), dtype='uint8')
    
    return final, corrected