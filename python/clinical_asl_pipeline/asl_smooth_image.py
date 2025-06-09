import numpy as np
from scipy.ndimage import convolve, gaussian_filter
import warnings

def asl_smooth_image(data, spatialdim, FWHM, voxelsize):
    # Smooth an ASL image using a Gaussian filter.
    # This function handles NaN values in the input data and applies a Gaussian smoothing
    # filter based on the specified full width at half maximum (FWHM) and voxel size.    
    # The function supports both 2D and 3D data, and uses a NaN-aware convolution method
    # to ensure that NaN values do not affect the smoothing process.
    #  Parameters:
    # - data: input array
    # - spatialdim: 2 or 3
    # - FWHM: full width at half maximum in mm
    # - voxelsize: voxel size (list or array)
    # Returns:
    # - output: smoothed image, dtype=float
    
    def nanconvn(a, k, *args):
        """NaN-aware N-dimensional convolution, integrated from nanconvn.m."""
        edge = False
        nanout = False
        shape = 'same'

        for arg in args:
            arg = arg.lower()
            if arg == 'edge':
                edge = True
            elif arg == 'noedge':
                edge = False
            elif arg == 'nanout':
                nanout = True
            elif arg == 'nonanout':
                nanout = False
            elif arg == 'same':
                shape = 'same'
            else:
                raise ValueError(f"Unknown argument: {arg}")

        sza = a.shape
        o = np.ones_like(a)
        on = np.ones_like(a)
        n = np.isnan(a)

        # Replace NaNs with 0
        a_filled = np.copy(a)
        a_filled[n] = 0
        on[n] = 0

        # Check that kernel has no NaNs
        if np.isnan(k).any():
            raise ValueError("Filter (k) contains NaN values.")

        # Flat function convolution → correction factor
        if np.any(n) or edge:
            flat = convolve(on, k, mode='reflect')
        else:
            flat = o

        if np.any(n) and not edge:
            flat /= convolve(o, k, mode='reflect')

        with np.errstate(invalid='ignore', divide='ignore'):
            c = convolve(a_filled, k, mode='reflect') / flat
        c[flat == 0] = np.nan
        
        if nanout:
            c[n] = np.nan

        return c

    def fspecial_gaussian_2d(size, sigma):
        """Generate a 2D Gaussian kernel."""
        x = np.arange(-size // 2 + 1, size // 2 + 1)
        x, y = np.meshgrid(x, x)
        g = np.exp(-(x**2 + y**2) / (2 * sigma**2))
        return g / g.sum()

    def fspecial3_gaussian(size, sigma):
        """Generate a 3D Gaussian kernel (fspecial3 equivalent)."""
        x = np.arange(-size // 2 + 1, size // 2 + 1)
        x, y, z = np.meshgrid(x, x, x, indexing='ij')
        g = np.exp(-(x**2 + y**2 + z**2) / (2 * sigma**2))
        return g / g.sum()

    # ---- Start of main smoothing logic ----

    sigma = FWHM / 2.355
    inplanevoxelsize = voxelsize[0]
    data_smooth = np.zeros_like(data)

    if spatialdim == 2:
        print(f"Smoothing 2D - can handle NaNs: FWHM(mm) = {FWHM}")
        filtWidth = 7
        filtSigma = sigma / inplanevoxelsize
        imageFilter = fspecial_gaussian_2d(filtWidth, filtSigma)

        for s in range(data.shape[2]):
            if np.isnan(data).any():
                data_smooth[:, :, s] = nanconvn(data[:, :, s], imageFilter, 'nanout')
            else:
                warnings.warn("Smoothing 2D without NaN handling")
                data_smooth[:, :, s] = gaussian_filter(data[:, :, s], sigma=filtSigma, mode='reflect')

    elif spatialdim == 3:
        if np.isnan(data).any():
            print(f"Smoothing 3D - can handle NaNs: FWHM(mm) = {FWHM}")
            filtWidth = 7
            filtSigma = sigma / inplanevoxelsize
            imageFilter = fspecial3_gaussian(filtWidth, filtSigma)
            data_smooth = nanconvn(data, imageFilter, 'nanout')
        else:
            print(f"Smoothing 3D: FWHM(mm) = {FWHM}")
            data_smooth = gaussian_filter(data, sigma=sigma / np.array(voxelsize), mode='reflect')

    output = data_smooth.astype(float)
    return output
