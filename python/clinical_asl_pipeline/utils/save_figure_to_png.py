"""
ClinicalASL - Clinical Arterial Spin Labeling processing pipeline

Figure saving utility module.

Repository: https://github.com/JSIERO/ClinicalASL

Author: Jeroen Siero
Institution: UMCU (University Medical Center Utrecht), The Netherlands
Contact: j.c.w.siero@umcutrecht.nl

Description:
    Utility function to save 3D data montages as PNG images with custom colormaps and labels.

License: BSD 3-Clause License
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from scipy.io import loadmat

def save_figure_to_png(data, mask, datarange, outputloc, suffix, label, colormap='viridis'):
    # Save a 3D data montage to PNG with black background and white labels.
    #
    # Args:
    #     data (np.ndarray): 3D array (H x W x S) of image data to visualize.
    #     mask (np.ndarray or None): 3D mask array (same shape as data), or None for no mask.
    #     datarange (tuple): (min, max) values for color scaling.
    #     outputloc (str): Directory to save the PNG file.
    #     suffix (str): Suffix for the output filename.
    #     label (str): Label for the colorbar.
    #     colormap (str): Name of colormap ('viridis', 'jet', 'vik', 'devon', etc.).
    #
    # The function:
    #   - Loads a colormap (including custom .mat colormaps if needed).
    #   - Sets the background color (zero values) to black.
    #   - Applies a mask if provided.
    #   - Builds a montage of all slices in the 3D data.
    #   - Scales and clamps data to the specified range.
    #   - Converts the montage to an RGB image using the colormap.
    #   - Plots the montage with a horizontal colorbar (white labels/ticks).
    #   - Saves the figure as a PNG with a black background.
    # -----------------------------------------------------------------------------
    current_dir = os.path.dirname(os.path.abspath(__file__))
    
    # Load colormap
    if colormap.lower() == 'viridis':
        cmap = plt.get_cmap('viridis')
    elif colormap.lower() == 'jet':
        cmap = plt.get_cmap('jet')
    elif colormap.lower() in ['vik', 'devon']:
        matfile = os.path.join(current_dir, f'{colormap.lower()}.mat')
        mat = loadmat(matfile)
        cmap_data = mat.get(colormap.lower())
        if cmap_data is None:
            raise ValueError(f"Could not find colormap data in {matfile}.")
        cmap = ListedColormap(cmap_data)
    else:
        raise ValueError(f"Unknown colormap: {colormap}")

    # Set first color to black (for background/zero values)
    cmap_array = cmap(np.linspace(0, 1, 256))
    cmap_array[0] = [0, 0, 0, 1]
    cmap = ListedColormap(cmap_array)

    # Label mapping
    label_dict = {
        'CBF': 'ml/100g/min',
        'CVR': 'ΔCBF ml/100g/min',
        'BOLDCVR': 'BOLD CVR (%)',
        'time': 'time (s)',
        'time_delta': 'Δtime (s)',
        'a.u.': 'a.u.',
        '%': '%'
    }
    label2 = label_dict.get(label, label)

    if mask is None:
        mask = np.ones_like(data)

    # Mask and clean data
    data_masked = data * mask
    data_masked[np.isnan(data_masked)] = 0  # ensure any NaNs become 0
    
    # Build montage
    ncols = 4
    n_slices = data.shape[2]
    nrows = int(np.ceil(n_slices / ncols))
    slice_h, slice_w = data.shape[:2]
    montage = np.zeros((nrows * slice_h, ncols * slice_w))

    for i in range(n_slices):
        sl = np.rot90(data_masked[:, :, i])
        r, c = divmod(i, ncols)
        montage[r*slice_h:(r+1)*slice_h, c*slice_w:(c+1)*slice_w] = sl

    # Clamp and scale
    montage = np.clip(montage, datarange[0], datarange[1])
    # Create scaled array, but preserve zeros as-is
    scaled = np.full_like(montage, 0)  # default: all black
    nonzero_mask = montage != 0

    # Scale only nonzero values
    scaled[nonzero_mask] = 255 * (montage[nonzero_mask] - datarange[0]) / (datarange[1] - datarange[0])
    scaled = np.clip(scaled, 0, 255).astype(np.uint8)
    
    # Convert to RGB
    rgb_img = cmap(scaled)
    rgb_img = (rgb_img[:, :, :3] * 255).astype(np.uint8)

    # Prepare figure
    h_px, w_px = rgb_img.shape[:2]
    dpi = 100  # or 150 if you want sharper output
    fig = plt.figure(figsize=(w_px / dpi, h_px / dpi), dpi=dpi)
    ax = fig.add_axes([0, 0.1, 1, 0.9])
    ax.imshow(rgb_img)
    ax.axis('off')
    fig.patch.set_facecolor('black')

    # Colorbar
    cbar_ax = fig.add_axes([0.2, 0.03, 0.6, 0.03])
    norm = plt.Normalize(datarange[0], datarange[1])
    cbar = plt.colorbar(plt.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cbar_ax, orientation='horizontal')

    # Force all labels/ticks to white
    cbar.set_label(label2, color='white')
    cbar.ax.tick_params(colors='white')
    cbar.ax.xaxis.label.set_color('white')
    for label in cbar.ax.get_xticklabels():
        label.set_color('white')
    cbar.outline.set_edgecolor('white')

    # Save
    os.makedirs(outputloc, exist_ok=True)
    save_path = os.path.join(outputloc, f"{suffix}.png")
    plt.savefig(save_path, dpi=dpi, bbox_inches='tight', pad_inches=0, facecolor='black')
    plt.close(fig)
    
