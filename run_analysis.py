import numpy as np
import matplotlib.pyplot as plt
import os
import glob

from src.utils import *
from src.data_parser import parse_data
from src.XRR_FFT_Analysis import XRR_FFT_Analysis

# files to analyze (accepting .brml and .xy files) and cif file for the film-material
cif_file = os.path.join("example_data", "YIG.cif")
files = glob.glob(os.path.join("example_data", "*.brml")) + glob.glob(
    os.path.join("example_data", "*.xy")
)

# Parameters for the analysis
# min_peak_prominnence is the minimum peak prominnence for the peak detection
min_peak_prominnence = 0.16
# relative_density is a factor applied to the theoretical density of the film material,
# so a relative density of 1.0 means the film is fully dense, in accordance with the
# theoretical density of the material. This influences the critical angle and the scaling of the
# layer thicknesses axis in the FFT-plot 
relative_density = 1.0

for file in files:
    # Laod data and Meta data from the Bruker BRML file
    # Meta Data contains information like the XRAY-Tube used and the wavelength
    # of the XRAY radiation. Some files can have no data, so we skip them
    data_raw, data_meta = parse_data(file)
    if len(data_raw) == 0:
        continue

    # Iterate over all datasets in the BRML file (If the experiment was scripted
    # via jobs, there will be multiple datasets in the BRML file)
    for id, dataset in enumerate(data_raw):
        if id == 0:
            fig, ax = plt.subplots(len(data_raw), 2, squeeze=False)
        # Get the data and the meta data for the current dataset. Some datasets can
        # have few data points, when the experiment was cancelled for some reason,
        # so we skip them
        data = data_raw[dataset]
        if np.size(data) < 100:
            continue

        # Analyze the XRR data using a Fourier analysis
        # Returns the Fourier transformed data, the layer thicknesses (ds)
        # and the meta data, such as the fitted model and peak data
        XRR_FT, ds, Meta = XRR_FFT_Analysis(
            data, data_meta, cif_file, min_peak_prominnence, relative_density
        )
        Meta["dataset"] = dataset
        plot_XRR(data, Meta, XRR_FT, ds, ax[id])

    sample = os.path.basename(file).split(".")[0]
    fig.suptitle(sample)
    fig.set_size_inches(14, len(data_raw) * 5)
    fig.tight_layout(pad=0.5)
    fig.show()
    # plt.savefig(os.path.join('results',sample+'.svg'), bbox_inches = 'tight', pad_inches = 0)
    plt.close()
