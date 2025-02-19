from Dans_Diffraction import Crystal
from Dans_Diffraction.functions_crystallography import wave2energy, molecular_refractive_index
import numpy as np
import scipy as scp
from scipy import *
import scipy.optimize as optimize
import scipy.signal.windows as fft_windows

def mol_refractive_index(cif_file, lambda_xray, relative_density=1.0):
    xtl = Crystal(cif_file)   
    energy = wave2energy(lambda_xray)   # keV
    density = xtl.Properties.density() * relative_density # g/cm^3
    n, Delta, Beta = molecular_refractive_index(xtl.Properties.xtl.cif['_chemical_formula_structural'], energy, density)
    Theta_crit = np.rad2deg(np.sqrt(2*Delta))*2
    
    return Delta, Theta_crit, xtl.Properties.xtl.cif['_chemical_formula_structural']

# https://opara.zih.tu-dresden.de/xmlui/handle/123456789/1804 
def XRR(data, crit_ang, lambda_xray=1.51):
    """
    Returns equally spaced data for the Fourier analysis of a typical XRR dataset.
    Data must be an array with 2 Theta angle in the first column and the intensity in the second column
    """

    # rescale x-axis to get at the end nm from the FFT (introduce diffraction vector: s_cor)
    with np.errstate(invalid='ignore'):
        s_cor = 2*np.sqrt((np.cos(np.pi*crit_ang/2/180))**2 -
                          (np.cos(np.pi*data[:, 0]/2/180))**2)/(lambda_xray/10)

    # mask all entrys which are 'nan' due to values below the critical angle
    mask = np.logical_not(np.isnan(s_cor))
    s_cor = s_cor[mask]

    # subtract a background of 'theta^4 * intensity'
    intensity = s_cor**4*data[mask, 1]

    # create new x wave with even spacing going from the lowest
    # to the highest diffraction vector
    x = np.linspace(s_cor.min(), s_cor.max(), 1024)
    f = scp.interpolate.interp1d(s_cor, intensity, kind='cubic')
    return x, f(x)

def FFT(x, y, d=None, window=2, n=None):
    """
    Performs a (real) FFT using no window (0), a hanning (1), a hamming (2) or flattop (otherwise) window.
    d is data spacing; if no d is given, the spacing of the first two data values in x is taken
    x data have to be equaly spaced!
    n = number of values taken into account, if n<len(y): data are cropped; if n>len(y): data are zero-padded
    returns only one half of the full FFT (1st quadrant)
    """
    if d is None:
        d = (x[1]-x[0])
    N = len(y)
    # declaration of window functions
    if 0 == window:
        window = np.ones(N)
    elif 1 == window:
        window = fft_windows.hann(N)
    elif 2 == window:
        window = fft_windows.hamming(N)
    else:
        window = fft_windows.flattop(N)
    if n is None:
        n = N
    # calculate fft with all correction factors:
    # 2/N: up front to renormalize the FFT
    # mean(window): window correction factor
    yf = 2/N*np.abs(scp.fftpack.fft(window*y/np.mean(window), n=n))
    # Calculate the frequency axis
    xf = scp.fftpack.fftfreq(n, d=d)
    return xf[:n//2], yf[:n//2]

def fit_peaks(x, y, p0, min_peak_prominnence=0.15):
    """
    Fits a peak to the data using the function func
    """    
    def funcGauss(p, a, pmax, w):
        # Gaussian peak function
        return a*np.exp(-np.log(2)*((pmax-p)/(w/2))**2)
    
    def func(p, *pr_list):
        # Create the peak functions in a loop
        # so that the number of peaks can be varied
        pr_list = np.array(pr_list).flatten()
        r = np.zeros(len(p))
        for pr in range(np.size(pr_list)//3):
            a, pmax, w,= pr_list[pr*3:pr*3+3]
            r += funcGauss(p, a, pmax, w)
        return r
    
    # Mask the data
    msk = np.logical_and(x > 1e-3, x < 100)
    # Find the peakfinding 
    msk_pk = np.logical_and(x > 3, x < 100) 
    # Find the peaks
    dx = x[1]-x[0]
    pk_id, meta = scp.signal.find_peaks(y*msk_pk,prominence=np.max(y*msk_pk)*min_peak_prominnence, distance=5/dx)
    #filter out peaks on left edge
    bool_edge_pk = pk_id == np.argmax(msk_pk)
    pk_id = pk_id[~bool_edge_pk]
    # Initial guess for the peak parameters and bounds
    p0 = []
    lb = []
    ub = []
    for i in range(len(pk_id)):
        pk_x = x[pk_id[i]]
        pk_y = y[pk_id[i]]
        p0.append([pk_y, pk_x, 0.25])
        lb.append([pk_y*0.5, pk_x-0.1, 0.1])
        ub.append([pk_y*1.5, pk_x+0.1, 5])
    # Append a peak at x=0 as a noise peak
    p0.append([np.max(y[msk]), 0.0, 0.1])
    lb.append([0.0, -1.0, 0.1])
    ub.append([np.max(y[msk]), 1, 10])
    
    # Fit the peaks
    pk, pcov = optimize.curve_fit(func,x[msk],y[msk],p0=np.array(p0).flatten(), bounds=(np.array(lb).flatten(),np.array(ub).flatten()))
    bkg_peak = funcGauss(x, pk[-3], pk[-2], pk[-1])
    fit_model = func(x, pk)
    R2 = 1 - np.sum((y[msk]-fit_model[msk])**2)/np.sum((y[msk]-np.mean(y[msk]))**2)
    # std_err = np.sqrt(np.diag(pcov))
    
    good_fit = R2 > 0.95
    d = list()
    fit_curves = list()
    for i in range(len(pk_id)):
        if good_fit:
            msk_pk = np.abs(pk[3*i+1]-x) < pk[3*i+2]
        else:
            msk_pk = np.abs(x[pk_id[i]]-x) < 0.5
            
        peak_R2 = 1 - np.sum((y[msk_pk]-fit_model[msk_pk])**2)/np.sum((y[msk_pk]-np.mean(y[msk_pk]))**2)
        good_peak = peak_R2 > 0.9
        if good_fit and good_peak:
            d.append(pk[3*i+1])
            bkg = bkg_peak[np.abs(x-pk[3*i+1]).argmin()]
            pk_mod = funcGauss(x, pk[3*i], pk[3*i+1], pk[3*i+2]) + bkg
            fit_curves.append((x[msk_pk], pk_mod[msk_pk]))
        else:
            msk_pk = np.abs(x[pk_id[i]]-x) < 0.5
            d.append(x[pk_id[i]])
            fit_curves.append((x[msk_pk], y[msk_pk]))
    return d, fit_curves, fit_model, R2


def XRR_FFT_Analysis(data, meta, cif_file, min_peak_prominnence=0.15, relative_density=0.0):
    """
    Analyzes XRR data using a Fourier analysis
    """
    XRR_FT = dict()
    Meta = dict()

    if np.size(data) < 100:
        return XRR_FT, Meta
    
    twottheta = data[:, 0]
    intensity = data[:, 1]
    if "alpha_average" in meta:
        lam = meta["alpha_average"]
    elif "alpha1" in meta:
        lam = meta["alpha1"]
    elif "Wavelength" in meta:
        lam = meta["Wavelength"]

    Meta['delt_th'], Meta['crit_ang'], Meta['chem_formula'] = mol_refractive_index(cif_file, lam, relative_density)
    
    Meta['I_at_crit'] = np.interp(Meta['crit_ang'], twottheta, intensity)

    XRR_FT['x'], XRR_FT['y'] = FFT(*XRR(data, Meta['crit_ang'], lam), window=2, n=2**14)
    ds, Meta['peaks'], Meta['model'], Meta['R2'] = fit_peaks(XRR_FT['x'], XRR_FT['y'], 0, min_peak_prominnence)
    return XRR_FT, ds, Meta
