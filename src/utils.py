import matplotlib.pyplot as plt
import numpy as np

ax_color = 'gray'
plt.rcParams.update({
    "figure.facecolor":  (1.0, 1.0, 1.0, 0.0),  
    "axes.facecolor":    (1.0, 1.0, 1.0, 0.0),  
    "savefig.facecolor": (1.0, 1.0, 1.0, 0.0),
    "axes.edgecolor":    ax_color, 
    "xtick.color":       ax_color, 
    "ytick.color":       ax_color,
    "text.color":        ax_color,
    "axes.labelcolor":   ax_color,
    "axes.titlecolor":   ax_color,
    "legend.facecolor":  (1.0, 1.0, 1.0, 0.8),
    "axes.grid":         True,
})

def plot_XRR(data, Meta, XRR_FT, ds, ax):
    """
    Analyzes XRR data using a Fourier analysis
    """
        
    ax[0].plot(data[:, 0], data[:, 1])
    ax[0].scatter(
        Meta['crit_ang'], Meta['I_at_crit'], color="tab:orange", label="Critical angle"
    )
    ax[0].annotate(
        r"$\theta_c({}) = $".format(Meta['chem_formula']) + str(np.round(Meta['crit_ang'], 3)[0]),
        (Meta['crit_ang'] * 1.2, Meta['I_at_crit']),
        ha="left",
        fontsize=12,
        color="k",
        bbox=dict(
            facecolor=(1.0, 1.0, 1.0, 0.7),
            edgecolor="gray",
            boxstyle="round,pad=0.5",
        ),
    )
    ax[0].set_xlabel(r"$2\theta$")
    ax[0].set_ylabel("Intensity")
    ax[0].set_yscale("log")
    ax[0].legend([Meta['dataset'], "Critical angle"], labelcolor="linecolor")
    ax[0].set_title("XRR data")
    ax[1].plot(XRR_FT['x'], XRR_FT['y'], label="Data", linewidth=1.5)
    ax[1].plot(XRR_FT['x'], Meta['model'], label="Model", linewidth=1.0)
    sc_max = np.max(XRR_FT['y'])
    for ip in range(len(ds)):
        ax[1].plot(
            Meta['peaks'][ip][0],
            Meta['peaks'][ip][1],
            label="Fit Peak " + str(ip + 1),
            linestyle=":",
        )
        ax[1].vlines(ds[ip], 0, np.max(Meta['peaks'][ip][1]) * 1.4, linestyles="dashed")
        ax[1].annotate(
            r"$d_{} = $".format(ip + 1) + str(np.round(ds[ip], 2)) + " nm",
            (ds[ip], np.max(Meta['peaks'][ip][1]) * 2),
            ha="left",
            va="bottom",
            fontsize=12,
            color="k",
            bbox=dict(
                facecolor=(1.0, 1.0, 1.0, 0.7),
                edgecolor="gray",
                boxstyle="round,pad=0.5",
            ),
        )
        sc_max = np.max([sc_max, np.max(Meta['peaks'][ip][1]) * 3])
    if len(ds) == 0:
        ax[1].set_xlim(0, 100)
    else:
        ax[1].set_xlim(0, np.max(ds) * 3)
    ax[1].set_xlabel("Thickness [nm]")
    ax[1].set_ylabel("Intensity")
    ax[1].set_ylim(0, sc_max)
    closest_index = np.argmin(np.abs(XRR_FT['y'] - XRR_FT['y'][0] * 0.9))
    x_at_max = XRR_FT['x'][closest_index]
    ax[1].annotate(
        r"$R^2_{\text{Fit}} = $" + str(np.round(Meta['R2'], 3)),
        (x_at_max + 1, np.max(XRR_FT['y']) * 0.9),
        ha="left",
        fontsize=12,
        color="k",
        bbox=dict(
            facecolor=(1.0, 1.0, 1.0, 0.7),
            edgecolor="gray",
            boxstyle="round,pad=0.5",
        ),
    )
    ax[1].legend(labelcolor="linecolor")
    ax[1].set_title("FFT of XRR data")
    

