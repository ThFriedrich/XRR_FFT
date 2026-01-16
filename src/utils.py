import plotly.graph_objects as go
import numpy as np

def plot_XRR_plotly(data, Meta, XRR_FT, ds, fig, row):
    """
    Analyzes XRR data using a Fourier analysis - Plotly version
    """
    
    # Left plot - XRR data
    fig.add_trace(
        go.Scatter(x=data[:, 0], y=data[:, 1], mode='lines', name=Meta['dataset'], showlegend=True),
        row=row, col=1
    )
    
    fig.add_trace(
        go.Scatter(
            x=Meta['crit_ang'], y=Meta['I_at_crit'], 
            mode='markers', marker=dict(color='orange', size=8),
            name='Critical angle', showlegend=True
        ),
        row=row, col=1
    )
    
    # Add annotation for critical angle
    fig.add_annotation(
        x=Meta['crit_ang'][0] * 1.2, y=Meta['I_at_crit'][0],
        text=f"θ<sub>c</sub>({Meta['chem_formula']}) = {np.round(Meta['crit_ang'], 3)[0]}",
        showarrow=False,
        xanchor='left',
        bgcolor='rgba(255, 255, 255, 0.7)',
        bordercolor='gray',
        borderwidth=1,
        row=row, col=1
    )
    
    fig.update_xaxes(title_text="2θ", row=row, col=1)
    fig.update_yaxes(title_text="Intensity", type="log", row=row, col=1)
    
    # Right plot - FFT
    fig.add_trace(
        go.Scatter(x=XRR_FT['x'], y=XRR_FT['y'], mode='lines', name='Data', line=dict(width=1.5)),
        row=row, col=2
    )
    
    fig.add_trace(
        go.Scatter(x=XRR_FT['x'], y=Meta['model'], mode='lines', name='Model', line=dict(width=1.0)),
        row=row, col=2
    )
    
    sc_max = np.max(XRR_FT['y'])
    
    for ip in range(len(ds)):
        peak_name = f'Peak {ip + 1}'
        legendgroup = f'peak_{row}_{ip}'
        
        fig.add_trace(
            go.Scatter(
                x=Meta['peaks'][ip][0], y=Meta['peaks'][ip][1],
                mode='lines', name=peak_name,
                line=dict(dash='dot'),
                legendgroup=legendgroup,
                showlegend=True
            ),
            row=row, col=2
        )
        
        # Add vertical line as a trace so it can be toggled with the legend
        peak_max = np.max(Meta['peaks'][ip][1])
        fig.add_trace(
            go.Scatter(
                x=[ds[ip], ds[ip]], 
                y=[0, peak_max * 1.4],
                mode='lines',
                line=dict(dash='dash', color='gray'),
                legendgroup=legendgroup,
                showlegend=False,
                hoverinfo='skip'
            ),
            row=row, col=2
        )
        
        # Add invisible trace at annotation position to control annotation visibility
        fig.add_trace(
            go.Scatter(
                x=[ds[ip]], 
                y=[peak_max * 2],
                mode='markers+text',
                marker=dict(size=0.1, opacity=0),
                text=f"d<sub>{ip + 1}</sub> = {np.round(ds[ip], 2)} nm",
                textposition='middle right',
                textfont=dict(size=12, color='black'),
                legendgroup=legendgroup,
                showlegend=False,
                hoverinfo='skip'
            ),
            row=row, col=2
        )
        
        sc_max = np.max([sc_max, peak_max * 3])
    
    if len(ds) == 0:
        x_max = 100
    else:
        x_max = np.max(ds) * 3
    
    fig.update_xaxes(title_text="Thickness [nm]", range=[0, x_max], row=row, col=2)
    fig.update_yaxes(title_text="Intensity", range=[0, sc_max], row=row, col=2)
    
    # Add R² annotation
    closest_index = np.argmin(np.abs(XRR_FT['y'] - XRR_FT['y'][0] * 0.9))
    x_at_max = XRR_FT['x'][closest_index]
    
    fig.add_annotation(
        x=x_at_max + 1, y=np.max(XRR_FT['y']) * 0.9,
        text=f"R²<sub>Fit</sub> = {np.round(Meta['R2'], 3)}",
        showarrow=False,
        xanchor='left',
        bgcolor='rgba(255, 255, 255, 0.7)',
        bordercolor='gray',
        borderwidth=1,
        row=row, col=2
    )

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
    

