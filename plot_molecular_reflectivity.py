import numpy as np
from Dans_Diffraction.functions_crystallography import molecular_reflectivity
import plotly.graph_objects as go

compounds = {
    "YIG": {
        "formula": "Y3Fe5O12",
        "density": 5.17  # g/cm^3
    },
    "Si": {
        "formula": "Si",
        "density": 2.33  # g/cm^3
    },
    "Al2O3": {
        "formula": "Al2O3",
        "density": 3.98  # g/cm^3
    },
    "Be": {
        "formula": "Be",
        "density": 1.85  # g/cm^3
    }
}
energy_kev = 8.048  # Cu K-alpha
angles = np.linspace(0, 5, 1000)  # incident angles in degrees

fig = go.Figure()
for name, props in compounds.items():
    chemical_formula = props["formula"]
    density = props["density"]
    R = molecular_reflectivity(chemical_formula, energy_kev, density, angles)

    fig.add_trace(go.Scatter(x=angles, y=R, mode='lines', name=name))

fig.update_layout(
    title=f'Reflectivity at {energy_kev} keV',
    xaxis_title='Incident Angle (degrees)',
    yaxis_title='Reflectivity',
    yaxis_type='log'
)
fig.show()

# fig.write_image("refractive_index_plot.svg")