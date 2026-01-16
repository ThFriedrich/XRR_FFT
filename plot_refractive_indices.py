import numpy as np
from Dans_Diffraction.functions_crystallography import molecular_refractive_index, energy2wave
from Dans_Diffraction import Crystal, Structures
import plotly.graph_objects as go


# Build a crystal from a packaged structure
s = Structures()
copper = s.Copper.build()
cobalt = s.Cobalt.build()
beryllium = s.Beryllium.build()
# Nickel: FCC structure, a = 3.524 Å, space group Fm-3m (225)
nickel = Crystal()
nickel.name = 'Nickel'
nickel.Cell.latt([3.524, 3.524, 3.524, 90, 90, 90])
nickel.new_atoms(u=[0.0], v=[0.0], w=[0.0], type='Ni')
nickel.Symmetry.load_spacegroup(225)
nickel.generate_structure()

edge_cu = copper.Properties.xray_edges()
# Emission lines are energy differences between shells:
# K_alpha1: L3 -> K transition (index 3 -> 0)
# K_alpha2: L2 -> K transition (index 2 -> 0)
# K_beta1:  M  -> K transition (index 4 -> 0, approximate with M1)
K_alpha1 = edge_cu[1][0] - edge_cu[1][3]  # Cu K-alpha1: K - L3 = 8.046 keV
K_alpha2 = edge_cu[1][0] - edge_cu[1][2]  # Cu K-alpha2: K - L2 = 8.027 keV
K_beta1 = edge_cu[1][0] - edge_cu[1][4]   # Cu K-beta1:  K - M1 ≈ 8.857 keV

energy_kev = np.linspace(5, 15, 1024)  # energy range in keV
abs_structs = {
    "Cu": copper,
    "Co": cobalt,
    "Be": beryllium,
    "Ni": nickel    
}

fig = go.Figure()
for name, props in abs_structs.items():
    density = props.Properties.density()  # g/cm^3
    chemical_formula = name
    n, Delta, Beta  = molecular_refractive_index(chemical_formula, energy_kev, density)
    lambda_xray = energy2wave(energy_kev)  # in Angstrom
    lin_abs_coeff = 4 * np.pi * Beta / lambda_xray * 10000 # in 1/um

    fig.add_trace(go.Scatter(x=energy_kev, y=lin_abs_coeff, mode='lines', name=name))
    # fig.add_trace(go.Scatter(x=energy_kev, y=Beta, mode='lines', name=name + ' Beta', line=dict(dash='dash')))


fig.add_vline(x=K_alpha1, line=dict(color='Red', dash='dash'), annotation_text='Cu K-alpha1', annotation_position='top right')
fig.add_vline(x=K_alpha2, line=dict(color='Orange', dash='dash'), annotation_text='Cu K-alpha2', annotation_position='top right')
fig.add_vline(x=K_beta1, line=dict(color='Green', dash='dash'), annotation_text='Cu K-beta1', annotation_position='top right')
fig.update_layout(
    title=f'Linear Absorption Coefficient',
    xaxis_title='Energy (keV)',
    yaxis_title='Linear Absorption Coefficient (1/µm)',
    yaxis_type='linear'
)
fig.show()

# fig.write_image("refractive_index_plot.svg")