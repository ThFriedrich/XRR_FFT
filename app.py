import streamlit as st
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os
from src.utils import *
from src.data_parser import parse_data_bytesIO
from src.XRR_FFT_Analysis import XRR_FFT_Analysis, density_from_cif, ANODE_WAVELENGTHS
import io

def main():
    st.set_page_config(
        page_title="XRR FFT Analysis Application",
        page_icon="imgs/Logo_RGB_farbig.png",
        layout="wide",
        initial_sidebar_state="auto",
        menu_items=None,
    )
    col01, col02 = st.columns(2)
    with col01:
        st.title("XRR FFT Analysis Application")
    with col02:
        st.image("imgs/Logo_lang_RGB_farbig.svg", width=300)

    with st.expander("ℹ️  How to use this tool", expanded=False):
        st.markdown(
            "1. Upload an XRR measurement as a `.brml` or `.xy` file.\n"
            "2. Specify the film material so the critical angle can be estimated:\n"
            "   - Enter a **chemical formula** (e.g. `Si`, `Fe2O3`, `Y3Fe5O12`) and a **density** in g/cm³, *or*\n"
            "   - Upload a **CIF file** (optional) to auto-fill formula and density from the crystal structure.\n"
            "3. Tune the **Minimum Peak Prominence** slider to control which FFT peaks are picked as layer thicknesses.\n\n"
            "Notes: The density value strongly affects the critical angle and therefore the extracted thicknesses — "
            "use a realistic film density (thin films are often less dense than the bulk crystal). "
            "Values from the CIF can be overridden manually at any time."
        )

    # File uploader for .brml and .xy files
    col1, col2 = st.columns(2)
    with col1:
        st.write("Please upload .brml or .xy files for analysis.")
        uploaded_file = st.file_uploader("Choose .brml or .xy files", type=["brml", "xy"], accept_multiple_files=False)
    with col2:
        st.write("Optionally upload a CIF file to derive formula and density.")
        cif_file = st.file_uploader("Choose a CIF file (optional)", type=["cif"], accept_multiple_files=False)

    cif_formula = None
    cif_density = None
    if cif_file is not None:
        cif_path = os.path.join("/tmp", cif_file.name)
        with open(cif_path, "wb") as f:
            f.write(cif_file.getbuffer())
        try:
            cif_formula, cif_density = density_from_cif(cif_path)
            st.success(f"From CIF: formula = {cif_formula}, density = {cif_density:.3f} g/cm³")
        except Exception as e:
            st.warning(f"Could not parse CIF: {e}")

    # Parse the data file early so we can pre-fill the anode/wavelength controls
    data_raw, data_meta = {}, {}
    file_lam = None
    detected_anode = None
    if uploaded_file is not None:
        file_content = uploaded_file.read()
        file_type = uploaded_file.name.split(".")[-1]
        data_raw, data_meta = parse_data_bytesIO(io.BytesIO(file_content), file_type)
        for key in ("alpha_average", "alpha1", "Wavelength"):
            if key in data_meta:
                file_lam = float(data_meta[key])
                break
        target = data_meta.get("target_type")
        if target and target in ANODE_WAVELENGTHS:
            detected_anode = target

    anode_options = list(ANODE_WAVELENGTHS.keys()) + ["Custom"]
    default_anode = detected_anode or "Cu"

    col_f, col_a, col_pp, col_d = st.columns([2, 1, 1, 1])
    with col_f:
        chem_formula = st.text_input(
            "Chemical formula (e.g. Si, Fe2O3, Y3Fe5O12)",
            value=cif_formula if cif_formula else "Si",
            key=f"formula_{cif_formula}",
        )
    with col_a:
        anode = st.selectbox(
            "X-ray anode / source",
            anode_options,
            index=anode_options.index(default_anode),
            help="Auto-selected from .brml metadata when available. "
                 "Choose 'Custom' for synchrotron or other sources.",
            key=f"anode_{detected_anode}",
        )
    with col_pp:
        min_peak_prominence = st.slider("Minimum Peak Prominence", 0.0, 1.0, 0.16)
    with col_d:
        density = st.slider(
            "Density (g/cm³)",
            0.1, 25.0,
            float(cif_density) if cif_density else 2.33,
            step=0.01,
            key=f"density_{cif_density}",
        )

    # Resolve the wavelength used for analysis
    if anode == "Custom":
        custom_default = float(file_lam) if file_lam else 1.0
        lam = st.number_input(
            "Custom X-ray wavelength λ (Å)",
            min_value=0.05, max_value=10.0,
            value=custom_default, step=0.0001, format="%.4f",
            help="Use this for synchrotron beamlines or non-standard sources.",
        )
    else:
        lam = float(file_lam) if file_lam else ANODE_WAVELENGTHS[anode]

    src_note = (
        f"from file ({uploaded_file.name})" if file_lam and anode != "Custom"
        else "custom value" if anode == "Custom"
        else f"Kα1 of {anode}"
    )
    st.caption(f"Using wavelength λ = {lam:.4f} Å ({src_note}).")

    if uploaded_file:
        if len(data_raw) == 0:
            st.warning(f"No data found in {uploaded_file.name}. Skipping.")
        else:
            sample = os.path.basename(uploaded_file.name).split(".")[0]
            
            for id, dataset in enumerate(data_raw):
                data = data_raw[dataset]
                if np.size(data) < 100:
                    st.warning(f"Dataset {id} in {uploaded_file.name} has insufficient data points. Skipping.")
                    continue

                XRR_FT, ds, Meta = XRR_FFT_Analysis(
                    data, data_meta, chem_formula, min_peak_prominence, density, anode, lam=lam
                )
                Meta["dataset"] = dataset
                
                # Create a separate figure for each dataset
                fig = make_subplots(
                    rows=1, cols=2,
                    subplot_titles=["XRR data", "FFT of XRR data"],
                    horizontal_spacing=0.12
                )
                
                plot_XRR_plotly(data, Meta, XRR_FT, ds, fig, 1)
                
                fig.update_layout(
                    title_text=f"{sample} - {dataset}",
                    height=500,
                    showlegend=True
                )

                # Configure plotly modbar with SVG download option
                config = {
                    'toImageButtonOptions': {
                        'format': 'svg',
                        'filename': f'xrr_fft_analysis_{sample}_{dataset}',
                        'scale': 2

                    }
                }

                # Display the plot in Streamlit
                st.plotly_chart(fig, use_container_width=True, key=f"xrr_plot_{id}", config=config)
    else:
        st.info("Please upload a data file to start the analysis.")

if __name__ == "__main__":
    main()