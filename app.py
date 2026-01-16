import streamlit as st
import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import os
from src.utils import *
from src.data_parser import parse_data_bytesIO
from src.XRR_FFT_Analysis import XRR_FFT_Analysis
import io

def main():
    st.set_page_config(
        page_title="XRR FFT Analysis Application",
        page_icon="imgs/Logo_RGB_farbig.png",
        layout="wide",
        initial_sidebar_state="auto",
        menu_items=None,
    )
    st.title("XRR FFT Analysis Application")

    # File uploader for .brml and .xy files
    col1, col2 = st.columns(2)
    with col1:
        st.write("Please upload .brml or .xy files for analysis.")
        uploaded_file = st.file_uploader("Choose .brml or .xy files", type=["brml", "xy"], accept_multiple_files=False)
    with col2:
        st.write("Please upload the CIF file for the film material.")
        cif_file = st.file_uploader("Choose a CIF file", type=["cif"], accept_multiple_files=False)
        if cif_file is not None:
            cif_path = os.path.join("/tmp", cif_file.name)
            with open(cif_path, "wb") as f:
                f.write(cif_file.getbuffer())
        
    col3, col4 = st.columns(2)
    # Sliders for parameters
    with col3:
        min_peak_prominence = st.slider("Minimum Peak Prominence", 0.0, 1.0, 0.16)
    with col4:
        relative_density = st.slider("Relative Density", 0.0, 2.0, 1.0)

    if uploaded_file and cif_file:
        # Read the content of the uploaded file
        file_content = uploaded_file.read()
        file_type = uploaded_file.name.split(".")[-1]
        data_raw, data_meta = parse_data_bytesIO(io.BytesIO(file_content),file_type)
        if len(data_raw) == 0:
            st.warning(f"No data found in {uploaded_file.name}. Skipping.")
        else:
            sample = os.path.basename(uploaded_file.name).split(".")[0]
            
            for id, dataset in enumerate(data_raw):
                data = data_raw[dataset]
                if np.size(data) < 100:
                    st.warning(f"Dataset {id} in {uploaded_file.name} has insufficient data points. Skipping.")
                    continue

                XRR_FT, ds, Meta = XRR_FFT_Analysis(data, data_meta, cif_path, min_peak_prominence, relative_density)
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

                # Display the plot in Streamlit
                st.plotly_chart(fig, use_container_width=True, key=f"xrr_plot_{id}")
    else:
        st.info("Please upload a CIF file and a data file to start the analysis.")

if __name__ == "__main__":
    main()