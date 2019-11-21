# MKT J170456.2-482100: manuscript plots

This repo contains the Jupyter Notebooks and Python scripts for producing Figures 2, 4, 5, 6, 8,  10, 11 and 12 in

[![DOI](https://zenodo.org/badge/doi/10.1093/mnras/stz3027.svg)](https://doi.org/10.1093/mnras/stz3027)

(and is also available open access <a href="http://arxiv.org/abs/1911.07713">on ArXiv</a>).

## The data required to run these notebooks and scripts can be found here:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3548868.svg)](https://doi.org/10.5281/zenodo.3548868)
reStructedText

### The notebook ands scripts

<ul>
  <li><em>Process_Optical_Data.ipynb</em> (or <em>Process_Optical_Data.py</em>) uses the functions in <em>Optical_Processing_Functions.py</em> to read and manipulate the ASAS, ASAS-SN and KELT data</li>
  <li><em>LombScargle_Analysis.ipynb</em> (or <em>LombScargle_Analysis.py</em>) uses the outputs of <em>Process_Optical_Data.ipynb</em> and the functions in <em>Optical_LombScargle_Functions.py</em> to perform the (final from the manuscript analysis) Lomb-Scargle analysis of the optical observations</li>
  <li><em>Optical_Plotting.ipynb</em> contains the code for making Figures 4, 5, 8, and 11</li>
  <li><em>Radio_Plotting.ipynb</em> contains the code for making Figures 2 and 10</li>
  <li><em>SpectralEnergyDistribution.ipynb</em> contains the code for making Figure 6</li>
  <li><em>Gudel_Plot.ipynb</em> contains the code for making Figure 12</li>
</ul>

### The files needed for each script

<ul>
  <li><em>Process_Optical_Data.ipynb</em> (or <em>Process_Optical_Data.py</em>) needs:
    <ul>
      <li><a href="https://zenodo.org/record/3548868/files/ASAS_data.tsv?download=1">ASAS_data.tsv</a></li>
      <li><a href="https://zenodo.org/record/3548868/files/KELT_S36_lc_027056_V01_west_tfa.dat?download=1">KELT_S36_lc_027056_V01_west_tfa.dat</a></li>
      <li><a href="https://zenodo.org/record/3548868/files/KELT_S36_lc_027057_V01_east_tfa.dat?download=1">KELT_S36_lc_027057_V01_east_tfa.dat</a></li>
      <li><a href="https://zenodo.org/record/3548868/files/ASASSN.csv?download=1">ASASSN.csv</a></li>
    </ul>  
    </li>
  <li><em>LombScargle_Analysis.ipynb</em> (or <em>LombScargle_Analysis.py</em>) needs:
  <ul>
      <li><a href="https://zenodo.org/record/3548868/files/TYC_optical_semesters.npy?download=1">TYC_optical_semesters.npy</a></li>
  </ul>
  </li>
  <li><em>Optical_Plotting.ipynb</em> needs:
  <ul>
    <li><a href="https://zenodo.org/record/3548868/files/TYC_optical_semesters.npy?download=1">TYC_optical_semesters.npy</a></li>
    <li><a href="https://zenodo.org/record/3548868/files/TYC_optical_binned.npy?download=1">TYC_optical_binned.npy</a></li>
    <li><a href="https://zenodo.org/record/3548868/files/TYC_optical_binned_noOutliers.npy?download=1">TYC_optical_binned_noOutliers.npy</a></li>
    <li><a href="https://zenodo.org/record/3548868/files/TYC_LS_periods.npy?download=1">TYC_LS_periods.npy</a></li>
    <li><a href="https://zenodo.org/record/3548868/files/TYC_LS_periodErrors.npy?download=1">TYC_LS_periodErrors.npy</a></li>
    <li><a href="https://zenodo.org/record/3548868/files/SALT_radial_velocities.npy?download=1">SALT_radial_velocities.npy</a></li>
    <li><a href="https://zenodo.org/record/3548868/files/LCO_radial_velocities.npy?download=1">LCO_radial_velocities.npy</a></li>
  </ul>
  
  </li>
  <li><em>Radio_Plotting.ipynb</em> needs:</li>
  <ul>
    <li><a href="https://zenodo.org/record/3548868/files/TYC_optical_semesters.npy?download=1">TYC_optical_semesters.npy</a></li>
    <li><a href="https://zenodo.org/record/3548868/files/TYC_LS_periods.npy?download=1">TYC_LS_periods.npy</a></li>
    <li><a href="https://zenodo.org/record/3548868/files/TYC_LS_periodErrors.npy?download=1">TYC_LS_periodErrors.npy</a></li>
    <li><a href="https://zenodo.org/record/3548868/files/SALT_radial_velocities.npy?download=1">SALT_radial_velocities.npy</a></li>
    <li><a href="https://zenodo.org/record/3548868/files/LCO_radial_velocities.npy?download=1">LCO_radial_velocities.npy</a></li>
  <li><a href="https://zenodo.org/record/3548868/files/TYC_MeerKAT_fluxes.npy?download=1">TYC_MeerKAT_fluxes.npy</a></li>
  <li><a href="https://zenodo.org/record/3548868/files/TYC_MeerKAT_local_RMS.npy?download=1">TYC_MeerKAT_local_RMS.npy</a></li>
  </ul>
  <li><em>SpectralEnergyDistribution.ipynb</em> needs:</li>
  <ul>
  <li><a href="https://zenodo.org/record/3548868/files/TYC_SED_datapoints.txt?download=1">TYC_SED_datapoints.txt</a></li>
  <li><a href="https://zenodo.org/record/3548868/files/TYC_SED_model.txt?download=1">TYC_SED_model.txt</a></li>
  </ul>
  <li><em>Gudel_Plot.ipynb</em> needs:</li>
  <ul>
  <li><a href="https://zenodo.org/record/3548868/files/TYC_MeerKAT_ScaledFlux.npy?download=1">TYC_MeerKAT_ScaledFlux.npy</a></li>
  </ul>
</ul>

For full information on how to correctly cite the datasets listed here, please see the Zenodo:

<a href="https://doi.org/10.5281/zenodo.3548868"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.3548868.svg" alt="DOI"></a>

The DOI for this github repo is:

[![DOI](https://zenodo.org/badge/222704083.svg)](https://zenodo.org/badge/latestdoi/222704083)
