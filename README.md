# MKT J170456.2-482100: manuscript plots

This repo contains the Jupyter Notebooks and Python scripts for producing Figures 2, 4, 5, 6, 8,  10, 11 and 12 in
[![DOI](https://zenodo.org/badge/doi/10.1093/mnras/stz3027.svg)](https://doi.org/10.1093/mnras/stz3027)
(and is also available open access <a href="http://arxiv.org/abs/1911.07713">on ArXiv</a>).

<ul>
  <li><em>Process_Optical_Data.ipynb</em> (or <em>Process_Optical_Data.py</em>) uses the functions in <em>Optical_Processing_Functions.py</em> to read and manipulate the ASAS, ASAS-SN and KELT data</li>
  <li><em>LombScargle_Analysis.ipynb</em> (or <em>LombScargle_Analysis.py</em>) uses the outputs of <em>Process_Optical_Data.ipynb</em> and the functions in <em>Optical_LombScargle_Functions.py</em> to perform the (final from the manuscript analysis) Lomb-Scargle analysis of the optical observations</li>
  <li><em>Optical_Plotting.ipynb</em> contains the code for making Figures 4, 5, 8, and 11</li>
  <li><em>Radio_Plotting.ipynb</em> contains the code for making Figures 2 and 10</li>
  <li><em>SpectralEnergyDistribution.ipynb</em> contains the code for making Figure 6</li>
  <li><em>Gudel_Plot.ipynb</em> contains the code for making Figure 12</li>
</ul>
