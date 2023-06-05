# schneida_tools

This package contains software used by Schneider et al. to analyze ice sheet 
climate. Within its core utilities are Python codes that use external libraries
(e.g., numpy, netCDF4, cartopy, pyproj, GDAL) to draw elevation
contours onto a Lambert azimuthal equal-area map projection using Greenland
and Antarctic ice sheet digital elevation model data from Howat et al. (2014)
and Bamber et al. (2009), respectively. Also included are analysis scripts
to inter-compare reanalysis precipitation over Greenland and Antarctica to
net accumulation rates compiled in the the Surface Mass Balance and Snow on Sea 
Ice Working Group (SUMup) dataset (Koenig & Montgomery, 2019).

## Mapping tutorial

### 0. Download and install the package and environment
	bash$ git clone https://github.com/amschne/lnd_srf_for.git; cd lnd_srf_for
	lnd_srf_for/$ conda env create -n [ENVNAME] --file environment.yml # this usually takes a few minutes
	lnd_srf_for/$ conda activate [ENVNAME]
	lnd_srf_for/$ pip install .

### 1. Test the installation - this should generate a short series of map demos
	 lnd_srf_for/$ pytest tests/

## Running the analyses of Schneider et al. (2023)

### 2. There are three main python scripts that reproduce figures in the manuscript
	lnd_srf_for/$ python schneida_tools/verify_precip_era5.py
	lnd_srf_for/$ python schneida_tools/analysis_era5.py
	lnd_srf_for/$ python schneida_tools/taylor.py

These scripts write the figure files  "p_cruncep_wfde5_sumup_gris_ais_precip.png," sumup_accum_locs.png" and "taylor_e5_we_g3_cn_m2_sumup.pdf" to the "results" directory and "pyplot_figure.png."

## Data access
All input data needed to run the analyses have been processed and made 
available within this repository. However, raw data should be obtained from the
original sources. If using these data as part of new analyses
using this repository, please attribute proper credit to the original
authors by citing the references as indicated below.

Regional climate model data used for estimating the net surface vapor fluxes 
are available at ftp://ftp.climato.be/fettweis/MARv3.5/Greenland/ 
(Fettweis et al., 2017), for Greenland, and at 
https://doi.org/10.5281/zenodo.4459259 (Kittel et al., 2021), for Antarctica. 
ERA5 and WFDE5 data are available at the Climate Data Store via the Copernicus 
program (DOI:10.24381/cds.f17050d7; DOI:10.24381/cds.20d54e34)
(Hersbach et al., 2020; Cucchi et al., 2020). CRUNCEP data are available via 
the National Center for Atmospheric Research Research Data Archive (DOI: 
10.5065/PZ8FF017) under the Creative Commons Attribution 4.0 International 
License (Viovy, 2018). The GSWP3 dataset (DOI: 10.20783/DIAS.501) is available 
from the Data Integration and Analysis System Program via the Japan Ministry of 
Education, Culture, Sports, Science and Technology (Kim, 2017). MERRA-2 monthly 
precipitation data are available via the U.S. National Aeronautics and Space 
Administration Global Modeling and Assimilation Office (DOI: 
10.5067/0JRLVL8YV2Y4) (Global Modeling And Assimilation Office, 2015).

Please also note that this repository contains data published by Bamber et al.
(2009), Howat et al. (2014), and Koenig & Montgomery (2019), the latter
which includes numerous datasets from ice core measurement therein. Please 
refer to Koenig & Montgomery (2019) or Schneider et al. (2023) for more details.

## References

Bamber, J., Gomez-Dans, J. & Griggs, J. (2009). _Antarctic 1 km Digital Elevation Model (DEM) from Combined ERS-1 Radar and ICESat Laser Satellite Altimetry, Version 1._ NASA National Snow and Ice Data Center DAAC. Retrieved 2020-11-17, from http://nsidc.org/data/NSIDC-0422/versions/1 (type: dataset) doi: 10.5067/H0FQ1KL9NEKM

Cucchi, M., Weedon, G. P., Amici, A., Bellouin, N., Lange, S., Mu¨ller Schmied, H., . . . Buontempo, C.	(2020).	WFDE5: bias-adjusted ERA5 reanalysis data for impact studies.	_Earth System Science Data, 12_(3), 2097–2120.	Retrieved from https://essd.copernicus.org/articles/12/2097/2020/ doi:10.5194/essd-12-2097-2020

Fettweis, X., Box, J. E., Agosta, C., Amory, C., Kittel, C., Lang, C., . . . Gall´ee, H. (2017, April).	Reconstructions of the 1900–2015 Greenland ice sheet surface mass balance using the regional climate MAR model. The _Cryosphere, 11_(2), 1015–1033.	Retrieved 2020-05-01, from https://www.the-cryosphere.net/11/1015/2017/ doi: 10.5194/tc-11-1015-2017

Global Modeling And Assimilation Office. (2015). _MERRA-2 tavgM 2d flx nx:2d,Monthly mean,Time-Averaged,Single-Level,Assimilation,Surface Flux Diagnostics V5.12.4_. NASA Goddard Earth Sciences Data and Information Services Center.	Retrieved 2023-03-30, from https://disc.gsfc.nasa.gov/datacollection/M2TMNXFLX 5.12.4.html (Type: dataset) doi: 10.5067/0JRLVL8YV2Y4

Hersbach, H., Bell, B., Berrisford, P., Hirahara, S., Hor´anyi, A., Mun˜oz-Sabater, J., . . . Th´epaut, J.-N. (2020). The ERA5 global reanalysis.	_Quarterly Journal of the Royal Meteorological Society, 146_ (730), 1999–2049.	Retrieved from https://rmets.onlinelibrary.wiley.com/doi/abs/10.1002/qj.3803 (eprint: https://rmets.onlinelibrary.wiley.com/doi/pdf/10.1002/qj.3803)	doi: https://doi.org/10.1002/qj.3803

Howat, I. M., Negrete, Al., & Smith, B. E. (2014, August). The Greenland Ice Mapping Project (GIMP) land classification and surface elevation data sets. _The Cryosphere, 8_(4), 1509-1518. Retrieved 2020-09-23, from https://tc.copernicus.org/articles/8/1509/2014/ doi: 105194/tc-8-1509-2014

Kim, H. (2017). _Global Soil Wetness Project Phase 3 Atmospheric Boundary Conditions (Experiment 1)_. Retrieved from https://doi.org/10.20783/DIAS.501 (Type: dataset) doi: 10.20783/DIAS.501

Kittel, C., Amory, C., Agosta, C., Jourdain, N. C., Hofer, S., Delhasse, A., . . . Fettweis, X. (2021). Diverging future surface mass balance between the Antarctic ice shelves and grounded ice sheet. _The Cryosphere, 15_(3), 1215–1236.	Retrieved from https://tc.copernicus.org/articles/15/1215/2021/ doi: 10.5194/tc-15-1215-2021

Lora Koenig, & Lynn Montgomery. (2019). _Surface Mass Balance and Snow Depth on Sea Ice Working Group (SUMup) accumulation on land ice subdataset, Greenland and Antarctica, 1987-2018._ Arctic Data Center. urn:uuid:ab62eb85-0417-48e6-89d2-9fee20d1d1a0.

Viovy, N. (2018). _CRUNCEP Version 7 - _Atmospheric Forcing Data for the Community Land Model_.	Research Data Archive at the National Center for Atmospheric Research, Computational and Information Systems Laboratory. Retrieved from http://rda.ucar.edu/datasets/ds314.3/