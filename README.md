# Extragalactic Hertzsprung gap stars and LRN progenitor candidates to predict and identify transients

### Catalog of HST yellow supergiants in nearby galaxies (Tranin et al. submitted to A&A, arxiv https://arxiv.org/abs/2409.11347)

Luminous red novae (LRNe) are intermediate-luminosity optical transients observed in the Milky Way and nearby galaxies, uniquely probing the evolution of binary stars, as they represent the ejection of their shared envelope of gas, leading to the shrinking of their orbit, or their merger. 

Here we use HST imaging of nearby galaxies to find their possible progenitors and precursors, making it possible to predict LRNe before their outburst and quickly identify the new transients matching their position.

Content of this folder:
* YSG_candidates_v0.1.zip: contains the fits file of the catalog of YSG candidates
* galaxy_summary_new.fits: catalog of galaxies considered in this study
* casjobs_query.py, retrieve_mast_data.py: codes to retrieve HST data from HSCv3 and/or MAST
* transform_mast_cat.py: code to transform MAST catalogs into HSCv3-like catalogs
* plot_cmd.py: code to analyze the CMD of a galaxy, select candidates, and plot the result
* stack_candidates.py: code to stack the candidates of all galaxies and flag contaminants
* study_gal.py: code to download and show HST science color images of a given galaxy

If you use this catalog please cite (Tranin et al. submitted to A&A, https://arxiv.org/abs/2409.11347)
Feel free to reach me at htranin[at]icc.ub.edu

Good luck with your transients ✨

**NB: material presented in this folder is preliminary, the final version will be uploaded after article acceptance**
