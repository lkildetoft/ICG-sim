# ICG-sim
ICG-SIM: A Monte-Carlo simulation of light propagation and fluorescence from  indocyanine-green (ICG) in tissues. Can easily be modified for other fluorophores. Specifically developed to simulate mal-perfusion in the anastomosis after esophagectomy, as such allows for three distinct areas with ICG-concentrations growing at different speeds.

Based on MCML by Steven Jaques et. al. and its extension to simulate 
fluorescence, but written and adapted from scratch. 
Saves the current "face" and depth image of the absorbed and fluorescence 
photons at each iteration. These files may become very large in the end 
depending on the amount of photons/iterations/binning, so tread 
carefully! 
