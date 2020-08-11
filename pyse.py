import numpy as np
from astropy.io import fits
from scipy.stats import sigmaclip
from astropy.stats import sigma_clip
from scipy.ndimage import label, generate_binary_structure

input_fits = 'SOURCESINSERTED_10Jy.FITS'
hdulist=fits.open(input_fits)
dimensions = (hdulist[0].header['NAXIS1'], hdulist[0].header['NAXIS2'])
# print()
# print("The dimensions of this image are: ",dimensions)

scidata=hdulist['PRIMARY'].data
kappa = 3

# clipped_data_around_mean, low_mean, upp_mean = sigmaclip(scidata, kappa, kappa)
# mean_from_clipping_around_mean = clipped_data_around_mean.mean()
# median_from_clipping_around_mean = np.med

clipped_data_around_median, low_median, upp_median = sigma_clip(scidata, sigma=kappa, masked=False, return_bounds=True)
mean_from_clipping_around_median = clipped_data_around_median.mean()
median_from_clipping_around_median = np.median(clipped_data_around_median)
TKP_mode_estimator = 2.5*median_from_clipping_around_median-1.5*mean_from_clipping_around_median

sci_clip = scidata >= 7.5782
structure=generate_binary_structure(2,2)
# sci_labels, sci_num = label(sci_clip)
sci_labels, sci_num = label(sci_clip, structure=structure)
print("Number of islands = {0}".format(sci_num))
