import numpy as np
from astropy.io import fits
# from scipy.stats import sigmaclip
from astropy.stats import sigma_clip
from scipy.ndimage import label, generate_binary_structure
from scipy.optimize import fsolve
from scipy.special import erf

input_fits = 'SOURCESINSERTED_100Jy.FITS'
# input_fits = 'UNCORRELATED_GAUSSIAN_NOISE_WITH_STD_1Jy.FITS'
hdulist=fits.open(input_fits)
dimensions = (hdulist[0].header['NAXIS1'], hdulist[0].header['NAXIS2'])
# print()
# print("The dimensions of this image are: ",dimensions)

scidata=hdulist['PRIMARY'].data
maximum = scidata.max()
minimum = scidata.min()
kappa = 1.5

unclipped_std = np.std(scidata)
unclipped_mean = np.mean(scidata)
# clipped_data_around_mean, low_mean, upp_mean = sigmaclip(scidata, kappa, kappa)
# mean_from_clipping_around_mean = clipped_data_around_mean.mean()
# median_from_clipping_around_mean = np.med

clipped_data_around_median, low_median, upp_median = sigma_clip(scidata, sigma=kappa, masked=False, return_bounds=True)
mean_from_clipping_around_median = clipped_data_around_median.mean()
median_from_clipping_around_median = np.median(clipped_data_around_median)
TKP_mode_estimator = 2.5*median_from_clipping_around_median-1.5*mean_from_clipping_around_median

clipped_std = np.std(clipped_data_around_median)

# sci_clip = scidata >= 7.5782
# structure=generate_binary_structure(2,2)
# # sci_labels, sci_num = label(sci_clip)
# sci_labels, sci_num = label(sci_clip, structure=structure)
# print("Number of islands = {0}".format(sci_num))

clip_limit=upp_median-median_from_clipping_around_median

def find_true_variance(sigma):
    help1=clip_limit/(sigma*np.sqrt(2))
    help2=np.sqrt(2*np.pi)*erf(help1)
    return sigma**2*(help2-2*np.sqrt(2)*help1*np.exp(-help1**2))-clipped_std**2*help2

true_variance=fsolve(find_true_variance, clipped_std)[0]
true_sigma=np.sqrt(true_variance)
print()