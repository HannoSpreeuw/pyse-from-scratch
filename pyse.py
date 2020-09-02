import numpy as np
from astropy.io import fits
# from scipy.stats import sigmaclip
from astropy.stats import sigma_clip
from scipy.ndimage import label, generate_binary_structure
from scipy.optimize import fsolve, root_scalar
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
kappa = 2.0

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

def find_true_std(sigma, clip_limit, clipped_std):
    help1=clip_limit/(sigma*np.sqrt(2))
    help2=np.sqrt(2*np.pi)*erf(help1)
    return sigma**2*(help2-2*np.sqrt(2)*help1*np.exp(-help1**2))-clipped_std**2*help2

true_sigma=fsolve(find_true_std, clipped_std, args=(clip_limit,clipped_std))[0]
check=find_true_std(true_sigma, clip_limit, clipped_std)
check_low=find_true_std(true_sigma*0.5, clip_limit, clipped_std)
check_high=find_true_std(true_sigma*1.5, clip_limit, clipped_std)


def var_helper(N, D, sigma_meas):
    """Correct for the fact the rms noise is computed from a clipped
    distribution.

    That noise will always be lower than the noise from the complete
    distribution.  The correction factor is a function of the computed
    rms noise only.
    """
    term1 = np.sqrt(2. * np.pi) * erf(N / np.sqrt(2.))
    term2 = 2. * N * np.exp(-N ** 2 / 2.)
    return (sigma_meas**2*term1 / (term1 - term2))-(D/N)**2

kappa = 1.98
test_clipped_std = 19.560791551187116
# test_clipped_std=0.7373492500752781
test_clip_limit = kappa*test_clipped_std
# test_clip_limit=1.460697620390798
test_sigma = fsolve(find_true_std, test_clipped_std, args=(test_clip_limit, test_clipped_std))[0]
test=find_true_std(test_sigma, test_clip_limit, test_clipped_std)
test_low=find_true_std(test_clipped_std, test_clip_limit, test_clipped_std)
test_high=find_true_std(test_clip_limit, test_clip_limit, test_clipped_std)


# test=var_helper(test_N, test_clip_limit, test_clipped_std)
# test_low=var_helper(10, test_clip_limit, test_clipped_std)
# test_upp=var_helper(0.1, test_clip_limit, test_clipped_std)

print()