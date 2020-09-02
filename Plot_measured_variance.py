# Plot sigma_measured**2 as a function of sigma following formula 2.25 of my thesis for various values of D
import numpy as np
import pylab as pyl
from scipy.special import erf

sigma=np.arange(1,100)
true_sigma_measured=19.5608
clip_limit=1.75*true_sigma_measured

help1=clip_limit/(sigma*np.sqrt(2))
help2=np.sqrt(2*np.pi)*erf(help1)

variance_measured=sigma**2*(help2-2*np.sqrt(2)*help1*np.exp(-help1**2))/help2

pyl.plot(variance_measured)
pyl.plot(true_sigma_measured**2*np.ones(len(sigma)))
pyl.show()