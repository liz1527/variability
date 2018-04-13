from __future__ import print_function, division
import numpy as np
from scipy import signal
import sys
import astropy.io.fits as pyfits
import copy
import scipy.fftpack as fftpack

# use: python RL_deconv.py psf1 psf2 tolerance outfile
# psf1 is the narrower

# the purpose is to do psf matching. So data is one psf, and psf the other (the narrower one). The output is the convolution kernel required.

# grad = signal.convolve2d(a, b, boundary='fill', mode='same')

# fourier deconv, with noise suppression.
def deconv(star, psf):
    star_fft = fftpack.fftshift(fftpack.fftn(star))
    psf_fft = fftpack.fftshift(fftpack.fftn(psf))
    return fftpack.fftshift(fftpack.ifftn(fftpack.ifftshift(star_fft/(psf_fft+0.1))))


def RLfunc(prev, data, psf, flip_psf):
    conv1 =  data / signal.convolve2d(psf, prev, mode='same')
    conv2 = signal.convolve2d(flip_psf, conv1, mode='same')
    cur = prev * conv2
    return cur

def difference(prev, cur, data, psf):
    a = signal.convolve2d(psf, cur, mode='same')
    b = signal.convolve2d(psf, prev, mode='same')
    diff = np.sum(np.power(a-b,2.))#/(a.shape[0]*a.shape[1])
    #print(diff)
    return diff

def iter_step(prev, data, psf, flip_psf):
    updated_value = RLfunc(prev, data, psf, flip_psf) #<function from Wikipedia>
    return updated_value

def iterate(initial_guess, tolerance, data, psf):
    n_iter = 0
    cur = initial_guess
    flip_psf = np.flipud(np.fliplr(psf))
    while True:
        n_iter += 1
        prev, cur = cur, iter_step(cur, data, psf, flip_psf)
        if abs(difference(prev, cur, data, psf)) <= tolerance:
            break
    print(n_iter, tolerance)
    return cur


# psf is the smaller of the PSFs
psf = pyfits.open(sys.argv[1])[0].data

# data is the broader PSF
data = pyfits.open(sys.argv[2])[0].data

# starting guess, assume same as narrow PSF
init = copy.copy(psf)

# alternate start, fourier deconv.
#init = deconv(data.astype(float), psf.astype(float)).real

# run the deconv.
kernel = iterate(init, float(sys.argv[3]), data, psf)

# write out
hdu = pyfits.PrimaryHDU(kernel/sum(kernel))
hdu.writeto(sys.argv[4], clobber=True)

# write out test file
test = signal.convolve2d(kernel, psf, mode='same')
hdu = pyfits.PrimaryHDU(test/np.sum(test))
hdu.writeto(sys.argv[1]+'_check.fits', clobber=True)
