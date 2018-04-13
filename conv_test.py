from __future__ import print_function, division
import astropy.io.fits as pyfits
from scipy import signal


main_im = pyfits.open('HRI.weight.fits')[0].data
main_head = pyfits.open('HRI.weight.fits')[0].header
kernel = pyfits.open('ch1mykernel.fits')[0].data

conv_im = signal.convolve2d(main_im, kernel, mode='same')

print(conv_im.shape)

prihdr = main_head
hdu = pyfits.PrimaryHDU(conv_im, header=prihdr)
hdulist = pyfits.HDUList([hdu])
hdulist.writeto('LRI_Kband.weight.fits', clobber=True)
