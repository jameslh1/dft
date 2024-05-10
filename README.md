# dft

Routines to take Discrete Fourier transform, including correct amplitude and absolute phase. Designed to work in THz time-domain spectroscopy.

Code: 
- maths.py provides a class to take the DFT including the correct amplitude and absolute phase, "takeFFT". 
- eos.py provides routines that perform a frequency-domain response function correction, required to convert an experimental electro-optic signal to the THz electric field.
- plot.py makes some standard plots, such as the time-domain, frequency-domain amplitude and phase spectra.
