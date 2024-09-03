# dft

Routines to take Discrete Fourier transform, including correct amplitude and absolute phase. Designed to work in THz time-domain spectroscopy.

Code: 
- maths.py provides a class to take the DFT including the correct amplitude and absolute phase, <b>takeFFT</b>, as well as a <b>timeShift</b> function that alters the linear phase in the frequency domain to shift a time-domain pulse via the Fourier time-shift theorem.
- eos.py provides routines that perform a frequency-domain response function correction, required to convert an experimental electro-optic signal to the THz electric field.
- plot.py makes some standard plots, such as the time-domain, frequency-domain amplitude and phase spectra.

Installation:
- on Linux you can download the source, then run
     $ sudo pip install .          - to use the default python interpreter
  or $ sudo pip3 install .         - to explicitly use python3

- Verify the installation by running "import dft" from inside a python environment
