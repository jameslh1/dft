from numpy import pi, linspace, exp, gradient, loadtxt, max, sqrt

import matplotlib.pyplot as plt
plt.style.use('nature')

##### import the DFT class library
from dft import *					# provides "takeFFT" class


def getdata(filename, t_offset=0.0e-12, normalise=False, expansion=1):
	d=loadtxt(filename,comments='%')

	if normalise:
		ETHz_raw=d[:,1]/max((d[:,1]))
	else:
		ETHz_raw=d[:,1]
	x_raw=d[:,0]

	x,ETHz = zeropad(x_raw,ETHz_raw,expansion=expansion)

	t=(x/0.15-0.0) * 1e-12	# mm to s including offset
	t=t-t_offset
	
	totalInt=integrateTHzPower(t,ETHz)		# 
	normTHzPower=ETHz**2/totalInt

	integrand = t * normTHzPower
	mu = trapz(integrand, x=t)
#	print(mu)

	integrand = (t-mu)**2 * normTHzPower
	variance = trapz(integrand,x=t)

	print(filename)
	print('-- pulse duration (std. dev.) = %4.3ffs' % (sqrt(variance)/1e-15))

	return t, ETHz, normTHzPower
	
	


plot_1=True
plot_2=False


############### EXAMPLE 1: experimental data for reference 7THz from December
if plot_1:

	THzFilename='purged-7THz-realigned.txt'
	t_offset=0.079e-12-4e-16
	t, y, normTHzPower = getdata(THzFilename,t_offset=t_offset,normalise=False)
	
	#### take the FFT including phase correction with t0 as the time from the start of the sampling to the peak of the pulse envelope
	out = takeFFT(t,y, fmax_THz=8.0)

	##### Do EOS response correction by dividing out the frequency-dependent response of the EO crystal.
	resp = EOSresponse(out.omega,lamGate=800e-9,l=200e-6)

	freqMask=abs(out.freq_THz)<7.6
	ETHzSpectrum_corrected=(out.y_fft/resp.response)[freqMask]
	
	masked_frequency=out.freq_THz[freqMask]

	plt.figure()
	ax1=plt.subplot(211)
	ax1.plot(out.time_ps-out.tgroup/1e-12,out.y)
	ax1.plot(0.0,0.0,'ko')
#	ax1.plot(out.tgroup/1e-12,0.0,'ko')

	ax1.set_xlabel('Time (ps)')
	ax1.set_ylabel('$E(t)$ (arb. units)')	
	ax1.set_xlim(t[0]/1e-12,t[-1]/1e-12)

	ax2=plt.subplot(223)
	ax2.semilogy(out.freq_THz, out.abs_y_fft)
	ax2.semilogy(masked_frequency, abs(ETHzSpectrum_corrected))	

	ax2.set_xlabel('Frequency (THz)')
	ax2.set_ylabel('$|E(f)|$ (arb. units)')		

	ax3=plt.subplot(224)
	ax3.plot(out.freq_THz, out.angle_y_fft)
	ax3.plot(masked_frequency, angle(ETHzSpectrum_corrected))	

	ax3.set_ylim(-pi,pi)
	ax3.set_yticks([-pi,-pi/2,0,pi/2,pi])
	ax3.set_yticklabels(['-$\pi$','-$\pi/2$','$0$','$\pi/2$','$\pi$'])
	ax3.set_xlabel('Frequency (THz)')
	ax3.set_ylabel('arg$[E(f)]$ (radians)')		

	plt.tight_layout()




############### EXAMPLE 2: experimental data for quartz from December
if plot_2:

	THzFilename='Quartz-1mm-purge-7THz.txt'
	t_offset=0.079e-12-4e-16
	t, y, normTHzPower = getdata(THzFilename,t_offset=t_offset,normalise=False)
	
	#### take the FFT including phase correction with t0 as the time from the start of the sampling to the peak of the pulse envelope
	out = takeFFT(t,y, fmax_THz=8.0)

	ax1.plot(out.time_ps,out.y)
	ax1.plot(out.tgroup/1e-12,0.0,'ko')
#	ax1.plot(out.time_ps-out.tgroup/1e-12,out.y)	
	ax2.semilogy(out.freq_THz, out.abs_y_fft)
	ax3.plot(out.freq_THz, out.angle_y_fft)



plt.show()
