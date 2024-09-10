from numpy import pi, linspace, exp, gradient, cosh, heaviside, sin, real, imag, angle, unwrap

import matplotlib.pyplot as plt
plt.style.use('nature')

##### import the DFT class library
import dft
from dft import i

plot_example_1=True
plot_example_2=True

############### EXAMPLE 1: Gaussian centred at time zero and a smaller one at 1.0ps. Valid for zero CEP
if plot_example_1:

	n=1000						# number of discrete samples
	t_start = -2e-12			# start of the time window
	t_end = 10e-12				# end of the time window

	t=linspace(t_start,t_end,n)

	#### pulse parameters
	t_pulse1 = 0e-12				# centre of the simulated pulse 1
	t_pulse2 = 1e-12				# centre of the simulated pulse 1	
	sigma=0.1e-12

	y=exp(-(t-t_pulse1)**2  / (2.0*sigma**2))			# simulated pulse
	y2=0.8*exp(-(t-t_pulse2)**2  / (2.0*sigma**2))		# simulated pulse	

	y=gradient(y)
	y2=gradient(y2)


	#### take the FFT including phase correction with t0 as the time from the start of the sampling to the peak of the pulse envelope
	out = dft.takeFFT(t,y, fmax_THz=8.0)
	out2 = dft.takeFFT(t,y2, fmax_THz=8.0)

	#### return the time-domain data from the FFT data
	input_omega=out.full_omega		 # omega
	input_E_omega=out.full_y_fft	 # E(omega)
	t0=out.t0						 # need to know this to get Ek back (from E(omega))
	dt=out.dt						 # need to know this to get amplitude scaling correct
	E_THz_recovered=dft.takeIFFT(input_omega,input_E_omega,t0=t0,dt=dt)


	plt.figure(figsize=(12,6))
	ax1=plt.subplot(211)
	ax1.plot(out.time_ps,out.y,label='model $E_1(t)$')
	ax1.plot(out2.time_ps,out2.y,label='model $E_2(t)$')
	
	ax1.plot(out.time_ps,E_THz_recovered,'r.',label='DFT$^{-1}$[$E_1(\omega)$]')	
	ax1.plot(out.tgroup/1e-12,0.0,'o',color='blue')
	ax1.plot(out2.tgroup/1e-12,0.0,'o',color='orange')
	
	ax1.set_xlabel('Time (ps)')
	ax1.set_ylabel('$E(t)$ (arb. units)')	
	ax1.set_xlim(-1,3)
	ax1.legend(loc=1)
	

	ax2=plt.subplot(223)
	ax2.semilogy(out.freq_THz, out.abs_y_fft)
	ax2.plot(out2.freq_THz, out2.abs_y_fft)	
	ax2.set_xlabel('Frequency (THz)')
	ax2.set_ylabel('$|E(f)|$ (arb. units)')		

	ax3=plt.subplot(224)
	ax3.plot(out.freq_THz, out.angle_y_fft)
	ax3.plot(out2.freq_THz, out2.angle_y_fft)	
	ax3.set_ylim(-pi,pi)
	ax3.set_yticks([-pi,-pi/2,0,pi/2,pi])
	ax3.set_yticklabels(['-$\pi$','-$\pi/2$','$0$','$\pi/2$','$\pi$'])
	ax3.set_xlabel('Frequency (THz)')
	ax3.set_ylabel('arg$[E(f)]$ (radians)')		

	plt.tight_layout()


############### EXAMPLE 2: Sech pulse with finite CEP.
if plot_example_2:

	n=2000						# number of discrete samples
	t_start = -20e-12			# start of the time window
	t_end = 20e-12				# end of the time window

	t=linspace(t_start,t_end,n)

	#### pulse parameters
	t_pulse1 = 0e-12				# centre of the simulated pulse 1
	t_pulse2 = 0e-12				# centre of the simulated pulse 1	
	sigma=0.1e-12					# pulse duration 100fs

	omega_c=2*pi*1e12				# angular carrier frequency
	tCEP=125e-15					# CEP time: 125fs is 1/8th of 1ps period, hence CEP=pi/4.

	y =1.0/cosh((t-t_pulse1)/sigma)*-1.0 * sin(omega_c*(t-tCEP))	# incorrect approach, putting tCEP into sin()
	y2=1.0/cosh((t-t_pulse2)/sigma)*-1.0 * sin(omega_c*t)			# correct approach, used to construct analytical representation of E(t)

	y=y/max(y)		# normalise the electric field
	y2=y2/max(y2)	# normalise the electric field

	#### take the FFT including phase correction with t0 as the time from the start of the sampling to the peak of the pulse envelope
	out  = dft.takeFFT(t,y , fmax_THz=8.0)
	out2 = dft.takeFFT(t,y2, fmax_THz=8.0)

	#### return the time-domain data from the FFT data
	input_omega=out.full_omega		 # omega
	input_E_omega=out.full_y_fft	 # E(omega)
	input_E_omega2=out2.full_y_fft	 # E(omega)	

	#### apply a mask to set all negative frequency components to zero
	mask=heaviside(out.full_omega,0.0)

	output_E_omega =2*mask*input_E_omega
	output_E_omega2=2*mask*input_E_omega2

	phi=omega_c*tCEP	
	t0=out.t0						 # need to know this to get Ek back (from E(omega))
	dt=out.dt						 # need to know this to get amplitude scaling correct
	E_THz_recovered =dft.takeIFFT(input_omega,output_E_omega,t0=t0,dt=dt)
	E_THz_recovered =real(E_THz_recovered)
	E_THz_recovered2=dft.takeIFFT(input_omega,output_E_omega2,t0=out2.t0,dt=out2.dt)
	E_THz_recovered2=real(E_THz_recovered2*exp(i*phi))			# change the phase

	fftOut  = dft.takeFFT(t,E_THz_recovered, fmax_THz=8.0)
	fftOut2 = dft.takeFFT(t,E_THz_recovered2, fmax_THz=8.0)	

	plt.figure(figsize=(12,6))
	ax1=plt.subplot(211)
#	ax1.plot(out.time_ps ,out.y ,label='model $E_1(t)$')
#	ax1.plot(out2.time_ps,out2.y,label='model $E_2(t)$')
	
	ax1.plot(out.time_ps ,E_THz_recovered /max(E_THz_recovered) ,label='simple model')	
	ax1.plot(out2.time_ps,E_THz_recovered2/max(E_THz_recovered2),label='correct model')		
	ax1.plot(out.tgroup/1e-12,0.0,'o',color='blue')
	ax1.plot(out2.tgroup/1e-12,0.0,'o',color='orange')
	
	ax1.set_xlabel('Time (ps)')
	ax1.set_ylabel('$E(t)$ (arb. units)')	
	ax1.set_xlim(-1,3)
	ax1.legend(loc=1)
	ax1.set_title('Incorrect (simple) CEP model v correct CEP model. sech envelope, 100fs, 1THz carrier freq.')

	ax2=plt.subplot(223)
#	ax2.semilogy(out.freq_THz, out.abs_y_fft)
#	ax2.plot(out2.freq_THz, out2.abs_y_fft)	

	ax2.semilogy(fftOut.freq_THz, fftOut.abs_y_fft,label='simple model')	
	ax2.plot(fftOut2.freq_THz, fftOut2.abs_y_fft,label='correct model')	
	
	ax2.set_xlabel('Frequency (THz)')
	ax2.set_ylabel('$|E(f)|$ (arb. units)')		
	ax2.legend(loc=1)


	ax3=plt.subplot(224)
#	ax3.plot(out.freq_THz, out.angle_y_fft,)
#	ax3.plot(out2.freq_THz, out2.angle_y_fft)
	
	ax3.plot(fftOut.freq_THz, fftOut.angle_y_fft,label='simple model')
	ax3.plot(fftOut2.freq_THz, fftOut2.angle_y_fft,label='correct model')
	
	ax3.set_ylim(-pi,pi)
	ax3.set_yticks([-pi,-pi/2,0,pi/2,pi])
	ax3.set_yticklabels(['-$\pi$','-$\pi/2$','$0$','$\pi/2$','$\pi$'])
	ax3.set_xlabel('Frequency (THz)')
	ax3.set_ylabel('arg$[E(f)]$ (radians)')		

	ax3.legend(loc=1)

	plt.tight_layout()


plt.show()
