from numpy import pi, linspace, exp, gradient

import matplotlib.pyplot as plt
plt.style.use('nature')

##### import the DFT class library
import dft


plot_example_1=True


############### EXAMPLE 1: Gaussian centred at time zero and a smaller one at 1.0ps
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
	t1=out.t0						 # need to know this to get Ek back (from E(omega))
		
	dt=out.dt						 # need to know this to get amplitude scaling correct
	E_THz_recovered=dft.takeIFFT(input_omega,input_E_omega,t0=t1,dt=dt)

	#### return the time-domain data from the FFT data
	input_omega2=out2.full_omega		 # omega
	input_E_omega2=out2.full_y_fft	 # E(omega)
	t2=out2.t0						 # need to know this to get Ek back (from E(omega))
		
	dt2=out2.dt						 # need to know this to get amplitude scaling correct
	tshift=out2.tgroup
	E_THz_recovered2=dft.takeIFFT(input_omega2,input_E_omega2,t0=t2-tshift,dt=dt2)
	print(out2.tgroup)

	plt.figure(figsize=(12,6))
	ax1=plt.subplot(211)
	ax1.plot(out.time_ps,out.y,label='model $E_1(t)$')
	ax1.plot(out2.time_ps,out2.y,label='model $E_2(t)$')
	
	ax1.plot(out.time_ps,E_THz_recovered,'b.',label='DFT$^{-1}$[$E_1(\omega)$]')	
	ax1.plot(out2.time_ps,E_THz_recovered2,'.',color='orange',label='DFT$^{-1}$[$E_2(\omega)$]')	
	ax1.plot(out.tgroup/1e-12,0.0,'ko')
	ax1.plot(out2.tgroup/1e-12,0.0,'ko')
	
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

plt.show()
