from numpy import fft,array,gradient, append,insert, nonzero, ravel, pi, argmax, unwrap, angle, mean, asarray, abs, argmin, exp, trapz, real
from scipy.optimize import leastsq
#from generic import *
################################################

C = 299792458.0
PS_TO_MM = 0.1499

hbar = 1.05457266E-34
h = hbar * 2 * pi
eV = 1.60217733E-19

i=complex(0,1)


def find(condition):
    res, = nonzero(ravel(condition))
    return res



def find_nearest(array, value):
	array = asarray(array)
	idx = (abs(array - value)).argmin()
	return idx, array[idx]


# extend length of "xdata" & "ydata" by a scale factor of "expansion"
def zeropad(xdata,ydata,expansion=2):
	if expansion>1:
#		print('zeropadded')
#		spacing=xdata[2]-xdata[1]		# old version
#		spacing=xdata[1]-xdata[0]		# option 1: fine for equally spaced xdata
		spacing=(xdata[-1]-xdata[0])/(len(xdata)-1)	# option 2: better if xdata is not equally spaced
		newx = xdata
		newy = ydata

		lastPoint=len(xdata)
		x0=xdata[-1]

		for x in range(1,expansion*len(xdata)):
			newx=append(newx,x0+x*spacing)
			newy=append(newy,0.0)

	else:
		newx = xdata
		newy = ydata


	return newx, newy



def integrateTHzPower(t,ETHz):
	integral=trapz(ETHz**2,x=t)

	return integral


def meanTime(t,ETHz):
	#### calculate the first moment of the intensity v time to work out the mean arrival time of the pulse.
	#  This method works okay if there is one main pulse in the time-domain data set.

	I0=trapz(ETHz**2,x=t)	# work out the normalisation for the intensity
	I =ETHz**2 / I0			# work out the intensity, normalised such that the integral w.r.t. time is 1
	
	meanTime=trapz(t*I,x=t)	# the mean time of the pulse is the first moment of the intensity, i.e. \int t I(t) dt.

	return meanTime			# return the mean time
	

def timeShift(fftdata,tshift=0.0e-12):

	#### return the time-domain data from the FFT data
	input_omega  =fftdata.full_omega	 # omega
	input_E_omega=fftdata.full_y_fft	 # E(omega)
	t1=fftdata.t0						 # need to know this to get Ek back (from E(omega))
		
	dt=fftdata.dt						 # need to know this to get amplitude scaling correct
	y_shifted=takeIFFT(input_omega,input_E_omega,t0=t1-tshift,dt=dt)

	t=fftdata.time

	return t, y_shifted






def takeIFFT(omega,E_omega,t0=0.0e-12,dt=0.01e-12):
    ###### WORK OUT Ek
    Ek = E_omega * exp(i*omega*t0)      
    # Note: - input to DFT routine needs to have phase corrected
    
    ###### SHIFT THE DFT 
    Ek = fft.ifftshift(Ek)
    # Note: this puts the sequence back into the right order
    
    Em = fft.fft(Ek)/(len(omega)*dt)
    
#    return real(Em)	# do we need to return the real part only?
    return Em


class takeFFT():
	"""
	class takeFFT() : usage: 
	  out = takeFFT(t,y)
	  -> returns class instance "out" with properties
	  		t = a copy of t
	  		y = a copy of y
	        n = number of points
	        dt= time step
	        freq_THz = frequency in THz
	        omega = angular frequency
			y_fft = array containing complex Fourier spectrum of y
			
	  **kwarg options: 
	   - t0 = 1.0e-12    -> sets t0 to 1ps
	   - fmax_THz = 7.0  -> sets the maximum frequency returned to 7THz
	   - xunit = {'s', 'mm', 'ps'} if "time" input is in s, mm or ps, respectively
	   - ret_neg_freq = {True, False [default]} 
	            -> if true, returns all discretely sampled frequencies from -fmax_THz to fmax_THz.
	            -> if false, returns all frequencies from 0 to fmax_THz.	            
	
	"""


	def __init__(self,time,ydata,fmax_THz=None,xunit='s', ret_neg_freq=False, t_offset=0.0, find_max_freq=True, assumedCarrierFreq=1e12,group_delay_option=0):
#		if xunit=='mm':
#			sf=1.0/0.1499
#			xlabel='Frequency (THz)'
#		el
		if xunit=='ps':
			sf=1e-12
			xlabel='Frequency (THz)'
		else:
			sf=1
			xlabel='Frequency (Hz)'	
			
		self.xlabel=xlabel
		self.time=time*sf			# time in SI units (s)
		self.time_ps=self.time/1e-12
		
		self.assumedCarrierFreq=assumedCarrierFreq	
		self.group_delay_option=group_delay_option
		
		self.n = len(time)
		self.dt= abs(self.time[1]-self.time[0])
		self.t0= t_offset-self.time[0]		# assume that the time-domain pulse envelope peaks at the time difference between 0.0 and the time of the first sample

		self.fmax_THz=fmax_THz
		

		####### CALCULATE THE DFT USING NUMPY
		fftdata = self.n * fft.ifft(ydata,n=self.n)   
		# Notes: - the ifft routine is used as it uses exp(i omega t).
		#        - numpy's ifft divides by n, so need to scale by n to be
		#          able to use this as the forward transform.

		####### WORK OUT THE FREQUENCY ELEMENTS
		freq=fft.fftfreq(n=self.n, d=self.dt)	
		# Notes: - freq runs from 0 -> f_max then -f_max -> -1/T, 
		#          where f_max = 1/dt for time step "dt".

		####### SHIFT THE DFT
		freq = fft.fftshift(freq)
		fftdata = fft.fftshift(fftdata)
		# Notes: - this makes sequence run from -f_max to +f_max 
		#          so that it is symmetric about zero frequency, like the CFT

		omega = 2*pi*freq 

		####### CORRECT THE DFT SO THAT IT AGREES WITH THE CFT
		fftdata = fftdata * self.dt * exp(-i * omega * self.t0)
		# Notes: - scaling factor "dt" needed for correct amplitude
		#        - phase factor exp(-i * omega * t0) accounts for pulse
		#          centred at time t0

		####### STORE THE USEFUL BITS IN THE CLASS INSTANCE
		self.t=1.0*time
		self.y=1.0*ydata
		self.full_y_fft= 1.0*fftdata		
		self.full_omega= 1.0*omega

		freq_THz = omega/(2*pi*1e12)

		if self.fmax_THz==None:
			self.freq_THz=freq_THz
			self.y_fft=fftdata			
		else:
			if ret_neg_freq:
				fmask=abs(freq_THz)<=self.fmax_THz
				self.freq_THz = freq_THz[fmask]
				self.y_fft = fftdata[fmask]
			else:
				fmask1=freq_THz>0.0
				fmask2=freq_THz[fmask1]<=self.fmax_THz			
				self.freq_THz=freq_THz[fmask1][fmask2]
				self.y_fft = fftdata[fmask1][fmask2]

		self.freq = self.freq_THz * 1e12	
		self.omega = 2*pi*self.freq
		self.abs_y_fft=abs(self.y_fft)
		self.angle_y_fft=angle(self.y_fft)

		self.get_tgroup()

	
		
	def get_tgroup(self):
		self.t_group_delay_dispersion = gradient(self.angle_y_fft, self.omega)


		if self.group_delay_option==1:		########### method 1: find the maximum of the THz pulse and take the group delay from that
			self.index_max_freq=argmax(self.abs_y_fft)			# index of the maximum frequency
			self.carrier_freq_THz=self.freq_THz[self.index_max_freq]
			self.tgroup=self.t_group_delay_dispersion[self.index_max_freq]

		elif self.group_delay_option==2:	########### method 2: find the group delay at a specific frequency
			idx,value=find_nearest(self.freq_THz*1e12,self.assumedCarrierFreq)
			self.index_max_freq=idx
			self.carrier_freq_THz=self.freq_THz[self.index_max_freq]
			self.tgroup=self.t_group_delay_dispersion[self.index_max_freq]

		elif self.group_delay_option==3:	########### method 3: find the group delay by averaging over a range of frequencies
			idx1,value=find_nearest(self.freq_THz*1e12,0.4e12)
			idx2,value=find_nearest(self.freq_THz*1e12,1.0e12)
			self.index_max_freq=int((idx2+idx1)/2)
			self.carrier_freq_THz=self.freq_THz[self.index_max_freq]
			self.tgroup=self.t_group_delay_dispersion[self.index_max_freq]

		elif self.group_delay_option==4:	########### method 4: find the group delay from the times of the positive and negative peak
			index_max_pos_time=argmax(self.y)
			index_max_neg_time=argmax(-self.y)
			index_max_time=int(round((index_max_pos_time+index_max_neg_time)/2.0))
			self.carrier_freq_THz=1.0/(abs(self.t[index_max_pos_time]-self.t[index_max_neg_time])*2)	# t_peak to peak is half of the carrier period, hence x2
			self.tgroup=self.t[index_max_time]

		else:		####### default method: use the mean time	
			self.tgroup=meanTime(self.t,self.y)

		if self.group_delay_option==1:
			print('-- carrier freq = %4.3fTHz' % (self.carrier_freq_THz))
			print('-- t_group(at %4.3fTHz) = %4.3ffs' % (self.carrier_freq_THz,self.tgroup/1e-15))
		elif self.group_delay_option==2:
			print('-- carrier freq = %4.3fTHz' % (self.assumedCarrierFreq/1e12))
			print('-- t_group(at %4.3fTHz) = %4.3ffs' % (self.assumedCarrierFreq/1e12,self.tgroup/1e-15))
		else:
			print('-- t_group(from time domain) = %4.3ffs' % (self.tgroup/1e-15))
#			print('-- "carrier freq" = 1/t_pkpk = %4.3fTHz' % (self.carrier_freq_THz/1e12))




class DualBeamData():
	def __init__(self,t1,ETHz1,t2,ETHz2,GDO=0,fmax_THz=7):
	
		self.t2=t2
		self.ETHz2=ETHz2
		self.t1=t1
		self.ETHz1=ETHz1

		#### take the FFT including phase correction with t0 as the time from the start of the sampling to the peak of the pulse envelope
		print('###### beam 1')
		self.beam1 = takeFFT(t1,ETHz1, fmax_THz=fmax_THz, group_delay_option=GDO)
		print('###### beam 2')
		self.beam2 = takeFFT(t2,ETHz2, fmax_THz=fmax_THz, group_delay_option=GDO)

		self.t1_shifted,self.ETHz1_shifted=timeShift(self.beam1,tshift=self.beam1.tgroup)
		self.t2_shifted,self.ETHz2_shifted=timeShift(self.beam2,tshift=self.beam2.tgroup)

		print('###### beam 1 shifted')
		self.beam1_shifted = takeFFT(self.t1_shifted,self.ETHz1_shifted, fmax_THz=6,group_delay_option=GDO)#, assumedCarrierFreq=1.5e12)
		print('###### beam 2 shifted')
		self.beam2_shifted = takeFFT(self.t2_shifted,self.ETHz2_shifted, fmax_THz=6,group_delay_option=GDO)#, assumedCarrierFreq=1.5e12)

