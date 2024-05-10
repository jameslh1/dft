from numpy import fft,array,gradient, append,insert, nonzero, ravel, pi, argmax, unwrap, angle, mean, asarray, abs, argmin, sin, cos, zeros, real, linspace, imag, sqrt, exp, max

i=complex(0,1)
hbar = 1.05457266E-34	# Planck's constant / 2pi
h = hbar * 2 * pi	# Planck's constant
e0 = 8.854187817e-12 	# permittivity of free space epsilon0 
c = 2.99792458e8	# speed of light in vacuum


def find_nearest(array, value):
	array = asarray(array)
	idx = (abs(array - value)).argmin()
	return idx, array[idx]


class EOSresponse():

	def __init__(self,omega,lamGate=800e-9,l=200e-6,fmax_THz=7.68):
		self.fmax=2.3
		self.R=500e-6
		self.l=l
		self.omega=omega	#[len(omega)/2+1:]
		self.lamGate=lamGate

		self.fmax_THz=fmax_THz # need to set a maximum frequency over which to apply the filter (before the pole in the response function at about 10THz for GaP), as otherwise division by zero can lead to artefacts
		
		self.nGr=3.556

		self.nEOcrystal=self.getn()
		self.r41=self.calc_r41()		
		self.calcResponse()


	def getn(self):
#		eel=9.075	# epsilon electronic = epsilon infinity
		eel=9.11	# epsilon electronic = epsilon infinity
		fTO=367.3/33.3 * 1e12
		fLO=403.0/33.3 * 1e12
		omegaTO=2*pi*fTO
		omegaLO=2*pi*fLO
		gamma=2*pi*4.3/33.3 * 1e12
#		epsZnTe=eel + (est*omTO**2/(omTO**2-self.omega**2-2*i*gamTO*self.omega))
#		print(omegaTO/(2*pi*1e12),omegaLO/(2*pi*1e12),gamma/(2*pi*1e12))
		A = (hbar*omegaLO)**2 - (hbar*omegaTO)**2
		B = (hbar*omegaTO)**2 - (hbar*self.omega)**2 - i * hbar**2 * gamma * self.omega
		frac = A/B
		nGaP = sqrt((1.0 + frac)*eel)
#		nZnTe=sqrt(epsZnTe)
		self.complexn=nGaP
		
		return nGaP


	def calc_r41(self):
#		re = 1.0e-12		# re = 1e-12 m/V in Casalbuoni
		re =5.0/2.2/1.075
		C = -0.53		# Leitenstorfer has C=-0.47
		
#		eel=9.11	# epsilon electronic = epsilon infinity
		eel=8.7		# epsilon electronic = epsilon infinity
		fTO=367.3/33.3 * 1e12
		fLO=403.0/33.3 * 1e12
		omegaTO=2*pi*fTO
		omegaLO=2*pi*fLO
		gamma=2*pi*4.3/33.3 * 1e12	
	
#		r41 = re * (1.0 + C / (1.0 - ((hbar*self.omega)**2-i*hbar*self.omega*gamma)/(hbar**2 * omegaTO**2)))
#		r41 = re * (1.0 + C * omegaTO**2 / (omegaTO**2 - self.omega**2 - i * gamma * self.omega))	# Casalbuoni
	
		r41 = re * (1.0 + C * omegaTO**2 / (omegaTO**2 - self.omega**2 + i * gamma * self.omega))	# Casalbuoni (opposite sign imaginary part)	
#		r41=r41/max(abs(r41))	
	
		return r41


	def calcResponse(self):
	
#		A = exp(-i*2*pi* self.omega* self.l*(self.nGr-self.complexn)/c) - 1.0
#		B = -i * 2 * pi* self.omega* self.l*(self.nGr - self.complexn)
#		response = 2.0/(self.complexn + 1.0) * c* A/B

		n=self.nEOcrystal	

		om0=2*pi*c/self.lamGate
		vg=c/self.nGr		# for GaP at 800nm

		dkP = n*self.omega/c - self.omega/vg	# omega=ang. freq at THz frequencies

		response=zeros((len(dkP)))
		omega_mask1=self.omega>0
		omega_mask2=self.omega<0		

		middle=len(response)//2		# //2 to floor the index

		response=(exp(i*dkP*self.l)-1.0)/(i*dkP*self.l)		# added an extra l on the bottom c.f. Gallot's result in order to have the dispersion function
#		response=(exp(i*dkP*self.l/2)-exp(-i*dkP*self.l/2))/(i*dkP*self.l)		# added an extra l on the bottom c.f. Gallot's result in order to have the dispersion function. Symmetrised it to be more like sin(dkP L)/(dkP L)

		tij=2.0/(self.complexn+1.0)
		response=response*tij


#		response[middle]=1.0+0.0*i
#		response[middle]=1.0+0.0*i

		freqTHz=self.omega/(2*pi*1e12)		
		idx,value=find_nearest(freqTHz,0.0)
		response[idx]=response[idx+1]

		response=response/max(abs(response[idx+2]))
	

		##### Since the "geometric" response function - group velocity mismatch term as defined in Gallot - has a slope to the phase, work out what this is so that we can correct for the time-shift this linear phase introduces.
		gdd_response=gradient(angle(response),self.omega)
		idx,value=find_nearest(freqTHz,1.0)
		tcorrect=gdd_response[idx]
		response=response * exp(-i*self.omega*tcorrect)		# correct the linear phase term so that it doesn't give a time shift.

		fmask1=freqTHz>self.fmax_THz
		fmask2=freqTHz<-self.fmax_THz

		response[fmask1]=1.0+0.0*i			
		response[fmask2]=1.0+0.0*i
		self.r41[fmask1]=1.0+0.0*i
		self.r41[fmask2]=1.0+0.0*i		
		
#		response=response*self.r41			# total response function including r41
				

		total_response=response*self.r41	# total response function including r41
		
#		total_response[fmask1]=1.0e3+0.0*i
#		total_response[fmask2]=1.0e3+0.0*i		
		
		self.geometric_response=response		# group velocity mismatch only
		self.response=total_response

#		print(max(abs(self.geometric_response)))		
#		print(max(abs(self.r41)))		
#		self.response=self.response/max(abs(self.response))	

