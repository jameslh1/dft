from numpy import pi, angle, min, max, abs
from pylab import figure, axis, subplot, tight_layout


def plotTimeDomain(axis,x,y,xunit='ps',**kwargs):
	if xunit=='s':
		sf=1.0
	elif xunit=='ps':
		sf=1.0/1e-12
		
	axis.plot(x*sf,y,**kwargs)
	axis.plot((x[0]*sf,x[-1]*sf),(0,0),'k',linewidth=0.5)

	axis.set_xlabel('Time (%s)' % (xunit))
	axis.set_ylabel('$E_{THz}(t)$ (Vm$^{-1}$)')	

	return
	
	
def plotSpectrum(axis,fftdata,option='amplitude',**kwargs):
	if option=='amplitude':
		axis.plot(fftdata.freq_THz,fftdata.abs_y_fft,**kwargs)
		axis.set_ylabel('$|E(\omega)|$ (Vm$^{-1}$)')
		axis.set_yscale('log')

	elif option=='phase':
		axis.plot(fftdata.freq_THz,fftdata.angle_y_fft,**kwargs)
		axis.set_ylabel('$<E(\omega)$ (radians)')	

		axis.plot((fftdata.freq_THz[0],fftdata.freq_THz[-1]),(0,0),'k',linewidth=0.5)
		
		axis.set_yticks([-pi,-pi/2,0,pi/2,pi])
		axis.set_yticklabels(['$-\pi$','$-\pi/2$','0','$\pi/2$','$\pi$'])			
		axis.set_ylim((-pi,pi))

	axis.set_xlabel('Frequency (THz)')

	return


def plotBaseline(axis,fft1,fft2):
	ratio=fft2.y_fft/fft1.y_fft

	l1,=axis.plot(fft1.freq_THz,abs(ratio))
	axis.set_ylabel('$|R(\omega)|=|E_2(\omega)/E_1(\omega)|$',color=l1.get_color())

	new_ax=axis.twinx()
	l2,=new_ax.plot(fft1.freq_THz,angle(ratio),'orange')
	new_ax.set_ylabel('$<R(\omega)$ (radians)',color=l2.get_color())	

	axis.plot((fft1.freq_THz[0],fft1.freq_THz[-1]),(0,0),'k',linewidth=0.5)

	axis.set_xlabel('Frequency (THz)')

	axis.set_ylim(0,1.5)
	new_ax.set_ylim(0,1.5)
	new_ax.set_yticks([0,pi/4,pi/2])
	new_ax.set_yticklabels(['0','$\pi/4$','$\pi/2$'])	

	return




def addLabel(axis,text,position='upper left'):
	if 'upper' in position:
		ypos=0.9
	else:
		ypos=0.1
		
	if 'left' in position:
		xpos=0.1
	else:
		xpos=0.9

	axis.text(xpos, ypos, text, horizontalalignment='center',
				verticalalignment='center', transform=axis.transAxes)
				
				
				
def plotDualBeamSpectra(data,spectrum_fmax=6.7,ratio_fmax=5.5,ETHz_y_unit='Vm$^{-1}$'):
	figure(figsize=(12,6))
	ax1=subplot(231)
	ax2=subplot(234)
	ax3=subplot(232)
	ax4=subplot(233)
	ax5=subplot(235)
	ax6=subplot(236)

	plotTimeDomain(ax1,data.t1,data.ETHz1,label='$E_1$')
	plotTimeDomain(ax1,data.t2,data.ETHz2,label='$E_2$')
	plotSpectrum(ax3,data.beam1,option='amplitude',label='$E_1$')
	plotSpectrum(ax3,data.beam2,option='amplitude',label='$E_2$')
	plotSpectrum(ax4,data.beam1,option='phase',label='$E_1$')
	plotSpectrum(ax4,data.beam2,option='phase',label='$E_2$')

	plotTimeDomain(ax2,data.t1_shifted,data.ETHz1_shifted,label='$E_1^*$')
	plotTimeDomain(ax2,data.t2_shifted,data.ETHz2_shifted,label='$E_2^*$')
	#plotSpectrum(ax5,beam1_shifted,option='amplitude')
	#plotSpectrum(ax5,beam2_shifted,option='amplitude')
	plotSpectrum(ax5,data.beam1_shifted,option='phase',label='$E_1^*$')
	plotSpectrum(ax5,data.beam2_shifted,option='phase',label='$E_2^*$')
	plotBaseline(ax6,data.beam1_shifted,data.beam2_shifted)


	#ax1.set_title('Raw data')
	ax1.set_xlim(-2,4)
	ax2.set_xlim(-0.7,1)
	ax3.set_xlim(0,spectrum_fmax)
	ax4.set_xlim(0,spectrum_fmax)
	ax5.set_xlim(0,ratio_fmax)
	ax6.set_xlim(0,ratio_fmax)

	ax1.set_ylabel('$E(t)$ (%s)'% (ETHz_y_unit))
	ax2.set_ylabel('$E(t)$ (%s)'% (ETHz_y_unit))	
	ax3.set_ylabel('$|E(\omega)|$ (%s)' % (ETHz_y_unit))
	ax5.set_ylabel('$|E(\omega)|$ (%s)' % (ETHz_y_unit))	

	addLabel(ax1,'(a)',position='upper left')
	addLabel(ax3,'(b)',position='upper right')
	addLabel(ax4,'(c)',position='upper left')
	addLabel(ax2,'(d)',position='upper left')
	addLabel(ax5,'(e)',position='upper left')
	addLabel(ax6,'(f)',position='upper left')

	ax1.legend()
	ax2.legend()
	ax3.legend(loc='lower left')
#	ax4.legend()
	ax5.legend()
	#ax1.legend()


	#ax2.set_title('Time shifted')
	tight_layout()

	all_axes=[ax1,ax3,ax4,ax2,ax5,ax6]
	
	return all_axes

