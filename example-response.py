import dft

from pylab import *

freq_THz=dft.linspace(0.001,20,1000)

freq=freq_THz*1e12
omega=2*dft.pi*freq

r=dft.EOSresponse(omega,lamGate=800e-9,l=200e-6)
r2=dft.EOSresponse(omega,lamGate=800e-9,l=220e-6)


ax1=subplot(111)
plot(freq_THz,abs(r.responseFn),'darkblue')
plot(freq_THz,abs(r2.responseFn),color='darkblue',linestyle='--')

xlim(0,6)

ax2=ax1.twinx()
plot(freq_THz,angle(r.responseFn),'r')
plot(freq_THz,angle(r2.responseFn),'r--')



show()
