import dft

from pylab import *

freq_THz=dft.linspace(0.001,20,1000)

freq=freq_THz*1e12
omega=2*dft.pi*freq

r=dft.EOSresponse(omega,lamGate=800e-9,l=200e-6)
#r2=dft.EOSresponse(omega,lamGate=800e-9,l=220e-6)


ax1=subplot(111)
plot(freq_THz,abs(r.response))
#plot(freq_THz,abs(r2.response),color='darkblue',linestyle='--')

xlim(0,10)
ylabel('$|G(\omega)r_{41}(\omega)|$',color='darkblue')
xlabel('Frequency (THz)')
title('GaP, 200$\mu$m thick')


ax2=ax1.twinx()
plot(freq_THz,angle(r.response),'orange')
#plot(freq_THz,angle(r2.response),'r--')

ylabel(r'$\angle G(\omega)r_{41}(\omega)$ (rad.)',color='orange')


show()
