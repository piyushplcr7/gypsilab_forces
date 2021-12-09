import numpy as np
import matplotlib.pyplot as plt
import sys

sys.path.append("/scratch/users/ppanchal/mpltools/mpltools")
from mpltools import annotation

data = np.loadtxt('sc3_0_0')
step = 1
tbem = data[::step,2]
tbdry = data[::step,4]
print("tbem shape", tbem.shape)
hbem = 28/data[::step,0]
print("hbem shape", hbem.shape)

# Using last value from the egg shell method as converged
#f = fv[-1]
# Using value from BEM computations
#t = -0.00647118123921061 # 5400 panels

t = -0.0567971810652226 # 7168 panels

fig = plt.figure()
ax = fig.add_subplot(111)

plt.loglog(hbem,abs(tbem-t)/abs(t),'.-',color='red',zorder = 2)
plt.loglog(hbem,abs(tbdry-t)/abs(t),'.-',color='orange',zorder = 2)
plt.legend(['Pullback approach','Stress tensor (BEM)'])
plt.xlabel('meshwidth h')
plt.ylabel('Relative error of net torque')

# Fitting lines to these plots
fit_bnd = np.polyfit(np.log(hbem)[1:],np.log(abs(tbdry-t)/abs(t))[1:],1)
fit_bem = np.polyfit(np.log(hbem),np.log(abs(tbem-t)/abs(t)),1)
y_bem = fit_bem[0] * np.log(hbem) + fit_bem[1]
y_bnd = fit_bnd[0] * np.log(hbem) + fit_bnd[1]

bestfitcolor = "silver"
plt.loglog(hbem,np.exp(y_bnd),'--',color=bestfitcolor,zorder = 1)
plt.loglog(hbem,np.exp(y_bem),'--',color=bestfitcolor,zorder = 1)

	
annotation.slope_marker((6e-2, 3e-4), (1.32, 1),ax = ax)

annotation.slope_marker((6e-2,3e-2), (0.31, 1),ax = ax)


print("boundary rate", fit_bnd[0])
print("BEM rate", fit_bem[0])

plt.savefig('sq_torque.eps')
plt.show()
