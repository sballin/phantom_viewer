# Discussion: http://stackoverflow.com/questions/643699/how-can-i-use-numpy-correlate-to-do-autocorrelation

import numpy as np
import matplotlib.pyplot as plt
import time

start_time = time.time()

# Generate data: sin(2pi*[f in Hz]) + noise
time_step = .0001
t = np.arange(0., 1., time_step)
y = np.sin(2*np.pi*10*t) + np.random.uniform(size=len(t)) - 0.5

yunbiased = y-np.mean(y)
ynorm = np.sum(yunbiased**2)
acor = np.correlate(yunbiased, yunbiased, "same")/ynorm
# use only second half
#acor = acor[len(acor)/2:]

print 'Computation time: %.4f seconds' % (time.time() - start_time)

plt.figure()
plt.suptitle('Direct autocorrelation')

plt.subplot(211)
plt.plot(t, y)
plt.xlabel('Time (s)')
plt.ylabel('Signal (V)')

plt.subplot(212)
plt.plot(t, acor)
plt.xlabel('Time offset (s)')
plt.ylabel('Autocorrelation')

plt.show()
