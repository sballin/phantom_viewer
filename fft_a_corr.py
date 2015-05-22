import numpy as np
import matplotlib.pyplot as plt
import time

start_time = time.time()

# Generate data: sin(2pi*[f in Hz]) + noise
time_step = .0001
t = np.arange(0., 1., time_step)
data = np.sin(2*np.pi*10*t) + np.random.uniform(size=len(t)) - 0.5

# a_corr <=> invfft(fft*conj(fft))

fft = np.fft.fft(data)
a_corr = np.fft.ifft(fft*np.conj(fft))

print 'Computation time: %.4f seconds' % (time.time() - start_time)

plt.figure()
plt.suptitle('Autocorrelation via FFT')

plt.subplot(211)
plt.plot(t, data)
plt.xlabel('Time (s)')
plt.ylabel('Signal (V)')

plt.subplot(212)
plt.plot(t, a_corr.real)
plt.xlabel('Time offset (s)')
plt.ylabel('Autocorrelation')

plt.show()
