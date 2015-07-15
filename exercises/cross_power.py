import numpy as np
import matplotlib.pyplot as plt
import time

start_time = time.time()

# Generate data: sin(2pi*[f in Hz]) + noise
time_step = 1e-4
t = np.arange(0., 1., time_step)
sig1 = np.sin(2*np.pi*3*t) + np.random.normal(size=len(t))
sig2 = np.sin(2*np.pi*7*t) + np.random.normal(size=len(t))

# a_corr <=> invfft(fft*conj(fft))

fft1 = np.fft.fft(sig1)
fft2 = np.fft.fft(sig2)

# Cross-power spectrum is the fourier transform of cross correlation
# x_corr <=> invfft(fft1*conj(fft2))
PS = np.abs(fft1*np.conj(fft2))
# Get frequency bin centers in cycles per unit of sample spacing
freqs = np.fft.fftfreq(sig1.size, time_step)
# Argument order corresponding to ascending frequencies
idx = np.argsort(freqs)

print 'Computation time: %.4f seconds' % (time.time() - start_time)

plt.figure()
plt.suptitle('Cross-power spectrum')

plt.subplot(211)
plt.plot(t, sig1, 'r')
plt.plot(t, sig2, 'b')
plt.xlabel('Time (s)')
plt.ylabel('Signal (V)')

plt.subplot(212)
plt.plot(freqs[idx], PS[idx])
plt.xlim([0, 10])
plt.xlabel('Frequency (Hz)')
plt.ylabel('Cross-poewr spectrum')

plt.show()

