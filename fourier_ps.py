import numpy as np
import matplotlib.pyplot as plt
import time


def power_analysis(signal, PS):
    pspower = 1./(signal.size)*np.sqrt(np.sum(PS))
    varpower = np.sqrt(np.var(signal))
    print '    sqrt of total PS(signal):'
    print '        1/N*sqrt(sum(PS(signal)))\t= ', pspower
    print '    sqrt of variance, rms fluct level, signal:'
    print '        sqrt(variance(signal))\t\t= ', varpower


start_time = time.time()

# Time
time_step = 1e-8
t = np.arange(0., 1., time_step)

# Signal 1: sin(2pi*[f in Hz])
freq1 = 300
signal1 = np.sin(2*np.pi*freq1*t)
PS1 = np.abs(np.fft.fft(signal1))**2
print '1V amplitude sine wave at %d Hz:' % freq1
power_analysis(signal1, PS1)

# Signal 2: noise
signal2 = np.random.normal(size=len(t))
PS2 = np.abs(np.fft.fft(signal2))**2
print '1V amplitude noise (normally distributed):'
power_analysis(signal2, PS2)

# Signal 3: sine + noise
signal3 = signal1 + signal2
PS3 = np.abs(np.fft.fft(signal3))**2
print '1V sine wave + 1V noise:'
power_analysis(signal3, PS3)

print 'Computation time: %.4f seconds' % (time.time() - start_time)

# Get frequency bin centers in cycles per unit of sample spacing
freqs = np.fft.fftfreq(signal1.size, time_step)
# Argument order corresponding to ascending frequencies
idx = np.argsort(freqs)

plt.figure()

plt.subplot(231).set_title('sin @ %d Hz' % freq1)
plt.plot(t, signal1)
plt.xlabel('Time (s)')
plt.ylabel('Signal (V)')
plt.xlim([0, 1./freq1])
plt.subplot(234)
plt.plot(freqs[idx], PS1[idx])
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power spectrum')
plt.xlim([0, freq1+10])

plt.subplot(232).set_title('noise')
plt.plot(t, signal2)
plt.xlabel('Time (s)')
plt.subplot(235)
plt.plot(freqs[idx], PS2[idx])
plt.xlabel('Frequency (Hz)')

plt.subplot(233).set_title('sin + noise')
plt.plot(t, signal3)
plt.xlabel('Time (s)')
plt.xlim([0, 1./freq1])
plt.subplot(236)
plt.plot(freqs[idx], PS3[idx])
plt.xlabel('Frequency (Hz)')
plt.xlim([0, freq1+10])

plt.show()
