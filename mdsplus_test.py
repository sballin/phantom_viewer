import numpy as np
import matplotlib.pyplot as plt
import MDSplus


def power_analysis(signal, PS):
    pspower = 1./(signal.size)*np.sqrt(np.sum(PS))
    varpower = np.sqrt(np.var(signal))
    print '    sqrt of total PS(signal):'
    print '        1/N*sqrt(sum(PS(signal)))\t= ', pspower
    print '    sqrt of variance, rms fluct level, signal:'
    print '        sqrt(variance(signal))\t\t= ', varpower


shot = 1101209029
tree = MDSplus.Tree('cmod', shot)
time_of_shot = tree.getNode('time_of_shot').data()
print time_of_shot
node = tree.getNode('spectroscopy.halph_dalph.analysis.h_to_d_ratio')
(time, signal) = (node.dim_of().data(), node.data())

PS = np.abs(np.fft.fft(signal))**2

# Get frequency bin centers in cycles per unit of sample spacing
freqs = np.fft.fftfreq(signal.size, 1./len(time))    # assuming 1s shot
# Argument order corresponding to ascending frequencies
idx = np.argsort(freqs)

plt.figure()

plt.subplot(211)
plt.plot(time, signal)
plt.xlabel('Time (s)')
plt.ylabel('Signal (V)')
plt.subplot(212)
plt.plot(freqs[idx], PS[idx])
plt.xlim([0, 10])
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power spectrum')

plt.show()
