import matplotlib.pyplot as plt
import numpy as np

# Define the number of points and the interval:
N = 100
T = 0.1

# Midpoint+1 is the most negative frequency subscript:
N21 = N/2 + 1

# The array of subscripts:
F = range(N)
# Insert negative frequencies in elements F(N/2 +1), ..., F(N-1):
# F[N21] = N21 -N + FINDGEN(N21-2)
F[N21:] = [f + N21 - N for f in F[N21:]]
# Compute T0 frequency:
F = np.array([f/(N*T) for f in F])

# Power spectrum is squared magnitude of Fourier transform
ps = np.abs(np.fft.fft(F))**2
# Get frequency bin centers in cycles per unit of sample spacing
freqs = np.fft.fftfreq(F.size, T)
# Argument order corresponding to ascending frequencies
idx = np.argsort(freqs)

plt.figure()
plt.plot(freqs[idx], ps[idx]) # keep in mind this is not a log plot
plt.xlabel('Frequency (Hz)')
plt.ylabel('Power spectrum')
plt.show()

"""
Original code
http://www.physics.nyu.edu/grierlab/idl_html_help/F4.html#wp676828

; Define the number of points and the interval:
N = 100
T = 0.1

; Midpoint+1 is the most negative frequency subscript:
N21 = N/2 + 1

; The array of subscripts:
F = INDGEN(N)
; Insert negative frequencies in elements F(N/2 +1), ..., F(N-1):
F[N21] = N21 -N + FINDGEN(N21-2)

; Compute T0 frequency:
F = F/(N*T)

; Shift so that the most negative frequency is plotted first:
PLOT, /YLOG, SHIFT(F, -N21), SHIFT(ABS(FFT(F, -1)), -N21)

"""
