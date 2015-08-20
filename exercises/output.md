Total power test
----------------

File: `examples/fourier_ps.py`

Calculations of total power from the power spectrum and the variance of the signal agree:

    1V amplitude sine wave at 300 Hz:
        sqrt of total PS(signal):
            1/N*sqrt(sum(PS(signal)))   =  0.707106781187
        sqrt of variance, rms fluct level, signal:
            sqrt(variance(signal))      =  0.707106781187
    Normally distributed noise:
        sqrt of total PS(signal):
            1/N*sqrt(sum(PS(signal)))   =  1.0014379912
        sqrt of variance, rms fluct level, signal:
            sqrt(variance(signal))      =  1.00143737174
    1V sine wave + noise:
        sqrt of total PS(signal):
            1/N*sqrt(sum(PS(signal)))   =  1.22616075395
        sqrt of variance, rms fluct level, signal:
            sqrt(variance(signal))      =  1.22616024802

![](resources/fourier_ps.png)


Cross-power spectrum
--------------------

File: `examples/cross_power.py`

Cross-power spectrum of sine waves with frequency 3 and 7 Hz:

![](resources/cross_power.png)


FFT autocorrelation 
-------------------

File: `examples/fft_a_corr.py`

![](http://mathworld.wolfram.com/images/equations/FourierTransform/NumberedEquation3.gif)

i.e. `autocorrelation = inverse_fft(fft*conjugate(fft))`

Taking the real part:

![](resources/fft_a_corr.png)

Faster than direct autocorrelation and seems accurate, although not normalized.


Direct autocorrelation 
----------------------

File: `examples/a_corr.py`

Autocorrelation computed directly. Output for `sin(2*pi*10*x) + noise`:

![](resources/a_corr.png)

It runs slowly compared to FFT when there are many points in the signal. Also not sure why there is attenuation around the edges.


Power spectrum
--------------

Power spectrum for `sin(2*pi*10.5*x) + sin(2*pi*75.5*x) + noise`. We see wide peaks around 10.5 and 75.5 as expected because the frequency bins are centered at the integers.

![](resources/ps.png)

