Total power test
----------------

File: `exercises/fourier_ps.py`

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

![](http://i.imgur.com/SWlodCA.png)


Cross-power spectrum
--------------------

File: `exercises/cross_power.py`

Cross-power spectrum of sine waves with frequency 3 and 7 Hz:

![](http://i.imgur.com/cXADgdq.png)


FFT autocorrelation 
-------------------

File: `exercises/fft_a_corr.py`

![](http://i.imgur.com/FmDXkP1.png)

i.e. `autocorrelation = inverse_fft(fft*conjugate(fft))`

Taking the real part:

![](http://i.imgur.com/RsMOz27.png)

Faster than direct autocorrelation and seems accurate, although not normalized.


Direct autocorrelation 
----------------------

File: `exercises/a_corr.py`

Autocorrelation computed directly. Output for `sin(2*pi*10*x) + noise`:

![](http://i.imgur.com/8aQptHa.png)

It runs slowly compared to FFT when there are many points in the signal. Also not sure why there is attenuation around the edges.


Power spectrum
--------------

Power spectrum for `sin(2*pi*10.5*x) + sin(2*pi*75.5*x) + noise`. We see wide peaks around 10.5 and 75.5 as expected because the frequency bins are centered at the integers.

![](http://i.imgur.com/v1fLOwE.png)

