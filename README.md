# Velocity-ACF
    The Python Script is one of my projects dealing with the Spectrum. It Based on Fast Fourier Transform (FFT) of the Velocity Autocorrelation Function (VACF).
    
This script first read the position file containing the cartesian coordinates with the extension of .xyz, whcih is generated from the CP2K/QuickStep simulations. And then calculate the time derivative of the Dipole, yielding the dipole prime (D_p in short). After computing the autocorrelation of the D_p, the DACF data array obtained. By performing the FFT on the DACF, the final IR spectrum produced. And then plotted on the graph panel by using the Matplotlib module.

Modules required:

- Numpy (version 1.9.1 or above)
- Scipy (version 0.17.0 or above)
- Matplotlib (version 1.4 or above)
