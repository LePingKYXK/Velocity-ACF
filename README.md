# Velocity-ACF
This Python (version 2.7) Script is one of my projects dealing with the Spectrum. It Based on Fast Fourier Transform (FFT) of the Velocity Autocorrelation Function (VACF).
    
This script first read the position file (with the extension of .xyz) containing the cartesian coordinates, whcih is generated from the CP2K/QuickStep simulations. And then calculate the time derivative of each component (x, y, z), yielding the velocities of each compoent (v_i, here i = x, y, z). After computing the autocorrelation of each velocity component, v_i, the VACF data array obtained. By performing the FFT on the VACF, the final spectrum produced. And then plotted on the graph panel by using the Matplotlib module.

The spectrum calculated from VACF would different from the IR spectrum obtained by DACF, since the IR spectrum obeys the selection rules which is intrinsic represented by total dipole moment data. On the contrary, the spectrum yielded from VACF does not obey the selection rules. As a result, the IR forbidden transitions would appear in VACF spectrum.

Modules required:

- Numpy (version 1.9.1 or above)
- Scipy (version 0.17.0 or above)
- Matplotlib (version 1.4 or above)

*Note:*
If the whole molecules or clusters shift or rotate when displaying the trajectory in VMD, it is strongly recommended to align the trajectory to the first frame using VMD before running this script. Using the aligned position file as the input file of the VACF script will get better result. The reason for performing the alignment is to reduce the influent bringing from the translation and rotation of the whole system to the finial VACF spectrum.
