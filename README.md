# Velocity-ACF
This Python (version 2.7) Script is one of my projects dealing with the Spectrum. It based on Fast Fourier Transform (FFT) of the Velocity Autocorrelation Function (VACF).  
    
This script first read the position file (with the extension of .xyz) containing the cartesian coordinates, whcih is generated from the CP2K/QuickStep simulations. And then calculate the time derivative of each component (x, y, z), yielding the velocities of each compoent (v_i, here i = x, y, z). After computing the autocorrelation of each velocity component, v_i, the VACF data array obtained. By performing the FFT on the VACF, the final spectrum produced. And then plotted on the graph panel by using the Matplotlib module.  

The spectrum calculated from VACF would different from the IR spectrum obtained by DACF, since the IR spectrum obeys the selection rules which is intrinsic represented by total dipole moment data. On the contrary, the spectrum yielded from VACF does not obey the selection rules. As a result, the IR forbidden transitions would appear in VACF spectrum.  

## Modules required:
- Numpy (version 1.9.1 or above)
- Scipy (version 0.17.0 or above)
- Matplotlib (version 1.4 or above)

**Note:**
If the whole molecules or clusters shift or rotate when displaying the trajectory in VMD, it is strongly recommended to align the trajectory to the first frame using VMD before running this script. Using the aligned position file as the input file of the VACF script will get better result. The reason for performing the alignment is to reduce the influent bringing from the translation and rotation of the whole system to the finial VACF spectrum.  

## E-mail address for contacting the authors:
`huan.wang@mail.huji.ac.il` or `wanghuan@iccas.ac.cn (China)`

## PLEASE READ THE FOLLOWING INSTRUCTIONS BEFORE RUNNING SCRIPT 
### The Format for Running This Script:  
`python` `VACF_KW.py` *`DIRECTORY_OF_YOUR_DATA`* *`INPUT_FILE_NAME`* *`DELTA_T`* *`OUTPUT_FILE_NAME`*  

### The values need to input manually when runing this script    
  1. `INPUT_FILE_NAME`: The POSITION.xyz file, which contains the coordinates of the system at each time step (frame).  
        + <sub>NOTE: do NOT need to re-split the Position file)</sub>
  2. `DELTA_T`: The Time_step set in simulation, in unit of fs
  3. `OUTPUT_FILE_NAME`: The Name of the Output File.
        + <sub>NOTE: do NOT need to type ">" sign!</sub>

### Main improvment in this version is that the script became more user-friendly
  After inputing the above mentioned values, the program will list the atoms and their corresponding indices (only the first 35 atoms will show if the system is too large). Then the program will ask the user to enter the type of mode, e.g. "s": stretch, "b": bend, "r": internal rotation, and "u": umbrella, and the indices of the group of atoms corresponding to the mode (That is the reason for the script called *partial* VACF).
- If you want to calculate the VACF of stretch mode between two atoms (e.g. the atom No. 5 and atom No. 7), please type  
    `s 5 7`

- If you want to calculate the VACF of bend mode between three atoms (e.g. the angle between atom No. 1, atom No. 2, atom No. 3), please type
    `b 1 2 3`

- If you want to calculate the VACF of internal rotation of within four atoms (e.g. the dihedral angle involve atom No. 8, atom No. 9 atom No. 10, and atom No. 11), please type 
    `r 8 9 10 11`

- If the user enter "all" or "-1", the program will choose all atoms to calcluate the total VACF.  

Consequently, the program will ask the user to enter the type of window function, e.g. "BlackmanHarris", "Gaussian", "Hamming "or "Hann", ect.  

After all steps of inputs finishing, the calculation begins.  

The result will save as a `.txt` file and the plot will be shown on the screen.

### ############################  Let's Try It! ############################ ###


A Python 3 verion of VACF was released on Mar. 18, 2018. 
