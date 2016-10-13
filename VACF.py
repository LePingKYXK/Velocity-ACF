#!/usr/bin/env python


from __future__ import division
from scipy import signal
import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import time

fname = sys.argv[1]
delta_t = float(sys.argv[2]) * 1.0e-15
window = sys.argv[3]
fout = sys.argv[4]



#### The functions will be used are listed as follows

def select_atoms():
    sel = raw_input("\n (NOTE: The input values are the IDs of the molecule," 
                    "    which should in range from 0 to the max N of the system.\n"
                    "    The total number of your inputs should in range 2 to 4 (including boundary).\n"
                "Here, 2 digitals (IDs) represent the diatomic stretch, \n"
                "      3 digitals (IDs) mean the bend mode (triatomic angle), \n"
                "      4 digitals (IDs) indicate the dihedral angle,\n"
                "\nPlease enter the IDs: \n")

    fmt = '\n The indices of the atoms you chose are: ' + '{:s} ' * len(sel)
    print fmt.format(*sel)
    if sel[0] == '-1':
        print "Program will deal with all atoms."
        return 
    elif 2 <= len(sel) <= 4:
        print "Program will deal with the selected atoms according to your inputs."
        return map(int,sel.split())
    else:
        print "Error! Please check the inputs."


def information(fname):
    ''' fetch the munber of atoms and the elements information
    from the first block of the input file. It looks duplicate,
    but it would be more convenient for the following steps.
    '''
    elements = []
    with open(fname,'r') as fo:
        Natom = int(fo.next())
        fo.next()
        for i in xrange(Natom):
            line = fo.next()
            info = line.split()[0]
            elements.append(info)
    return (Natom, elements)


def screen_print(Natom,elements):
    species = set(elements)
    fmt1 = 'There are {:d} kinds of elements, which are '\
            + '{:s}, ' * len(species)
    print "\nNumber of atoms in this system: ", Natom
    print fmt1.format(len(species), *species)

    elements = np.asarray(elements)
    indices = np.where(elements)[0]
    if len(elements) <= 30:
        print "\nThe elements and the corresponding indices are:"
        fmt2 = 'Elements:' + '{:>3s}, ' * len(elements)
        fmt3 = 'Indices: ' + '{:>3d}, ' * len(indices)
        print fmt2.format(*elements)
        print fmt3.format(*indices)

        
def choose_mode():
    mode = raw_input("\n Please choose the mode: \n"
                     "step: read data step by step, the finial data is a 2D array. \n"
                     "      suit for large file (larger than 1 GB)\n" 
                     "once: read data at once, the finial data is a 3D array. \n"
                     "      suit for smal file (less than 1 GB)\n")
    return mode


def read_3d_data(fname, num_atoms):
    ''' read data at once from input file.
    The final data is a 3D array.
    '''
    time_step = 0
    coords = []
    with open(fname, 'r') as fo:
        for line in fo:
            try:
                fo.next()
            except StopIteration:
                break
            for n in xrange(num_atoms):
                line = fo.next()
                info = line.split()[1:]
                coords.append(info)
            time_step += 1
    coords np.array(coords, np.float64).reshape(time_step, num_atoms, 3)
    
    if 2 <= len(sel) <= 4:
        # read data according to the order of the inputs.
        return coords[:,np.array(sel),:]


def calc_2D_derivative(array_2D, delta_t):
    dy = np.zeros(np.shape(array_2D))
    dx = np.repeat(delta_t, len(array_2D[:,0]))
    for i in xrange(3):
        dy[:,i] = np.gradient(array_2D[:,i], dx, edge_order=2)
    return dy


def calc_bond(data):
    return np.linalg.norm((data[:,0,:] - data[:,1,:]), axis=1)


def calc_angle(data):

    v1 = data[:,1,:] - data[:,0,:]
    v2 = data[:,2,:] - data[:,0,:]
    
    dot = (v1 * v2).sum(axis=1)
# calculate the dot product using Einstein summation
#    dot = np.einsum("ij,ij->i", v1,v2)

    norm1 = np.linalg.norm(v1,axis=1)
    norm2 = np.linalg.norm(v2,axis=1)
    
    theta = np.arccos(dot/(norm1 * norm2))
    
    return np.degrees(theta)


def calc_dihedral(data):

    v1 = data[:,1,:] - data[:,0,:]
    v2 = data[:,2,:] - data[:,1,:]
    v3 = data[:,3,:] - data[:,2,:]

    n1 = np.cross(v1,v2)
    n2 = np.cross(v2,v3)

    dot = (n1 * n2).sum(axis=1)

    norm1 = np.linalg.norm(n1,axis=1)
    norm2 = np.linalg.norm(n2,axis=1)
    
    phi = np.arccos(dot/(norm1 * norm2))
    
    return np.degrees(phi)


def choose_window(data, kind=window):
    if kind == "Gaussian":
        sigma = 2 * math.sqrt(2 * math.log(2))
        window = signal.gaussian(len(data), std=4000/sigma, sym=False)
    elif kind == "BH":
        window = signal.blackmanharris(len(data), sym=False)
    elif kind == "Hamming":
        window = signal.hamming(len(data), sym=False)
    elif kind == "Hann":
        window = signal.hann(len(data), sym=False)
    return window


def zero_padding(sample_data):
    ''' A series of Zeros will be padded to the end of the dipole moment
    array (before FFT performed), in order to obtain a array with the
    length which is the "next power of two" of numbers.
    
    #### Next power of two is calculated as: 2**math.ceil(math.log(x,2))
    #### or Nfft = 2**int(math.log(len(data_array)*2-1, 2))
    '''
    return int(2 ** math.ceil(math.log(len(sample_data), 2)))


def calc_2D_ACF(array_2D):
    yunbiased = array_2D - np.mean(array_2D, axis=0)
    ynorm = np.sum(np.power(yunbiased,2), axis=0)
    print "the average value of input data array", ynorm
    autocor = np.zeros(np.shape(array_2D))
    for i in xrange(3):
        autocor[:,i] = signal.fftconvolve(array_2D[:,i],
                                          array_2D[:,i][::-1],
                                          mode='full')[len(array_2D)-1:] / ynorm[i]
    print "check point 03: shape of the signal.FFTcorrelate()", np.shape(autocor)
    return autocor


def VACF_total(data):
    for i in xrange(len(data)):
        autocor_i = calc_ACF(data[i])
        if i == 0:
            accum_autocor = autocor_i
        else:
            accum_autocor += autocor_i
    return accum_autocor


def calc_FFT_array_2D(array_2D, window):
    '''
    This function is for calculating the "intensity" of the ACF at each frequency
    by using the discrete fast Fourier transform.
    '''
####
#### http://stackoverflow.com/questions/20165193/fft-normalization
####
    window = choose_window(array_2D, "Gaussian")
    WE = sum(window) / len(array_2D)
    wf = window / WE
    # convolve the blackman-harris window function.
    sig = array_2D * wf[None,:].T

    # A series of number of zeros will be padded to the end of the DACF array before FFT.
    N = zero_padding(sig)

    yfft = np.fft.fft(sig, N, axis=0) / len(sig)
#    yfft = np.fft.fft(segment, n=int(N_fft), axis=0)/len(segment)    # without window function
    
    return np.square(np.absolute(yfft))


def save_results(fout, wavenumber, intensity):
    with open(fout, "w") as fw:
        title = ("Wavenumber", "IR Intensity", "cm^-1", "a.u.")
        np.savetxt(fout, np.c_[wavenumber[0:4000],intensity[0:4000]],
                   fmt="%10.5f %15.5e",
                   header="{0:>10}{1:>16}\n{2:^11}{3:^20}".format(*title),
                   comments='')


def visualization(derivative, ACF, wavenumber, intensity):
    matplotlib.style.use('ggplot')
    derivative = derivative * delta_t
    plt.subplot(3,1,1)
    L1 = np.arange(len(derivative))
    plt.plot(L1, derivative, color='red', linewidth=1.5)
    plt.axis([0, len(derivative), 
              1.1*np.min(derivative), 1.1*np.max(derivative)], fontsize=15)
    plt.xlabel("Data Points", fontsize=15)
    plt.ylabel("Derivative of Dipole (a.u.)", fontsize=15)

    plt.subplot(3,1,2)
    L2 = np.arange(len(ACF))
    plt.plot(L2, ACF, color='red', linewidth=1.5)
    plt.axis([0, len(ACF), 1.1*np.min(ACF), 1.1*np.max(ACF)], fontsize=15)
    plt.xlabel("Data Points", fontsize=15)
    plt.ylabel("DACF (a.u.)", fontsize=15)

    plt.subplot(3,1,3)
    plt.plot(wavenumber, intensity, color='black', linewidth=1.5)
    plt.axis([0, 4000,
             -1.1*np.min(intensity), 1.1*np.max(intensity)],
             fontsize=15)
    plt.xlabel("Wavenumber (cm$^{-1}$)", fontsize=15)
    plt.ylabel("Intensity (a.u.)", fontsize=15)
    plt.subplots_adjust(hspace = 0.5)
    plt.show()



########################
delta_t = 0.5 * 1e-15
#Fs = 1/delta_t
T = 50.0
#scaling_factor = 0.968

######## The constants will be used in this script ########
c = 2.9979245899e10 # speed of light in vacuum in [cm/s], from Wikipedia.
kB = 0.6950347      # Boltzman constant in [cm^-1/K], from Wikipedia.
h_bar = 6.283185    # Reduced Planck constant in atomic unit, where h == 2*pi
beta = 1.0/(kB * T) #                        




######## The main program ########
if __name__ == "__main__":
    Natom, elements = information(fname)
    screen_print(Natom,elements)
    sel = select_atoms()

    start = time.clock()
    
    data = read_data(fname, sel)
#    print "check point 01: data", data, np.shape(data)
    # data is a 3-D array contained coordinates of selected atoms in whole trajectory.
    
    if len(sel) == 1 and sel[0] == -1:
        print "Program will deal with all atoms."
        for i in range(len(data[0,:,:])):
            print '\n %d===' %i
            derivative = calc_2D_derivative(data[:,i,:], delta_t)
            VACF = calc_2D_ACF(derivative)
            yfft_i = calc_FFT_array_2D(VACF, window)
            if i == 0:
                yfft = yfft_i
            else:
                yfft += yfft_i
        print "\ncheck point 04: yfft = \n", yfft, np.shape(yfft)
        
    elif len(sel) == 2:
        values = calc_bond(data)
    elif len(sel) == 3:
        values = calc_angle(data)
    elif len(sel) == 4:
        values = calc_dihedral(data)
    else:
        print "Error!"
    derivative = calc_derivative(values, delta_t)
    print "derivative of the data =", np.shape(derivative)

    ACF = calc_ACF(derivative)
    yfft = calc_FFT(ACF)
    
    wavenumber = np.fft.fftfreq(len(yfft), delta_t * c)[0:int(len(yfft) / 2)]
    intensity = np.sum(yfft, axis=1)[0:int(len(yfft)/2)]

#### Normalized the intensity
#    norm_intensity = intensity / max(intensity)
    save_results(fout,wavenumber,intensity)
    finish = time.clock()
    print "Work Completed! Used time: %.5f second." %(finish - start)
    visualization(derivative, ACF, wavenumber, intensity)
    

    
'''
def read_2d_data(fname, Natom, sel):
    ''' Using iterator to read the file to a list line by line,
    then convert this list into ndarray.
    The way of conversion dependends one the input numbers.
    The order of the inputs also took into consideration.
    '''
    coords = []
    for n in xrange(Natom):
        line = fo.next()
        xyz = line.split()[1:]
        coords.append(xyz)
    coords = np.asfarray(coords,dtype=np.float64).reshape(timestep,Natom,3)

    if 2 <= len(sel) <= 4:
        # read data according to the order of the inputs.
        coords = coords[:,np.array(sel),:]
    return coords
'''
