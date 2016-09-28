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
def read_data(fname, sel):
    with open(fname, 'r') as fo:
        line_count = 0
        block_count = 0
        atom_num = int(fo.readline())
        line_per_block = atom_num + 2
        coords = []

        sel_atom_num = atom_num
        #read the rest of the 1st block.
        fo.next()
        for i in xrange(sel_atom_num):
            line = fo.next()
            xyz = line.split()[1:]
            coords.append(xyz)
        block_count += 1
            
        #read the rest blocks.
        for line in fo:
            if line_count > 1:
                xyz = line.split()[1:]
                coords.append(xyz)
            line_count += 1
            if line_count == line_per_block:    #already read a block
                line_count = 0
                block_count += 1

        shape = (block_count, sel_atom_num, 3)
    return np.asfarray(coords,dtype=np.float64).reshape(shape)


def calc_2D_derivative(array_2D, delta_t):
    dy = np.zeros(np.shape(array_2D))
    dx = np.repeat(delta_t, len(array_2D[:,0]))
    for i in xrange(3):
        dy[:,i] = np.gradient(array_2D[:,i], dx, edge_order=2)
    return dy


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


def plot(wavenumber, intensity):
    matplotlib.style.use('ggplot')
    plt.plot(wavenumber,intensity, color='black', linewidth=1.5)
    plt.axis([0,4000, np.min(intensity),1.1*np.max(intensity)], fontsize=15)
    plt.xlabel("Wavenumber (cm$^{-1}$)", fontsize=15)
    plt.ylabel("Intensity (a.u.)", fontsize=15)
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
    # Test. Input 3 numbers in range from 0 to 8.
    sel = raw_input("(NOTE: The input values should in range of the cluster's IDs.\n"
                    "The total number of your inputs should in range of 2 to 4 (including).\n"
                    "Here, 2 digits means diatomic stretch, 3 digits means triatomic angle bend,\
                     4 digits means dihedral angle,\n"
                    "while -1 will read all atoms in the system.)\n"
                    "\nPlease enter numbers: \n")
    sel = map(int, sel.split())
    print "Your inputs are:", sel

    start = time.clock()
    
    data = read_data(fname, sel)
    print "check point 01: data", data, np.shape(data)
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
        
    wavenumber = np.fft.fftfreq(len(yfft), delta_t * c)[0:int(len(yfft) / 2)]
    intensity = np.sum(yfft, axis=1)[0:int(len(yfft)/2)]

#### Normalized the intensity
#    norm_intensity = intensity / max(intensity)
    save_results(fout,wavenumber,intensity)
    print "Work Completed! Used time: %.5f second." %(time.clock() - start)
    plot(wavenumber, intensity)
