import numpy as np
import matplotlib.pyplot as plt


if __name__=='__main__':

    # create array from text file
    evals = np.genfromtxt("e_values.txt")
    
    # neglecting first 3 values as they are negative and non-zero
    # otherwise we would have imaginary frequencies
    freq = np.sqrt(evals[3:])/(2*np.pi)

    # create histogram with 40 bins
    plt.hist(freq, bins=40)
    
    # set x limits to be the same as those for the power spectrum
    plt.xlim(0, 2.5)

    # save and plot
    plt.savefig("eval_spectrum.jpg")
    plt.show()
