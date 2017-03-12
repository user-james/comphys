import numpy as np
import matplotlib.pyplot as plt
if __name__=='__main__':
    data = np.genfromtxt("equil_mode_frequencies")
   
    widths = (0.05, 0.1, 0.2)

    bin_width = float(input("Enter Bin Width: "))
    mn, mx = min(data[:, 1]), max(data[:, 1])
    nbins = int((mx - mn)/bin_width)

    plt.hist(data[:, 1], bins=nbins)
    plt.title("Histogram of Equilibrium Mode Frequency")
    plt.xlabel("Frequency")

    if bin_width in widths:
        plt.savefig("freq_histogram_{:.2f}.pdf".format(bin_width))

    plt.show()

