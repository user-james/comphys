import numpy as np
import matplotlib.pyplot as plt


if __name__=='__main__':

    waves = []

    for i in range(6):
        waves.append(np.genfromtxt("H_wavefn_{:d}".format(i)))

    #wavenum = int(input("Enter wavenumber (0-5):"))

    #while (wavenum > 5 or wavenum < 0):
    #    print("Wave number not within 0-5 range")
    #    wavenum = int(input("Enter wavenumber (0-5):"))

    #wavefn = waves[wavenum]

    count = 0
    for wavefn in waves:
        plt.plot(wavefn[:, 0], wavefn[:, 1]**2, label="$n={:d}$".format(count))
        count += 1
    
    plt.xlabel("Radius")
    plt.ylabel("Probability")
    plt.legend(loc=1)
    plt.show()
