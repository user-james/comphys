import numpy as np
import matplotlib.pyplot as plt

if __name__=='__main__':
    
    count = 0

    while count < 4:
        num = input("Enter spectrum file number: ")
        filename = 'spectrum_' + num
        data = np.genfromtxt(filename)
        start_pt = int(len(data[:, 0])*0.05/2.5)
        plt.plot(data[start_pt:, 0], data[start_pt:, 1], label='Spectrum {:s}'.format(num))
        
        count += 1
    
    plt.legend(fontsize=11)
    plt.savefig("power_spectrum.jpg")
    plt.show()
