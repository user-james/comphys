import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':

    data = np.genfromtxt("gruneisen.txt")
    
    scaling = data[:, 0]
    pressure = -1*data[:, 1]
    grun = data[:, 2]
    temperature = pressure/grun

    slope = (scaling[36]- scaling[31])/(temperature[36] - temperature[31])
    print("Thermal Expansivity =", slope)

    plt.plot(temperature, scaling)    
    plt.plot(temperature[31], scaling[31], 'r.')
    plt.plot(temperature[36], scaling[36], 'r.')
    plt.xlabel("Temperature")
    plt.ylabel("Scale Factor")
    plt.show()
