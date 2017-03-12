import numpy as np
import matplotlib.pyplot as plt


def wave_n(n):
    ''' Returns analytic funtion for given n '''
    return lambda x: np.sqrt(2)*np.sin(n*np.pi*x)

if __name__=='__main__':

    # written so that this can have a max of 3
    num_wavefns = 3 

    # just used to define subplots
    base_plt = num_wavefns*100 + 20
    subplt = 0

    # define x space to plot analytical function
    x = np.linspace(0,1, 100)

    for n in range(num_wavefns):
        simwave = np.genfromtxt("box_wavefn_{:d}".format(2*n))
         
        # simulation plot
        subplt += 1
        plt.subplot(base_plt + subplt)
        plt.plot(simwave[:, 0], simwave[:, 1])
        plt.title(r'n = {:d}'.format(2*n+1))
        plt.ylabel(r'$\psi_{n}(x)$')
        plt.xlabel(r'$x$')

        # analytic plot
        subplt += 1
        wavefn = wave_n(2*n+1)
        plt.subplot(base_plt + subplt)
        plt.plot(x, wavefn(x))
        plt.title(r'n = {:d}'.format(2*n+1))
        plt.ylabel(r'$\psi_{n}(x)$')
        plt.xlabel(r'$x$')

    
    # format, save and plot
    plt.tight_layout()
    plt.savefig("wave_comp_1000.jpg")
    plt.show()
