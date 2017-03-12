import numpy as np
import matplotlib.pyplot as plt

def avg_err(n):
    '''
    Calculates average error across all energy eigenvalues in a file
    given the number of grids used to calculate the eigenvalues
    '''

    # specifies filename and creates data set from it
    filename = "e_{:d}points".format(n)
    data = np.genfromtxt(filename)
    err = 0

    # length of each column in file
    size = len(data[:, 0])

    # calculates sum of errors
    for n in range(size):
        err += abs((n+1)**2 * np.pi**2 - data[n, 1])

    # returns average error
    return err/size

def expected_curve(x):
    return 1/x**2

def err_curve(x):
    return 1/x

if __name__=='__main__':

    # creates range from 100-1000 (including 1000) in steps of 50 
    grids = np.arange(100, 1001, 50) 
    errors = []

    # x space for a error plot comparison
    x = np.linspace(100, 1000, 1000)

    # calculates error for each number of grid points we tested
    for num in grids:
        errors.append(avg_err(num))

    # plots curve for simulated error
    plt.subplot(311)
    plt.plot(grids, errors)
    plt.ylim(.25, 3)
    plt.ylabel("Average Error")
    plt.xlabel("# Grid Points")
    plt.title("Error Comparison")

    # plots curve that we expect to look similar to our error
    plt.subplot(312)
    plt.plot(x, expected_curve(x), 'r-')
    plt.ylabel(r'$\frac{1}{x^2}$')
    plt.xlabel(r'$x$')

    # plots curve that we expect to look similar to our error
    plt.subplot(313)
    plt.plot(x, err_curve(x), 'k-')
    plt.ylabel(r'$\frac{1}{x}$')
    plt.xlabel(r'$x$')
    

    # format, save, plot
    plt.tight_layout()
    plt.savefig("error_comp.jpg")
    plt.show()
