import numpy as np
import sympy as sym
import matplotlib.pyplot as plt

r = sym.symbols('r')

if __name__=='__main__':
    
    # Potential definition
    V = 1/r**12 - 1/r**6

    # calculates omega and d(omega)/dr symbolically
    omega_sq = sym.Derivative(sym.Derivative(V, r).doit()).doit()
    omega = sym.sqrt(omega_sq)    
    d_omega = sym.Derivative(omega, r).doit()
    
    # transforms symbolic function to numerical function so it can be plotted
    log_deriv = sym.lambdify(r, d_omega/omega, "numpy")
    N_omega = sym.lambdify(r, omega, "numpy")

    print("Freq. for LJ-Pot. @ Equilibrium: ", N_omega(1.122))


                ############
                # PLotting #
                ############
    

    data = np.genfromtxt('equil_mode_frequencies')
   
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx().twiny()
    

    r = np.linspace(0, 20, 100)
    ax2.plot(r, log_deriv(r), 'r', label="LJ-Potential")
    ax1.plot(data[:, 1], data[:, 2], label="Sim. Data")
    ax2.set_ylim(-40, 5)
    
    ax1.set_xlabel("Frequency")
    ax2.set_xlabel("Atomic Separation")

    ax1.set_ylabel("Log Derivative")
    ax2.set_ylabel("Log Derivative (LJ Potential)")

    ax1.legend(loc=3)
    ax2.legend(loc=4)
    plt.savefig("l_deriv_comp.pdf")
    plt.show()
