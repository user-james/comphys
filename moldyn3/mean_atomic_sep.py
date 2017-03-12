import numpy as np
import matplotlib.pyplot as plt


def weighted_mean(values, weights):
    '''
    This function takes in two 'numpy' arrays corresponding to values to be averaged and
    their respective weights.

    A weighted average is then calculated and returned
    '''
    
    numerator, denominator = 0, 0

    for i in range(len(weights)):
            numerator += values[i]*weights[i]
            denominator += values[i]

    return numerator/denominator


if __name__ == '__main__':

    w_means = []

    for dist_num in range(1, 19):

        #dist_num = input("Enter distribution number interested in: ")
        pair_dist = np.genfromtxt("pair_distribution_" + str(dist_num))

        w_means.append(weighted_mean(pair_dist[:, 0], pair_dist[:, 1]))

    # temperatures were inputted manually from file 'heat_cap_out.txt'
    temperatures = np.array([0.003741, 0.007202, 0.010944, 0.013278, 0.016081, 0.016796, 0.020078, 0.023263, 0.026217, 0.029421, 0.032574, 0.034132, 0.038733, 0.041114, 0.043668, 0.045586, 0.047795, 0.050532])
    w_means = np.array(w_means)

    #####################################################################################################
    # This next section uses a Python package to fit a line to the atomic separetion vs. temperature data
    #####################################################################################################
    temperature_matrix = np.vstack([temperatures, np.ones(len(temperatures))]).T
    slope, constant = np.linalg.lstsq(temperature_matrix, w_means)[0]
    line = slope*temperatures + constant    # line looks like y = m*x + constant

    # thermal expansivity
    therm_exp = slope/w_means[0]
    print("Thermal Expansivity =", therm_exp)

    # Plotting
    plt.plot(temperatures, w_means, '.')
    plt.plot(temperatures, line, 'r', label='Slope: {:.4f}'.format(slope))
    plt.legend(loc = 4)
    plt.xlabel("Temperature")
    plt.ylabel("Average Atomic Distance")
    plt.savefig("separation_vs_temp.jpg")
    plt.show()
