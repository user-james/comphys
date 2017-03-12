from subprocess import call
import numpy as np
import matplotlib as plt



if __name__ == '__main__':

    v1 = [1.75, 2, 2.25, 2.5]
    v2 = [3.75, 4, 4.25, 4.5]

    data = []
        
    for i in v1:
        for j in v2:
            num1, num2 = str(i), str(j)
            call(["../LJ_Quantum_Dimer/NFE_wave_functions/NFE_wavefunctions", num1, num2])
            
            bands = np.genfromtxt("NFE_bands")

            band1 = min(bands[:, 2]) - max(bands[:, 1])
            band2 = min(bands[:, 3]) - max(bands[:, 2])

            err1 = abs(i - band1)
            err2 = abs(j - band2)

            data.append([v1, err1, v2, err1])

#    np.vstack(data)
#    np.savetxt("band_gap_comp", data, delimiter=',')
    print(data)


    
