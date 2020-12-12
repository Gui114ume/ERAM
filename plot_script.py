import matplotlib.pyplot as plt
import numpy as np

#x,y = np.loadtxt('Data/Perf_ERAM_750_750',  unpack = True)
#sum_x = sum(x)
#sum_y = sum(y)
#moy_exec = sum_y / sum_x
#
#print(moy_exec)

x,y = np.loadtxt("Data/ERAM_50_50_Convergence.txt", unpack = True)
#f.write(str(moy_exec)+"EOF")

plt.plot(x, y)
plt.title(label = 'Convergence matrice 100x100 ERAM')
plt.xlabel('Iteration')
plt.ylabel('Residu')
plt.show()