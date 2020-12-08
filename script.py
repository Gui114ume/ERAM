import matplotlib.pyplot as plt
import numpy as np

x,y = np.loadtxt('Data/Perf_10x10_seq',  unpack = True)
sum_x = sum(x)
sum_y = sum(y)
moy_exec = sum_y / sum_x

print(moy_exec)

f = open("Data/Perf_10x10_seq_moy_exec", "w+")
f.write(str(moy_exec)+"EOF")

plt.scatter(x, y)

plt.xlabel('Nombre d\'executions')
plt.ylabel('Temps d\'execution')
plt.savefig("Data/plot/fig1.png")