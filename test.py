import numpy as np
import matplotlib.pyplot as plt
x = np.loadtxt('./build/Release/tdist.dat').reshape((200,200,1))
fig, ax = plt.subplots()
cb = ax.contourf(x[:,:,0], levels=20)
fig.colorbar(cb)
plt.show()