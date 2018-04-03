import matplotlib.pyplot as plt
from SMCPy.SMCUKMO import NC2SMC

smc = NC2SMC('./test.nl')
smc.run()
plt.show()
