import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from SMCPy.Grid import NC2SMC

smc = NC2SMC('./test_global.nl')
smc.run()
smc.plot_patches(proj=ccrs.PlateCarree(central_longitude=180))
plt.show()
