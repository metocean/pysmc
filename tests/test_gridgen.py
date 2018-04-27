import matplotlib.pyplot as plt
import cartopy.crs as ccrs
from SMCPy.Grid import GRIDGEN2SMC

smc = GRIDGEN2SMC('./test_gridgen.nl')
smc.run()
smc.plot_patches(proj=ccrs.PlateCarree(central_longitude=180))
plt.show()
