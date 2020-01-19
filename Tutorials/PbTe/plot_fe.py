from numpy import *
from matplotlib.pyplot import *


fig, axarr= subplots(ncols = 3, nrows = 1, dpi = 100, figsize = (10,4),
                     sharex = True)

x = linspace(-1, 1, 1000)
f_m = 2*x**4 - x**2
f_c = x**4
f_M = 0.5*(x**2 + x**4)
axarr[0].set_title(r"$ T < T_c$")
axarr[0].plot(x, f_m)
axarr[0].set_xlabel("Q [order parameter]")
axarr[0].set_ylabel("Free energy")

axarr[1].set_title(r"$T = T_c$")
axarr[1].plot(x, f_c)
axarr[1].set_xlabel("Q [order parameter]")

axarr[2].set_title(r"$T > T_c$")
axarr[2].plot(x, f_M)
axarr[2].set_xlabel("Q [order parameter]")


fig.tight_layout()
fig.savefig("second_order.png")
show()
