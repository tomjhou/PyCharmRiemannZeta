import numpy as np
# import matplotlib
import matplotlib.pyplot as plt
from scipy.special import gamma, loggamma
import platform
import time

import riemann_math as rm
# matplotlib.use("TkAgg")

if __name__ == "__main__":
    print("Done importing")

height = 4000
offset = 5
gap = 0.1

if rm.RIEMANN_ITER_LIMIT < height:
    rm.RIEMANN_ITER_LIMIT = height

rm.precompute_coeffs()
print("Done computing coefficients, iteration limit = " + str(rm.RIEMANN_ITER_LIMIT))

# Input vector along critical line Re[s] = 0.5
imag_part = np.arange(1, height, gap) * 1j
s = imag_part + 0.5
ax = np.imag(s)  # this is what will be used for x-axis of plot

# Input vector along line Re[s] = 5
s2 = imag_part + offset

# Make figure now so user doesn't have to stare at blank screen too much longer
plt.figure()
# Reduce margins
plt.tight_layout()
# Make even smaller margins
plt.subplots_adjust(left=0.05, right=0.99, top=0.95, bottom=0.05)

# plt.subplot(2,1,1)
plt.axhline(color='k')
plt.title('Normalized to gamma magnitude')


# plt.show() will pause the GUI, so we use plt.draw() followed by plt.pause() instead.
def show():
    plt.draw()
    plt.pause(0.001)


# This will make a blank window show up so user knows something is about to happen
show()

#
# Calculate riemann values at line Re[s] = 0.5
#
# We want to calculate this:
#     y = -np.real((np.pi ** (offset / 2)) * rm.riemann_symmetric(s) / abs(gamma(s / 2)))
# but because gamma(s) rounds to 0 for very large s, we calculate logs instead. Note that
# log(abs(gamma(s)) = Re[log(gamma(s))], as given by derivation below:
#
# For any complex c = a+bi
#     log(c) = log(abs(c)) + i*Arg(c)
# Hence, for any complex c,
#     log(abs(c)) = Re(log(c))

t0 = time.time()
y_log = rm.riemann_symmetric(s, use_log=True) - np.real(loggamma(s / 2))
y = np.real(np.exp(y_log) / (ax * ax + 0.25))
plt.plot(ax, y, linewidth=1)
delay = time.time() - t0
print("Plotted values along critical line Re[s]=0.5 in %1.2f seconds" % delay)
show()

# Calculate at line Re[s] = offset
y2_log = np.log(np.pi) * (offset / 2) + rm.riemann_symmetric(s2, use_log=True) - np.real(loggamma(s2 / 2))
y2 = np.real(np.exp(y2_log) / (ax * ax + offset * offset))
plt.plot(ax, -y2, linewidth=1)
delay = time.time() - t0 - delay
print("Plotted values along Re[s] = %1.1f in %1.2f seconds" % (offset, delay))
show()

if __name__ == "__main__":
    # This keeps main window open until dismissed. Don't call this if we came here from elsewhere.
    print("Close window to exit program")
    plt.show()
