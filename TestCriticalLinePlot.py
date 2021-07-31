import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import riemann_math as rm
from scipy.special import gamma

# matplotlib.use("TkAgg")

print("Done importing")

rm.RIEMANN_ITER_LIMIT = 200
rm.precompute_coeffs()
print("Done computing coefficients")

height = 200
offset = 5
num_points = height * 10  # density of points

# Input vector along critical line Re[s] = 0.5
s = np.linspace(0j, 1j * height, num_points) + 0.5
ax = np.imag(s)  # this is what will be used for x-axis of plot

# Input vector along line Re[s] = 5
s2 = np.linspace(0j, 1j * height, num_points) + offset

# Make figure now so user doesn't have to stare at blank screen too much longer
plt.figure()
plt.subplot(2,1,1)
plt.axhline(color='k')
plt.title('Normalized to gamma magnitude')

# plt.show() creates a second mainloop() that causes everything to hang if there is
# already a mainloop running. Use plt.draw() followed by plt.pause() instead.
def show():
    plt.draw()
    plt.pause(0.001)

# Each plot may take a while.
y = np.real(rm.riemann_symmetric(s) / abs(gamma(s / 2)))
plt.plot(ax, y)
show()

y2 = -np.real((np.pi ** (offset / 2)) * rm.riemann_symmetric(s2) / abs(gamma(s2 / 2)))
plt.plot(ax, y2)
show()

plt.subplot(2,1,2)
plt.title('Normalized to gamma magnitude+phase')
y = np.real(rm.riemann_symmetric(s) / gamma(np.conj(s) / 2))
plt.plot(ax, y)
show()

y2 = np.real((np.pi ** (offset / 2)) * rm.riemann_symmetric(s2) / gamma(s2 / 2))
plt.plot(ax, y2)
plt.axhline(color='k')
show()

if __name__ == "__main__":
    # This keeps main window open until dismissed. Don't call this if we came here from elsewhere.
    plt.show()
