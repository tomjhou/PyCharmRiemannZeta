import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import riemann_math as rm
from scipy.special import gamma, loggamma

# matplotlib.use("TkAgg")

print("Done importing")

if rm.RIEMANN_ITER_LIMIT < 1000:
    rm.RIEMANN_ITER_LIMIT = 1000

rm.precompute_coeffs()
print("Done computing coefficients, limit = " + str(rm.RIEMANN_ITER_LIMIT))

height = 2000
offset = 5
num_points = height * 10  # density of points

# Input vector along critical line Re[s] = 0.5
s = np.linspace(0j, 1j * height, num_points) + 0.5
ax = np.imag(s)  # this is what will be used for x-axis of plot

# Input vector along line Re[s] = 5
s2 = np.linspace(0j, 1j * height, num_points) + offset

# Make figure now so user doesn't have to stare at blank screen too much longer
plt.figure()
#plt.subplot(2,1,1)
plt.axhline(color='k')
plt.title('Normalized to gamma magnitude')


# plt.show() creates a second mainloop() that causes everything to hang if there is
# already a mainloop running. Use plt.draw() followed by plt.pause() instead.
def show():
    plt.draw()
    plt.pause(0.001)


#
# Calculate riemann values at line Re[s] = 0.5
# Each plot may take a while.
#
# For any complex c = a+bi
#     log(c) = log(abs(c)) + i*Arg(c)
# Hence, for any complex c,
#     log(abs(c)) = Re(log(c))

# We want to calculate this:
#     y = -np.real((np.pi ** (offset / 2)) * rm.riemann_symmetric(s) / abs(gamma(s / 2)))
# but because gamma(s) rounds to 0 for very large s, we calculate logs instead. Note that
# log(abs(gamma(s)) = Re[log(gamma(s))]
#
y_log = rm.riemann_symmetric(s, use_log=True) - np.real(loggamma(s/2))
y = np.real(np.exp(y_log))
plt.plot(ax, y, linewidth=1)
show()

# Calculate at line Re[s] = offset
y2_log = np.log(np.pi) * (offset / 2) + rm.riemann_symmetric(s2, use_log=True) - np.real(loggamma(s2 / 2))
y2 = np.real(np.exp(y2_log))
plt.plot(ax, -y2, linewidth=1)
show()



if __name__ == "__main__":
    # This keeps main window open until dismissed. Don't call this if we came here from elsewhere.
    plt.show()
