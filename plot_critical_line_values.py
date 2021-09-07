import numpy as np
# import matplotlib
import matplotlib.pyplot as plt
from scipy.special import gamma, loggamma
import platform
import time

import riemann_math as rm
# matplotlib.use("TkAgg")

def show_critial_line():
    imag_start = 1
    imag_end = 30000
    imag_gap = 0.01
    show_off_critical = False

    real_offset = 5

    iter_limit = imag_end + 2000

    if rm.RIEMANN_ITER_LIMIT < iter_limit:
        rm.RIEMANN_ITER_LIMIT = iter_limit  # int(height * 3)

    rm.precompute_coeffs()
    print("Done computing coefficients, iteration limit = " + str(rm.RIEMANN_ITER_LIMIT))

    # Input vector along critical line Re[s] = 0.5
    imag_part = np.arange(imag_start, imag_end, imag_gap) * 1j
    s = imag_part + 0.5
    ax = np.imag(s)  # this is what will be used for x-axis of plot

    # Make figure now so user doesn't have to stare at blank screen too much longer
    plt.figure()
    # Reduce margins
    plt.tight_layout()
    # Make even smaller margins
    plt.subplots_adjust(left=0.05, right=0.99, top=0.95, bottom=0.05)

    # plt.subplot(2,1,1)
    plt.axhline(color='k')
    plt.title('Normalized to gamma magnitude')

    # When running in GUI, plots don't show up right away. Normally, you can call plt.show(),
    # but this pauses GUI, so we use plt.draw() followed by brief plt.pause() instead.
    def show():
        plt.draw()
        plt.pause(0.001)


    # This will make a blank window show up so user knows something is about to happen
    show()

    #
    # Calculate riemann values at line Re[s] = 0.5
    #
    # We want to calculate this:
    #     y = -np.real[rm.riemann_symmetric(s) / abs(gamma(s / 2))]
    #
    # Note that this is slightly wasteful because riemann_symmetric multiplies by gamma and s(s-1), and
    # we end up dividing that back out.
    #
    # Also, because gamma(s) rounds to 0 for very large s, we calculate logs instead. Note that
    # log(abs(gamma(s)) = Re[log(gamma(s))], as given by derivation below:
    #
    # For any complex c = a+bi
    #     log(c) = log(abs(c)) + i*Arg(c)
    # Hence, for any complex c,
    #     log(abs(c)) = Re(log(c))

    t0 = time.time()
    rm.make_powers(s)
    t1 = time.time()
    print("Calculated k^s along critical line Re[s]=0.5 in %1.2f seconds" % (t1 - t0))

    y = np.real(rm.riemann_real(s, is_vertical=True))

    t2 = time.time()
    print("Calculated values along critical line Re[s]=0.5 in %1.2f seconds" % (t2 - t1))

    plt.plot(ax, y, linewidth=1)
    show()

    # Input vector along line Re[s] = 5
    if show_off_critical:
        s2 = imag_part + real_offset

        # Calculate at line Re[s] = offset
        y2 = np.real((np.pi ** (real_offset / 2)) * rm.riemann_real(s2, is_vertical=True))
        plt.plot(ax, -y2, linewidth=1)
        t3 = time.time()
        print("Plotted values along Re[s] = %1.1f in %1.2f seconds" % (real_offset, t3 - t2))
        show()

    #
    #  Now generate histograms
    #
    list_minima = ((y[1:-1] < y[:-2]) & (y[1:-1] <= y[2:])).nonzero()[0]
    list_maxima = ((y[1:-1] > y[:-2]) & (y[1:-1] >= y[2:])).nonzero()[0]

    def plot_hist_pair(mins, maxes, numbins=200):
        values_min = y[mins[1:] + 1]  # Exclude first minimum. (Not necessary)
        values_max = y[maxes[1:] + 1]  # Exclude first maximum, as it is the only one below zero

        hist_min = np.histogram(values_min, bins=numbins, density=True)
        hist_max = np.histogram(values_max, bins=numbins, density=True)

        plt.plot(hist_min[1][1:], hist_min[0])
        plt.plot(hist_max[1][1:], hist_max[0])

    plt.figure()

    print("Found " + str(len(list_minima)) + " minima and " + str(len(list_maxima)) + " maxima")
    # Plot histogram for first 500 minima/maxima
    plot_hist_pair(list_minima[0:500], list_maxima[0:500], numbins=40)
    # Plot histogram for first 10k minima/maxima
    plot_hist_pair(list_minima[0:10000], list_maxima[0:10000], numbins=120)
    plot_hist_pair(list_minima, list_maxima)
    plt.title("Riemann minima and maxima " + str(imag_end))

    show()

if __name__ == "__main__":
    print("Done importing")

    show_critial_line()

    # This keeps main window open until dismissed. Don't call this if we came here from elsewhere.
    print("Close window to exit program")
    plt.show()
