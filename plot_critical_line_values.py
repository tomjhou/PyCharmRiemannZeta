import numpy as np
from matplotlib import use
import matplotlib.pyplot as plt
import time
import pickle
from tkinter import filedialog
import riemann_math as rm

use("TkAgg")


# When running in GUI, plots don't show up right away. Normally, you can call plt.show(),
# but this pauses GUI, so we use plt.draw() followed by brief plt.pause() instead.
def show():
    plt.draw()
    plt.pause(0.001)


def show_critical_line(show_graph1=True, show_off_critical=False, save_full=True):
    imag_start = 1
    imag_end = 120000
    imag_points = (imag_end - imag_start) * 100 + 1

    real_offset = 5

    iter_limit = imag_end + 2000

    if rm.RIEMANN_ITER_LIMIT < iter_limit:
        rm.RIEMANN_ITER_LIMIT = iter_limit  # int(height * 3)

    print(f"Precomputing coefficients, iteration limit = {rm.RIEMANN_ITER_LIMIT}. Please wait...")
    rm.precompute_coeffs()

    # Documentation says to use linspace for floats, arange for integers. Ok, whatever.
    # Notably, linspace includes endpoint, arange does not.
    imag_part = np.linspace(imag_start, imag_end, imag_points)

    # Create input vector along critical line Re[s] = 0.5
    s = imag_part * 1j + 0.5
    ax = np.imag(s)  # this is what will be used for x-axis of plot

    # Make figure now so user doesn't have to stare at blank screen too much longer
    if show_graph1:
        fig1 = plt.figure()
        # Reduce margins
        plt.tight_layout()
        # Make even smaller margins
        plt.subplots_adjust(left=0.05, right=0.99, top=0.95, bottom=0.05)

        # Default window is small and short. Make it wide and short
        plotsize = plt.gcf().get_size_inches()
        fig1.set_size_inches(plotsize[0] * 4, plotsize[1])

        fig1.canvas.manager.window.wm_geometry("+25+50")

        # plt.subplot(2,1,1)
        plt.axhline(color='k')
        plt.title('Normalized to gamma magnitude')

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

    print("Calculating Riemann values along critical line. Please wait ...")
    y: np.ndarray = np.real(rm.riemann_real(s, is_vertical=True))

    t2 = time.time()
    print("Calculated values along critical line Re[s]=0.5 in %1.2f seconds" % (t2 - t1))

    if show_graph1:
        plt.plot(ax, y, linewidth=1)
        show()

        # Input vector along line Re[s] = 5
        if show_off_critical:
            s2 = imag_part * 1j + real_offset

            # Calculate at line Re[s] = offset
            y2 = np.real((np.pi ** (real_offset / 2)) * rm.riemann_real(s2, is_vertical=True))
            plt.plot(ax, -y2, linewidth=1)
            t3 = time.time()
            print("Plotted values along Re[s] = %1.1f in %1.2f seconds" % (real_offset, t3 - t2))
            show()

    #
    #  nonzero() returns indices of local min and max.
    #
    list_minima: np.ndarray = ((y[1:-1] < y[:-2]) & (y[1:-1] <= y[2:]))
    list_minima = list_minima.nonzero()[0]
    list_maxima: np.ndarray = ((y[1:-1] > y[:-2]) & (y[1:-1] >= y[2:]))
    list_maxima = list_maxima.nonzero()[0]

    print("Found " + str(len(list_minima)) + " minima and " + str(len(list_maxima)) + " maxima")

    if save_full:
        #  Save results including raw values y. These files can be very large
        filename = 'riemann_extrema_' + str(imag_end) + '_full.pkl'
        with open(filename, 'wb') as f:  # Python 3: open(..., 'wb')
            pickle.dump([imag_start, imag_end, y, list_minima, list_maxima], f)

    #  Save max/min values only, excluding raw values y (much smaller file)
    filename = 'riemann_extrema_' + str(imag_end) + '.pkl'
    values_min = y[list_minima + 1]
    values_max = y[list_maxima + 1]

    with open(filename, 'wb') as f:  # Python 3: open(..., 'wb')
        pickle.dump([imag_start, imag_end, values_min, values_max], f)

    show_histograms(imag_start, imag_end, values_min, values_max)


def show_histograms(imag_start=1, imag_end=np.nan, values_min=None, values_max=None, load_full=False):

    if values_min is None:
        values_min = []
    if values_max is None:
        values_max = []
    if len(values_min) == 0:
        #            file_path = 'riemann_extrema_120000_full.pkl'
        file_path = filedialog.askopenfilename()
        with open(file_path, 'rb') as f:
            if load_full:
                # Load raw values into vector y. These files might be very large
                imag_start, imag_end, y, list_minima, list_maxima = pickle.load(f)
                values_min = y[list_minima + 1]
                values_max = y[list_maxima + 1]
            else:
                # Don't load raw values. These files are much smaller
                imag_start, imag_end, values_min, values_max = pickle.load(f)

    def plot_hist_pair(mins, maxes, numbins=360, log=False):
        v_min = mins[1:]  # Exclude first minimum. (Not necessary)
        v_max = maxes[1:]  # Exclude first maximum, as it is the only one below zero

        hist_min = np.histogram(v_min, bins=numbins, density=True)
        hist_max = np.histogram(v_max, bins=numbins, density=True)

        if log:
            plt.plot(hist_min[1][1:], np.log10(hist_min[0]))
            plt.plot(hist_max[1][1:], np.log10(hist_max[0]))
        else:
            plt.plot(hist_min[1][1:], hist_min[0])
            plt.plot(hist_max[1][1:], hist_max[0])

    fig2 = plt.figure()
    fig2.canvas.manager.window.wm_geometry("+25+625")

    print("Plotting " + str(len(values_min)) + " minima and " + str(len(values_max)) + " maxima")

    # Plot histogram for first 10k minima/maxima
    plot_hist_pair(values_min[0:10000], values_max[0:10000], numbins=120)
    # Plot all min/max
    plot_hist_pair(values_min, values_max, numbins=240)
    plt.title(f"Distribution of Riemann local min/max for Imag(s) = {imag_start}-{imag_end}")

    # Plot logs
    fig3 = plt.figure()
    fig3.canvas.manager.window.wm_geometry("+625+625")
    plot_hist_pair(values_min[0:10000], values_max[0:10000], numbins=120, log=True)
    plot_hist_pair(values_min, values_max, numbins=240, log=True)

    plt.title(f"Log distribution of Riemann local min/max for Imag(s) = {imag_start}-{imag_end}")
    show()


if __name__ == "__main__":
    print("Done importing")

    show_critical_line()

    # This keeps main window open until dismissed. Don't call this if we came here from elsewhere.
    print("Close window to exit program")
    plt.show()
