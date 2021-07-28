#
# Plots complex functions as heatmaps
#
# Domain coloring code is copied from here:
# https://nbviewer.jupyter.org/github/empet/Math/blob/master/DomainColoring.ipynb
#
import tkinter
import typing

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import hsv_to_rgb
from scipy.special import gamma
import time

import riemann_math as rm


print('Done importing libraries')

mesh_points = 0  # Number of mesh points. Used to convert progress calls to percent
qFlag = False    # Becomes true when user clicks "quit" button. Causes computations to wind down gracefully and stop
update_progress_callback: typing.Callable[[float], None]  # = None

# Counts number of points processed during long computations. Used to update progress bar
points_processed = 0


class Settings:
    pass


# Static variables (this class is never instantiated)
Settings.SCALE = .98       # Plot size as a function of screen height
Settings.MESH_DENSITY = 0.35  # Set # of mesh points as a function of "standard"
Settings.oversample = False
Settings.REUSE_FIGURE = False
Settings.phase_only = False
Settings.top_only = True
Settings.last_selection = -1  # Save last graph selection so we can recalculate
Settings.auto_recalculate = False
Settings.oversample = False
Settings.parameterA = 1.0        # Parameter for formula
Settings.parameterB = 0.0


def make_fig_plot(_event, _id):
    if _id < 0:
        # Negative ID is the recalculate button
        _id = Settings.last_selection
        Settings.REUSE_FIGURE = True
    else:
        Settings.last_selection = _id
    make_plot(_id)
    plt.pause(.001)


def do_quit(_event):
    global qFlag

    rm.quit_computation_flag = True
    qFlag = True


def update_computation_status(increment=1):
    global points_processed

    points_processed += increment

    if update_progress_callback is not None:
        update_progress_callback(points_processed * 100 / mesh_points)
    else:
        if np.mod(points_processed, 20000) == 0:
            print("   \r" + str(int(points_processed / 1000)) + "k ", end='')  # Print status every 10k samples


mpl.use('TkAgg')  # Generally I prefer this. It is faster, and seems more polisehd than Qt5Agg
# mpl.use('Qt5Agg')  # Need install PyQt5. Temporary window won't close in MacOS. Windows resize when moved.

import mpl_figure_manager as mfm
fig_mgr: mfm.MplFigureManager   # = None

# New button-based GUI selection method
plot_list = {
    0: "Pinwheel",
    1: "Riemann, strip",
    2: "Riemann, square20",
    4: "Symmetric Riemann, strip",
    5: "Symmetric Riemann, square2",
    19: "Symmetric Riemann, square20",
    7: "Gamma, strip",
    8: "Gamma, square4",
    18: "Gamma, square20",
    10: "sin(az) + b",
    11: "cos(az) + b",
    12: "Exponential a^z + b",
    13: "Power z^(a+bi)",
    14: "Exp: pi^(-z/2)",
    15: "Riemann partial sum",
    16: "1 - 2 ^ (1-z)",
    17: "Dirichlet eta",
    20: "Riemann * Gamma"
}

checkbox_list = ["Auto recalculate",
                 "Im>0 only",
                 "Phase only"]

if __name__ == "__main__":
    print("\nPlease run main.py instead of this file. Thank you!!")


if False:

    # Old text-based selection method
    print('Select plot type:\n'
          '3. Standard Riemann, zoomed critical strip\n'
          '6. Symmetric Riemann, zoomed critical strip\n'
          '9. Gamma, zoomed critical strip\n'
          '14. Exponential: pi^(z/2)\n')


# Computes hue corresponding to complex number z. Return value is in range 0 - 1
def Hcomplex(z):

    h = np.angle(z) / (2 * np.pi) + 1  # Add 1 to avoid negative values.
    # mod(h,1) returns remainder, i.e. portion after decimal point
    return np.mod(h, 1)


def g(x):
    return (1 - 1 / (1 + x ** 2)) ** 0.2


def fmt(n: float):
    if n >= 500000:
        return '{:1.2f}'.format(n / 1000000) + 'M'

    if n >= 1000:
        return '{:1.2f}'.format(n / 1000) + 'k'

    return str(n)


# Evaluates complex function f at nodes of grid determined by re, im, and N
#   re=(a, b) and im=(c, d), are pairs of real numbers that define limits of rectangular region.
#   N is a real number specifying grid points per unit interval.
def eval_grid(f, re, im, n):
    global mesh_points

    w = re[1] - re[0]
    h = im[1] - im[0]
    res_w = n * w  # horizontal resolution
    res_h = n * h  # vertical resolution
    x = np.linspace(re[0], re[1], int(res_w))
    y = np.linspace(im[0], im[1], int(res_h))
    x, y = np.meshgrid(x, y)
    z = x + 1j * y

    print('Mesh has ' + str(int(res_w)) + " x " + str(int(res_h)) + " = " + fmt(np.size(z)) + ' points')

    mesh_points = np.size(z)

    return f(z), z


# This is the color-generating function that is passed to plot_domain().
# Maps complex number to HSV.
# Hue is based on phase angle,
# Saturation is constant,
# Intensity is based on magnitude, 0 = black, 1 = 0.5, infinity is white.
def color_to_HSV(w, max_saturation):  # Classical domain coloring
    # w is the  array of values f(z)
    # s is the constant saturation

    global Settings

    H = Hcomplex(w)  # Determine hue
    # S = s * np.ones(H.shape)  # Saturation = 1.0 always

    mag = np.absolute(w)  #
    mag = np.where(mag < 1e-323, 1e-323, mag)  # To avoid divide by zero errors, impose a min magnitude

    if Settings.phase_only:
        v = 1
    else:
        elas = 0.1  # Elasticity ... higher numbers give faster transition from dark to light
        intensity_range = 0.97  # Diff between max, min intensity. Should be <1.0, otherwise low values go black

        # Create steps in V, so we will get magnitude contours representing factors of 10 change in magnitude
        step_size = 10
        mag = np.log(mag) / np.log(step_size)
        mag = np.floor(mag)
        mag = np.power(step_size, mag)

        # V = mag ^ 0.2 / (1 + mag ^ 0.2). This has the property that V(1/m) = 1 - V(m).
        # Hence, V[f^-1(s)] = 1 - V[f(s)]
        v = 1 - intensity_range / (1 + mag ** elas)

    # Generate white "spokes" to accentuate pinwheel effect
    num_spokes = 6
    spokes = np.mod(H * num_spokes, 1)
    spokes = np.where(spokes > 0.5, spokes - 1, spokes)
    spokes = np.abs(spokes)

    # Create "spokes" variable that ranges from 0 to 1, with 1 being max spoke effect
    spoke_width = 0.067
    spokes = 1 - spokes / spoke_width
    # Remove negative values
    spokes = np.where(spokes < 0, 0, spokes)

    # Squared value causes spoke effect to turn on more gradually at edges
    spokes = spokes * spokes

    # Saturation becomes zero when spokes == 1
    sat = (1 - spokes) * max_saturation
    # Intensity becomes 1 when spokes == 1
    v = np.where(spokes <= 1, v + (1 - v) * spokes, v)

    # V = np.ones(H.shape)
    # the points mapped to infinity are colored with white; hsv_to_rgb(0, 0, 1)=(1, 1, 1)=white

    HSV = np.dstack((H, sat, v))
    RGB = hsv_to_rgb(HSV)
    return RGB


# Creates figure, then calls plot_domain with standard options. When done, calculates
# points per second calculation speed
def plot_domain2(f, re=(-1, 1), im=(-1, 1), title=''):  # Number of points per unit interval)
    global Settings, mesh_points, points_processed

    if not hasattr(Settings, 'last_figure'):
        # Last figure never created. May have just launched program.
        Settings.REUSE_FIGURE = False
    elif not plt.fignum_exists(Settings.last_figure.number):
        # Last figure no longer exists. May have been closed by user.
        Settings.REUSE_FIGURE = False

    if Settings.REUSE_FIGURE:
        plt.figure(Settings.last_figure.number)
        # Clear flag so that we only reuse once. If we want to
        # reuse again, caller must set flag again
        Settings.REUSE_FIGURE = False
    else:
        # Create new figure, and bind to keypress handler
        aspect = abs(re[1] - re[0]) / (im[1] - im[0])
        fig = fig_mgr.make_plot_fig(Settings.SCALE * aspect, Settings.SCALE)
        fig.canvas.mpl_connect('key_press_event', on_keypress)
        Settings.last_figure = fig
        fig.canvas.draw()

    points_processed = 0
    if update_progress_callback is not None:
        update_progress_callback(0)

    t1 = time.time()

    density = fig_mgr.screen_y_pixels / abs(im[1] - im[0]) * Settings.MESH_DENSITY
    if Settings.oversample:
        density = density * 4

    # Mesh points is calculated here. Prior to calling this, it will be zero.
    mesh_points = plot_domain(color_to_HSV, f, re, im, title, 1, density, True)

    # Report time delay
    if mesh_points > 0:
        delay = time.time() - t1
        print('Completed domain coloring plot in ' + str(delay) + ' seconds')
        rate = mesh_points / delay
        print('Rate: ' + fmt(rate) + ' points per second')

        try:
            fig_mgr.canvas_plot.draw()
#            fig_mgr.canvas_plot.flush_events()   # Doesn't seem to be needed
        except tkinter.TclError:
            print("Previous canvas (fig_mgr.canvas_plot) invalid. Did you close previous window?")
            # If we recalculated but closed the original plot, then there will be a new figure generated
            # automatically, but fig_mgr.canvas_plot will not be valid, and we need plt.pause() to show figure
#            plt.pause(0.001)

            # Get current figure and save it, so that next recalculate will work normally.
            fig_mgr.fig_plot = plt.gcf()
            fig_mgr.canvas_plot = fig_mgr.fig_plot.canvas

            # Show figure
            fig_mgr.fig_plot.show()

        if update_progress_callback is not None:
            update_progress_callback(100)


def plot_domain(color_func, f, re=(-1, 1), im=(-1, 1), title='',
                saturation=0.9,  # Saturation
                density=200,  # Number of points per unit interval
                show_axis=None):  # Whether to show axis

    # Evaluate function f on a grid of points
    w, mesh = eval_grid(f, re, im, density)

    if rm.quit_computation_flag:
        print("Calculation cancelled by user after " + str(points_processed) + " of " +
              str(mesh_points) + " points = " + '{:1.2f}'.format(100 * points_processed / mesh_points) + " %")
        return 0

    # Convert results to color
    domc = color_func(w, saturation)
    plt.xlabel("$\Re(z)$")
    plt.ylabel("$\Im(z)$")
    plt.title(title)
    if show_axis:
        plt.imshow(domc, origin="lower",
                   extent=[re[0], re[1], im[0], im[1]])

    else:
        plt.imshow(domc, origin="lower")
        plt.axis('off')

    return np.size(mesh)


def on_keypress(event):
    global qFlag
    if event.key == "x":
        qFlag = True
        print("User pressed x on keyboard")


# This works whether s is single number or vector
def riemann_partial_sum(s, partial_sum_limit):

    r_sum = 1

    if partial_sum_limit <= 1:
        return 1

    for x in range(1, partial_sum_limit):

        if np.mod(x, 2) == 0:
            # Note that even though index is EVEN, the denominator is ODD,
            # since it is 1/((x+1)^s)
            r_sum += np.power(x + 1, -s)
        else:
            # ODD index for EVEN denominators
            r_sum -= np.power(x + 1, -s)

        if update_progress_callback is not None:
            update_progress_callback(100*(x+1)/partial_sum_limit)

        if rm.quit_computation_flag:
            print("Partial sum cancelled by user after " + str(x) + " of " + str(partial_sum_limit) + " iterations")
            break

    return r_sum


def make_plot(_selection):

    x_center = 0
    y_center = 0
    rm.quit_computation_flag = False

    if (0 < _selection <= 6) or (_selection == 17) or (_selection == 19):
        rm.precompute_coeffs()

    print('Please wait while computing heatmap... (this may take a minute or two)')

    y_max = 30
    if Settings.top_only:
        y_min = 0
    else:
        y_min = -30

    a = Settings.parameterA
    b = Settings.parameterB

    if _selection == 0:
        mesh_size = 10
        plot_domain2(lambda z: z, re=[-mesh_size, mesh_size], im=[-mesh_size, mesh_size], title='$z$')
    elif _selection == 1:
        # Standard version, critical strip centered at (0,0)
        plot_domain2(lambda z: rm.riemann(z),
                     re=[x_center - 1, x_center + 2],
                     im=[y_min, y_max],
                     title='Riemann($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT))
    elif _selection == 2:
        # Standard version, square
        mesh_size = 20
        plot_domain2(lambda z: rm.riemann(z),
                     re=[x_center - mesh_size, x_center + mesh_size],
                     im=[y_center - mesh_size, y_center + mesh_size],
                     title='Riemann($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT))
    elif _selection == 3:
        # Standard, zoomed into 0.5 + 50j
        mesh_size = 10
        x_center = 0.5
        y_center = 50
        plot_domain2(lambda z: rm.riemann(z),
                     re=[x_center - mesh_size/5, x_center + mesh_size/5],
                     im=[y_center - mesh_size, y_center + mesh_size],
                     title='Riemann($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT))
    elif _selection == 4:
        # Symmetric version, critical strip
        plot_domain2(lambda z: rm.RiemannSymmetric(z),
                     re=[x_center - 1, x_center + 2],
                     im=[y_min, y_max],
                     title='RiemannSymmetric($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT))
    elif _selection == 5:
        # Symmetric version, square
        mesh_size = 2
        plot_domain2(lambda z: rm.RiemannSymmetric(z),
                     re=[x_center - mesh_size, x_center + mesh_size],
                     im=[y_center - mesh_size, y_center + mesh_size],
                     title='RiemannSymmetric($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT))
    elif _selection == 19:
        # Symmetric version, square
        mesh_size = 20
        plot_domain2(lambda z: rm.RiemannSymmetric(z),
                     re=[x_center - mesh_size, x_center + mesh_size],
                     im=[y_center - mesh_size, y_center + mesh_size],
                     title='RiemannSymmetric($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT))
    elif _selection == 6:
        # Symmetric version, zoomed into 0.5 + 50j
        mesh_size = 10
        x_center = 0.5
        y_center = 50
        plot_domain2(lambda z: rm.RiemannSymmetric(z),
                     re=[x_center - mesh_size/5, x_center + mesh_size/5],
                     im=[y_center - mesh_size, y_center + mesh_size],
                     title='RiemannSymmetric($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT))
    elif _selection == 7:
        # Gamma(s/2) function, vertical strip centered at (0,0)
        plot_domain2(lambda z: gamma(z / 2),
                     re=[x_center - 1, x_center + 2],
                     im=[y_min, y_max],
                     title='gamma($z/2$)')
    elif _selection == 8:
        # Gamma(s/2) function, square ±4, ±4
        mesh_size = 4
        plot_domain2(lambda z: gamma(z),
                     re=[x_center - mesh_size, x_center + mesh_size],
                     im=[y_center - mesh_size, y_center + mesh_size],
                     title='gamma($z$)')
    elif _selection == 18:
        # Gamma(s/2) function, square ±20, ±20
        mesh_size = 20
        plot_domain2(lambda z: rm.gamma_with_progress(z),
                     re=[x_center - mesh_size, x_center + mesh_size],
                     im=[y_center - mesh_size, y_center + mesh_size],
                     title='gamma($z$)')
    elif _selection == 9:
        # Gamma(s/2) function, zoomed in, centered at 0.5 + 50j
        mesh_size = 2
        x_center = 0.5
        y_center = 50
        plot_domain2(lambda z: gamma(z / 2),
                     re=[x_center - mesh_size, x_center + mesh_size],
                     im=[y_center - mesh_size * 5, y_center + mesh_size * 5],
                     title='gamma($z/2$)')
    elif _selection == 10:
        # sine function, square
        mesh_size = 10
        plot_domain2(lambda z: np.sin(a * z) + b,
                     re=[x_center - mesh_size, x_center + mesh_size],
                     im=[y_center - mesh_size, y_center + mesh_size],
                     title='sin(a*$z$)+b')
    elif _selection == 11:
        # cosine function, square
        mesh_size = 20
        plot_domain2(lambda z: np.cos(a * z) + b,
                     re=[x_center - mesh_size, x_center + mesh_size],
                     im=[y_center - mesh_size, y_center + mesh_size],
                     title='cos(a*$z$)+b')
    elif _selection == 12:
        # complex exponential function, square
        mesh_size = 40
        plot_domain2(lambda z: np.power(a, z) + b,
                     re=[x_center - mesh_size, x_center + mesh_size],
                     im=[y_center - mesh_size, y_center + mesh_size],
                     title='{:1.2f}'.format(a) + '^$z$' + ' + {:1.2f}'.format(b))
    elif _selection == 13:
        # complex power function, square
        mesh_size = 40
        plot_domain2(lambda z: np.power(z, a - 1j*b),
                     re=[x_center - mesh_size, x_center + mesh_size],
                     im=[y_center - mesh_size, y_center + mesh_size],
                     title='$z$^(' + '{:1.2f}'.format(a) + ' + ' + '{:1.2f}'.format(b) + 'j)')
    elif _selection == 14:
        #  pi ^ (z/2), one of the two functions multiplied by Riemann to get symmetric Riemann
        mesh_size = 40
        plot_domain2(lambda z: np.power(np.pi, -z/2),
                     re=[x_center - mesh_size, x_center + mesh_size],
                     im=[y_center - mesh_size, y_center + mesh_size],
                     title='3^$z$')
    elif _selection == 15:
        # Partial summation of Dirichlet eta function (alternating Riemann)
        # This sum converges for Re(s) > 0

        partial_sum_limit = 2
        for x in range(0, 20):

            plot_domain2(lambda z: riemann_partial_sum(z, partial_sum_limit),
                         re=[x_center - 1, x_center + 3],
                         im=[y_min, y_max],
                         title='Riemann partial sum ' + str(partial_sum_limit))

            Settings.REUSE_FIGURE = True

            partial_sum_limit = partial_sum_limit * 2

            # plt.pause(0.2)

            if rm.quit_computation_flag:
                break

        Settings.REUSE_FIGURE = False

    elif _selection == 16:
        # This is the function that converts Riemann zeta to Dirichlet eta function
        # 1 - 2 ^ (1-s)
        mesh_size = 30
        plot_domain2(lambda z: 1 - np.power(2, 1 - z),
                     re=[x_center - mesh_size, x_center + mesh_size],
                     im=[y_center - mesh_size, y_center + mesh_size],
                     title='1-2^(1-$z$)')
    elif _selection == 17:
        #  Dirichlet Eta instead of Riemann zeta
        mesh_size = 20
        plot_domain2(lambda z: rm.riemann(z, do_eta=True),
                     re=[x_center - mesh_size, x_center + mesh_size],
                     im=[y_center - mesh_size, y_center + mesh_size],
#                     re=[x_center - 1, x_center + 3],
#                     im=[y_min, y_max],
                     title='Eta($z$)')
    elif _selection == 20:
        mesh_size = 60
        plot_domain2(lambda z: rm.RiemannGamma(z),
                     re=[x_center - mesh_size, x_center + mesh_size],
                     im=[y_center - mesh_size, y_center + mesh_size],
                     title='Riemann(z) * Gamma($z$)')

    rm.quit_computation_flag = False

    print()  # Riemann prints updates periodically - add newline in case Riemann did not
