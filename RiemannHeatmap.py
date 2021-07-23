#
# Code is copied from here: https://nbviewer.jupyter.org/github/empet/Math/blob/master/DomainColoring.ipynb
#
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import hsv_to_rgb
import RiemannMath as rm
from scipy.special import gamma
import time

rm.RIEMANN_ITER_LIMIT = 80

print('Done importing libraries')

# Computes hue corresponding to complex number z. Return value is in range 0 - 1
def Hcomplex(z):
    H = np.angle(z) / (2 * np.pi) + 1   # Add 1 to avoid negative values.
    return np.mod(H, 1)  # mod returns remainder. This is equivalent to taking decimal portion


def g(x):
    return (1 - 1 / (1 + x ** 2)) ** 0.2

# Evaluates complex function f at nodes of grid determined by re, im, and N
#   re=(a, b) and im=(c, d), are pairs of real numbers that define limits of rectangular region.
#   N is a real number specifying grid points per unit interval.
def eval_grid(f, re, im, N):

    l = re[1] - re[0]
    h = im[1] - im[0]
    resL = N * l  # horizontal resolution
    resH = N * h  # vertical resolution
    x = np.linspace(re[0], re[1], int(resL))
    y = np.linspace(im[0], im[1], int(resH))
    x, y = np.meshgrid(x, y)
    z = x + 1j * y

    rm.callCount = 0
    print('Mesh has ' + str(np.size(z)/1000) + 'k points')
    return f(z)


# This is the color-generating function that is passed to plot_domain().
# Maps complex number to HSV.
# Hue is based on phase angle,
# Saturation is constant,
# Intensity is based on magnitude, 0 = black, 1 = 0.5, infinity is white.
def color_to_HSV(w, s):  # Classical domain coloring
    # w is the  array of values f(z)
    # s is the constant saturation

    H = Hcomplex(w)  # Determine hue
    # S = s * np.ones(H.shape)  # Saturation = 1.0 always

    mag = np.absolute(w)  #
    mag = np.where(mag < 1e-323, 1e-323, mag)  # To avoid divide by zero errors, impose a min magnitude
    elas = 0.1 # Elasticity ... higher numbers give faster transition from dark to light
    range = 0.97  # Diff between max and min intensity. Should be slightly less than 1.0, otherwise low values will go completely black

    # Create steps in V, so we will get magnitude contours representing factors of 10 change in magnitude
    stepSize = 10
    mag = np.log(mag) / np.log(stepSize)
    mag = np.floor(mag)
    mag = np.power(stepSize, mag)

    # V = mag ^ 0.2 / (1 + mag ^ 0.2). This has the property that V(1/m) = 1 - V(m).
    # Hence, V[f^-1(s)] = 1 - V[f(s)]
    V = 1 - range / (1 + mag ** elas)

    # Generate white "spokes" to accentuate pinwheel effect
    numSpokes = 6
    spokes = np.mod(H*numSpokes,1)
    spokes = np.where(spokes > 0.5, spokes - 1, spokes)
    spokes = np.abs(spokes)

    # Create "spokes" variable that ranges from 0 to 1, with 1 being max spoke effect
    spokeWidth = 0.067
    spokes = 1 - spokes/spokeWidth
    # Remove negative values
    spokes = np.where(spokes < 0, 0, spokes)

    # Squared value causes spoke effect to turn on more gradually at edges
    spokes = spokes * spokes

    # Saturation becomes zero when spokes == 1
    S = 1 - spokes
    # Intensity becomes 1 when spokes == 1
    V = np.where(spokes <= 1, V + (1 - V) * spokes, V)

    # V = np.ones(H.shape)
    # the points mapped to infinity are colored with white; hsv_to_rgb(0, 0, 1)=(1, 1, 1)=white

    HSV = np.dstack((H, S, V))
    RGB = hsv_to_rgb(HSV)
    return RGB

def plot_domain2(f, re=[-1, 1], im=[-1, 1], title='', N=200):  # Number of points per unit interval)
    plot_domain(color_to_HSV, f, re, im, title, 1, N, True)

def plot_domain(color_func, f, re=[-1, 1], im=[-1, 1], title='',
                s=0.9,  # Saturation
                N=200,  # Number of points per unit interval
                daxis=None): # Whether to show axis

    # Evaluate function f on a grid of points
    w = eval_grid(f, re, im, N)
    domc = color_func(w, s)
    plt.xlabel("$\Re(z)$")
    plt.ylabel("$\Im(z)$")
    plt.title(title)
    if (daxis):
        plt.imshow(domc, origin="lower",
                   extent=[re[0], re[1], im[0], im[1]])

    else:
        plt.imshow(domc, origin="lower")
        plt.axis('off')


qFlag = False
nextFlag = False

def on_keypress(event):
    global qFlag, nextFlag
    if event.key == "x":
        qFlag = True
    if event.key == "n":
        nextFlag = True

def RiemannPartial(s, partialSum):

    sum = 1

    if partialSum <= 1:
        return 1

    for x in range(1,partialSum):
        if np.mod(x,2) == 0:
            # Note that even though index is EVEN, the denominator is ODD,
            # since it is 1/((x+1)^s)
            sum = sum + np.power(x + 1, -s)
        else:
            # ODD index for EVEN denominators
            sum = sum - np.power(x + 1, -s)
    return sum


def make_plot(selection, screen_y):
    xCenter = 0
    yCenter = 0

    if selection > 0 and selection <= 6:
        rm.precompute_coeffs()

    print('Please wait while computing heatmap... (this may take a minute or two)')

    if selection == 0:
        rsize = 10
        plot_domain2(lambda z: z, re=[-rsize, rsize], im=[-rsize, rsize], title='$z$', N=screen_y / 2 / rsize)
    elif selection == 1:
        # Standard version, critical strip centered at (0,0)
        rsize = 10;
        plot_domain2(lambda z: rm.Riemann(z),
                     re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize * 5, yCenter + rsize * 5],
                     title='Riemann($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT),
                     N=screen_y/rsize/5)
    elif selection == 2:
        # Standard version, square
        rsize = 40;
        plot_domain2(lambda z: rm.Riemann(z),
                     re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize, yCenter + rsize],
                     title='Riemann($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT),
                     N=screen_y/rsize)
    elif selection == 3:
        # Standard, zoomed into 0.5 + 50j
        rsize = 2;
        xCenter = 0.5;
        yCenter = 50;
        plot_domain2(lambda z: rm.Riemann(z),
                     re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize * 5, yCenter + rsize * 5],
                     title='Riemann($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT),
                     N=screen_y/rsize/5)
    elif selection == 4:
        # Symmetric version, critical strip
        rsize = 30;
        plot_domain2(lambda z: rm.RiemannSymmetric(z), re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize * 5, yCenter + rsize * 5],
                     title='RiemannSymmetric($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT),
                     N=screen_y/rsize/5)
    elif selection == 5:
        # Symmetric version, square
        rsize = 2;
        plot_domain2(lambda z: rm.RiemannSymmetric(z), re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize, yCenter + rsize],
                     title='RiemannSymmetric($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT),
                     N=screen_y/rsize)
    elif selection == 6:
        # Symmetric version, zoomed into 0.5 + 50j
        rsize = 2;
        xCenter = 0.5;
        yCenter = 50;
        plot_domain2(lambda z: rm.RiemannSymmetric(z),
                     re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize * 5, yCenter + rsize * 5],
                     title='RiemannSymmetric($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT),
                     N=screen_y/rsize/5)
    elif selection == 7:
        # Gamma(s/2) function, vertical strip centered at (0,0)
        rsize = 30;
        plot_domain2(lambda z: gamma(z / 2),
                     re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize * 5, yCenter + rsize * 5],
                     title='gamma($z/2$)',
                     N=screen_y/rsize/5)
    elif selection == 8:
        # Gamma(s/2) function, square
        rsize = 4;
        plot_domain2(lambda z: gamma(z),
                     re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize, yCenter + rsize],
                     title='gamma($z$)',
                     N=screen_y * 1.5/rsize)
    elif selection == 9:
        # Gamma(s/2) function, zoomed in, centered at 0.5 + 50j
        rsize = 2;
        xCenter = 0.5;
        yCenter = 50;
        plot_domain2(lambda z: gamma(z / 2),
                     re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize * 5, yCenter + rsize * 5],
                     title='gamma($z/2$)',
                     N=screen_y/rsize/5)
    elif selection == 10:
        # sine function, square
        rsize = 10;
        plot_domain2(lambda z: np.sin(z / 2),
                     re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize, yCenter + rsize],
                     title='sin($z$)',
                     N=screen_y/rsize)
    elif selection == 11:
        # cosine function, square
        rsize = 20;
        plot_domain2(lambda z: np.cos(z / 2),
                     re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize, yCenter + rsize],
                     title='cos($z$)',
                     N=screen_y/rsize)
    elif selection == 12:
        # complex exponential function, square
        rsize = 40;
        plot_domain2(lambda z: np.power(3, z),
                     re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize, yCenter + rsize],
                     title='3^$z$',
                     N=screen_y/rsize)
    elif selection == 13:
        # complex power function, square
        rsize = 40;
        plot_domain2(lambda z: np.power(z, 3-8j),
                     re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize, yCenter + rsize],
                     title='$z$^(3-8j)',
                     N=screen_y / rsize)
    elif selection == 14:
        # pi ^ s/2
        rsize = 40;
        plot_domain2(lambda z: np.power(np.pi,z/2),
                     re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize, yCenter + rsize],
                     title='$pi$^(z/2)',
                     N=screen_y/rsize)
    elif selection == 15:
        # Partial summation of Dirichlet eta function (alternating Riemann)
        rsize = 20;
        global nextFlag, qFlag

        partialSum = 2
        for x in range(1, 100):
            plot_domain2(lambda z: RiemannPartial(z, partialSum),
                     re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize, yCenter + rsize],
                     title='Riemann partial sum ' + str(partialSum),
                     N=screen_y/4/rsize)

            partialSum = partialSum * 2
            plt.pause(0.2)

            continue

    elif selection == 16:
        # This is the function that converts Riemann zeta to Dirichlet eta function
        # 1 - 2 ^ (1-s)
        rsize = 20;
        plot_domain2(lambda z: 1 - np.power(2, 1 - z),
                     re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize, yCenter + rsize],
                     title='1-2^(1-$z$)',
                     N=screen_y * 1.5 / rsize)

    print()   # Riemann prints updates periodically - add newline in case Riemann did not

print("Available backends:")
print(mpl.rcsetup.all_backends)
mpl.use('TkAgg')    # Generally I prefer this. It is faster, and temporary window seems to close more reliably (on Win and MacOS)
# mpl.use('Qt5Agg')  # Sometimes this is slower. Need to install PyQt5 for this to work. Temporary window won't close under MacOS
backend = mpl.get_backend()
print("Matplotlib backend is: " + backend) # Returns Qt5Agg after installing PyQt5 ... if you don't have Qt5, I think it returns TkAgg something

# Get screen height in pixels. Somehow, this always requires creating a temporary window. Is there something more elegant?
mgr = plt.get_current_fig_manager()
if backend == "Qt5Agg":
    mgr.full_screen_toggle()  # This creates a temporary full-screen window. Works for Qt5Agg but not TkAgg
    screen_y = mgr.canvas.height()
elif backend == "TkAgg":
    screen_x, screen_y = mgr.window.wm_maxsize()  # Creates temporary full-screen window. Works for TkAgg but not Qt5Agg
else:
    # Make default assumption
    print('Unable to determine screen resolution. Will default to 1024')
    screen_y = 1024

plt.close()  # Close temporary window. On MacOS and Qt5, this doesn't seem to work until next window is created. On Windows or TkAgg, it is fine.
plt.pause(0.01)  # This is not necessary.

print('Screen y resolution is: ' + str(screen_y))
screen_y = screen_y - 50  # Subtract a small amount or else the toolbar at bottom will mess things up.

while True:
    print('Select plot type:\n'
          '0. Pinwheel\n'
          '1. Standard Riemann, critical strip (30 sec)\n'
          '2. Standard Riemann, large square (31 sec)\n'
          '3. Standard Riemann, zoomed critical strip\n'
          '4. Symmetric Riemann, critical strip (47 sec)\n'
          '5. Symmetric Riemann, large square (31 sec)\n'
          '6. Symmetric Riemann, zoomed critical strip\n'
          '7. Gamma, critical strip\n'
          '8. Gamma, large square\n'
          '9. Gamma, zoomed critical strip\n'
          '10. sin function, large square\n'
          '11. cos function, large square\n'
          '12. Complex exponential\n'
          '13. Complex power spirals\n'
          '14. Exponential: pi^(z/2)\n'
          '15. Riemann partial sum\n'
          '16. 1 - 2 ^ (1-s)')
    selection = input("Select type: ")
    selection = int(selection)

    t1 = time.time()

    # Create blank figure
    plt.rcParams['figure.figsize'] = 12, 12
    fig2 = plt.figure()
    fig2.canvas.mpl_connect('key_press_event', on_keypress)
#    plt.gca().figure.canvas.mpl_connect('key_press_event', on_keypress)

    # Make smaller margins
    plt.tight_layout()
    # Make even smaller margins
    plt.subplots_adjust(left=0.05, right=0.99, top=0.95, bottom=0.05)

    # Draw plot
    make_plot(selection, screen_y)

    # Report time delay
    print('Completed domain coloring plot in ' + str(time.time() - t1) + ' seconds')

    print("Press x to exit, n for next plot (focus must be on Riemann window)")
    nextFlag = False
    while True:
      if qFlag:
          break
      if nextFlag:
          break
      plt.pause(0.05)

    if qFlag:
        break
