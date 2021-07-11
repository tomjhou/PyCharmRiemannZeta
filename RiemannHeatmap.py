#
# Code is copied from here: https://nbviewer.jupyter.org/github/empet/Math/blob/master/DomainColoring.ipynb
#
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
    S = s * np.ones(H.shape)  # Saturation = 1.0 always

    mag = np.absolute(w)  #
    elas = 0.2 # Elasticity ... higher numbers give faster transition from dark to light
    range = 0.97  # Diff between max and min intensity. Should be slightly less than 1.0, otherwise low values will go completely black
    V = 1 - range / (1 + mag ** elas)

    # Generate white "spokes" to accentuate pinwheel effect
    numSpokes = 6
    spokes = np.mod(H*numSpokes,1)
    spokes = np.where(spokes > 0.5, spokes - 1, spokes)
    spokes = np.abs(spokes)

    # Create "spokes" variable that ranges from 0 to 1, with 1 being max spoke effect
    spokeWidth = 0.20
    spokes = 1 - spokes/spokeWidth
    # Remove negative values
    spokes = np.where(spokes < 0, 0, spokes)

    # Squared value causes spoke effect to turn on more gradually at edges
    spokes = spokes * spokes

    # Saturation becomes zero at spokes
    S = 1 - spokes
    # Intensity becomes 1 at spokes
    V = np.where(spokes < 1, V + (1 - V) * spokes, V)

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
        sum = sum + np.power(x + 1, -s)
    return sum


def make_plot(selection, fig2):
    print('Please wait while computing heatmap... (this may take a minute or two)')

    xCenter = 0
    yCenter = 0

    if selection <= 6:
        rm.precompute_coeffs()

    if selection == 1:
        # Standard version
        rsize = 30;
        plot_domain2(lambda z: rm.Riemann(z),
                     re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize * 5, yCenter + rsize * 5],
                     title='Riemann($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT),
                     N=120/rsize)
    elif selection == 2:
        # Standard version, square
        rsize = 40;
        plot_domain2(lambda z: rm.Riemann(z),
                     re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize, yCenter + rsize],
                     title='Riemann($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT),
                     N=240/rsize)
    elif selection == 3:
        # Standard, zoomed into 0.5 + 50j
        rsize = 2;
        xCenter = 0.5;
        yCenter = 50;
        plot_domain2(lambda z: rm.Riemann(z),
                     re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize * 5, yCenter + rsize * 5],
                     title='Riemann($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT),
                     N=50/rsize)
    elif selection == 4:
        # Symmetric version, critical strip
        rsize = 30;
        plot_domain2(lambda z: rm.RiemannSymmetric(z), re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize * 5, yCenter + rsize * 5],
                     title='RiemannSymmetric($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT),
                     N=200/rsize)
    elif selection == 5:
        # Symmetric version, square
        rsize = 20;
        plot_domain2(lambda z: rm.RiemannSymmetric(z), re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize, yCenter + rsize],
                     title='RiemannSymmetric($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT),
                     N=300/rsize)
    elif selection == 6:
        # Symmetric version, zoomed into 0.5 + 50j
        rsize = 2;
        xCenter = 0.5;
        yCenter = 50;
        plot_domain2(lambda z: rm.RiemannSymmetric(z),
                     re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize * 5, yCenter + rsize * 5],
                     title='RiemannSymmetric($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT),
                     N=50/rsize)
    elif selection == 7:
        # Gamma(s/2) function
        rsize = 30;
        plot_domain2(lambda z: gamma(z / 2),
                     re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize * 5, yCenter + rsize * 5],
                     title='gamma($z/2$)',
                     N=100/rsize)
    elif selection == 8:
        # Gamma(s/2) function, square
        rsize = 40;
        plot_domain2(lambda z: gamma(z / 2),
                     re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize, yCenter + rsize],
                     title='gamma($z/2$)',
                     N=400/rsize)
    elif selection == 9:
        # Gamma(s/2) function, zoomed into 0.5 + 50j
        rsize = 2;
        xCenter = 0.5;
        yCenter = 50;
        plot_domain2(lambda z: gamma(z / 2),
                     re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize * 5, yCenter + rsize * 5],
                     title='gamma($z/2$)',
                     N=50/rsize)
    elif selection == 10:
        # sine function, square
        rsize = 10;
        plot_domain2(lambda z: np.sin(z / 2),
                     re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize, yCenter + rsize],
                     title='sin($z$)',
                     N=200/rsize)
    elif selection == 11:
        # cosine function, square
        rsize = 20;
        plot_domain2(lambda z: np.cos(z / 2),
                     re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize, yCenter + rsize],
                     title='cos($z$)',
                     N=200/rsize)
    elif selection == 12:
        # complex exponential function, square
        rsize = 40;
        plot_domain2(lambda z: np.power(2, z),
                     re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize, yCenter + rsize],
                     title='2^$z$',
                     N=200/rsize)
    elif selection == 13:
        # complex power function, square
        rsize = 40;
        plot_domain2(lambda z: np.power(z, 3-1j),
                     re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize, yCenter + rsize],
                     title='$z$^(3-1j)',
                     N=200 / rsize)
    elif selection == 14:
        # pi ^ s/2
        rsize = 40;
        plot_domain2(lambda z: np.power(np.pi,z/2),
                     re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize, yCenter + rsize],
                     title='$pi$^(z/2)',
                     N=200/rsize)
    elif selection == 15:
        # pi ^ s/2
        rsize = 60;
        global nextFlag, qFlag

        for partialSum in range(20, 20, 5):
            plot_domain2(lambda z: RiemannPartial(z, partialSum),
                     re=[xCenter - rsize, xCenter + rsize],
                     im=[yCenter - rsize, yCenter + rsize],
                     title='Riemann partial sum ' + str(partialSum),
                     N=160/rsize)

            plt.pause(0.5)

            continue

            nextFlag = False
            while True:
                if qFlag:
                    break
                if nextFlag:
                    break
                plt.pause(0.05)

plt.rcParams['figure.figsize'] = 6, 6
plt.tight_layout()

rsize = 10
plot_domain2(lambda z: z, re=[-rsize, rsize], im=[-rsize, rsize], title='$z$', N=160/rsize)
plt.gca().figure.canvas.mpl_connect('key_press_event', on_keypress)

plt.pause(0.001)

while True:
    print('Select plot type:\n'
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
          '13. Complex power\n'
          '14. Exponential: pi^(z/2)\n'
          '15. Riemann partial sum')
    selection = input("Select type: ")
    selection = int(selection)

    t1 = time.time()

    # Create blank figure
    plt.rcParams['figure.figsize'] = 12, 12
    fig2 = plt.figure()
    fig2.canvas.mpl_connect('key_press_event', on_keypress)

    # Make smaller margins
    plt.subplots_adjust(left=0.04, right=0.99, top=0.99, bottom=0.05)

    # Draw plot
    make_plot(selection,fig2)

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
