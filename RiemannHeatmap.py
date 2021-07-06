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
rm.precompute_coeffs()


# Computes hue corresponding to complex number z. Return value is in range 0 - 1
def Hcomplex(z):
    H = np.angle(z) / (2 * np.pi) + 1   # Add 1 to avoid negative values.
    return np.mod(H, 1)  # mod returns remainder. This is equivalent to taking decimal portion


def g(x):
    return (1 - 1 / (1 + x ** 2)) ** 0.2

# Evaluates complex function f at nodes of grid determined by re, im, and N
#   re=(a, b) and im=(c, d), are pairs of real numbers that define limits of rectangular region.
#   N is a real number specifying grid points per unit interval.
def func_vals(f, re, im, N):

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
# Darkness value based on magnitude, 0 = black, 1 = 0.5, infinity is white.
def color_to_HSV(w, s):  # Classical domain coloring
    # w is the  array of values f(z)
    # s is the constant saturation

    H = Hcomplex(w)  # Determine hue
    S = s * np.ones(H.shape)  # Saturation = 1.0 always
    mag = np.absolute(w)  #
    elas = 0.2 # Elasticity ... higher numbers give faster transition from dark to light
    V = (1.0 - 1.0 / (1 + mag ** elas))
    # V = np.ones(H.shape)
    # the points mapped to infinity are colored with white; hsv_to_rgb(0, 0, 1)=(1, 1, 1)=white

    HSV = np.dstack((H, S, V))
    RGB = hsv_to_rgb(HSV)
    return RGB

def plot_domain2(f, re=[-1, 1], im=[-1, 1], title='', N=200):  # Number of points per unit interval)
    plot_domain(color_to_HSV, f, re, im, title, 0.5, N, True)

def plot_domain(color_func, f, re=[-1, 1], im=[-1, 1], title='',
                s=0.9,  # Saturation
                N=200,  # Number of points per unit interval
                daxis=None): # Whether to show axis
    w = func_vals(f, re, im, N)
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


def on_keypress(event):
    global qFlag
    if event.key == "x":
        qFlag = True

plt.rcParams['figure.figsize'] = 4, 4
plt.tight_layout()

rsize = 10
plot_domain2(lambda z: z, re=[-rsize, rsize], im=[-rsize, rsize], title='$z$', N=100)
plt.gca().figure.canvas.mpl_connect('key_press_event', on_keypress)

plt.pause(0.001)


print('Select plot type:')
print('1. Standard Riemann, critical strip (30 sec)')
print('2. Standard Riemann, large square (31 sec)')
print('3. Standard, zoomed critical strip')
print('4. Symmetric, critical strip (47 sec)')
print('5. Symmetric, large square (31 sec)')
print('6. Symmetric, zoomed critical strip')
print('7. Gamma, critical strip')
print('8. Gamma, large square')
print('9. Gamma, zoomed critical strip')
print('10. sin function, large square')
print('11. cos function, large square')
print('12. x^(3-1j)')
print('13. pi^(z/2)')
selection = input("Select type: ")
selection = int(selection)

t1 = time.time()
plt.rcParams['figure.figsize'] = 12, 12

print('Please wait while computing heatmap... (this may take a minute or two)')
fig2 = plt.figure()

xCenter = 0
yCenter = 0

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
    rsize = 60;
    plot_domain2(lambda z: rm.Riemann(z),
                 re=[xCenter - rsize, xCenter + rsize],
                 im=[yCenter - rsize, yCenter + rsize],
                 title='Riemann($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT),
                 N=160/rsize)
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
    rsize = 40;
    plot_domain2(lambda z: rm.RiemannSymmetric(z), re=[xCenter - rsize, xCenter + rsize],
                 im=[yCenter - rsize, yCenter + rsize],
                 title='RiemannSymmetric($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT),
                 N=120/rsize)
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
    rsize = 125;
    plot_domain2(lambda z: gamma(z / 2),
                 re=[xCenter - rsize, xCenter + rsize],
                 im=[yCenter - rsize, yCenter + rsize],
                 title='gamma($z/2$)',
                 N=200/rsize)
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
                 N=100/rsize)
elif selection == 11:
    # cosine function, square
    rsize = 20;
    plot_domain2(lambda z: np.cos(z / 2),
                 re=[xCenter - rsize, xCenter + rsize],
                 im=[yCenter - rsize, yCenter + rsize],
                 title='cos($z$)',
                 N=200/rsize)
elif selection == 12:
    # complex power function, square
    rsize = 40;
    plot_domain2(lambda z: np.power(z, 5j),
                 re=[xCenter - rsize, xCenter + rsize],
                 im=[yCenter - rsize, yCenter + rsize],
                 title='$z$^(5j)',
                 N=200/rsize)
elif selection == 13:
    # pi ^ s/2
    rsize = 40;
    plot_domain2(lambda z: np.power(np.pi,z/2),
                 re=[xCenter - rsize, xCenter + rsize],
                 im=[yCenter - rsize, yCenter + rsize],
                 title='$pi$^(z/2)',
                 N=200/rsize)
print('Completed domain coloring plot in ' + str(time.time() - t1) + ' seconds')

fig2.canvas.mpl_connect('key_press_event', on_keypress)

print("Press x to exit (focus must be on Riemann window)")
while not qFlag:
    plt.pause(0.05)
