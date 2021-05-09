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


def Hcomplex(z):# computes the hue corresponding to the complex number z
    H = np.angle(z) / (2*np.pi) + 1
    return np.mod(H, 1)


def g(x):
    return (1- 1/(1+x**2))**0.2


def func_vals(f, re, im, N):  # evaluates the complex function at the nodes of the grid
    # re and im are  tuples, re=(a, b) and im=(c, d), defining the rectangular region
    # N is the number of discrete points per unit interval

    l = re[1] - re[0]
    h = im[1] - im[0]
    resL = N * l  # horizontal resolution
    resH = N * h  # vertical resolution
    x = np.linspace(re[0], re[1], int(resL))
    y = np.linspace(im[0], im[1], int(resH))
    x, y = np.meshgrid(x, y)
    z = x + 1j * y
    return f(z)

# Maps 0 to black, and infinity to white. Intermediate colors are
def domaincol_c(w, s):  # Classical domain coloring
    # w is the  array of values f(z)
    # s is the constant saturation

    H = Hcomplex(w) # Determine hue
    S = s * np.ones(H.shape) # Saturation = 1.0 always
    modul = np.absolute(w)   #
#    V = (1.0 - 1.0 / (1 + modul ** 2)) ** 0.2
    V = np.ones(H.shape)
    # the points mapped to infinity are colored with white; hsv_to_rgb(0, 0, 1)=(1, 1, 1)=white

    HSV = np.dstack((H, S, V))
    RGB = hsv_to_rgb(HSV)
    return RGB


def plot_domain(color_func, f, re=[-1, 1], im=[-1, 1], title='',
                s=0.9, N=200, daxis=None):
    w = func_vals(f, re, im, N)
    domc = color_func(w, s)
    plt.xlabel("$\Re(z)$")
    plt.ylabel("$\Im(z)$")
    plt.title(title)
    if (daxis):
        plt.imshow(domc, origin="lower", extent=[re[0], re[1], im[0], im[1]])

    else:
        plt.imshow(domc, origin="lower")
        plt.axis('off')

qFlag = False

def on_keypress(event):
    global qFlag
    if event.key == "x":
        qFlag = True


print('Select plot type:')
print('1. Standard Riemann, critical strip (30 sec)')
print('2. Standard Riemann, large square')
print('3. Symmetric, critical strip')
print('4. Symmetric, large square')
print('5. Standard, zoomed critical strip')
print('6. Symmetric, zoomed critical strip')
print('7. Gamma, critical strip')
print('8. Gamma, large square')
print('9. Gamma, zoomed critical strip')
selection = input("Select type: ")
selection = int(selection)

plt.rcParams['figure.figsize'] = 4, 4
plt.tight_layout()
rsize = 1
plot_domain(domaincol_c, lambda z:z, re=[-rsize, rsize], im=[-rsize, rsize], title='$z$', daxis=True, N = 100)
plt.gca().figure.canvas.mpl_connect('key_press_event', on_keypress)

plt.pause(0.001)

t1 = time.time()
plt.rcParams['figure.figsize'] = 12, 12


print('Please wait while computing heatmap... (this may take a minute or two)')
fig2 = plt.figure()

xCenter = 0
yCenter = 0

if selection == 1:
    # Standard version
    rsize = 30; plot_domain(domaincol_c, lambda z:rm.Riemann(z), re=[xCenter-rsize, xCenter+rsize], im=[yCenter - rsize*5, yCenter + rsize*5], title='Riemann($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT), daxis=True, N = 4)
elif selection == 2:
    # Standard version, square
    rsize = 125; plot_domain(domaincol_c, lambda z:rm.Riemann(z), re=[xCenter-rsize, xCenter+rsize], im=[yCenter - rsize, yCenter + rsize], title='Riemann($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT), daxis=True, N = 2)
elif selection == 3:
    # Symmetric version
    rsize = 30; plot_domain(domaincol_c, lambda z:rm.RiemannSymmetric(z), re=[xCenter-rsize, xCenter+rsize], im=[yCenter - rsize*5, yCenter + rsize*5], title='RiemannSymmetric($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT), daxis=True, N = 5)
elif selection == 4:
    # Symmetric version, square
    rsize = 125; plot_domain(domaincol_c, lambda z:rm.RiemannSymmetric(z), re=[xCenter-rsize, xCenter+rsize], im=[yCenter - rsize, yCenter + rsize], title='RiemannSymmetric($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT), daxis=True, N = 2)
elif selection == 5:
    # Standard, zoomed into 0.5 + 50j
    rsize = 2; xCenter = 0.5; yCenter = 50; plot_domain(domaincol_c, lambda z:rm.Riemann(z), re=[xCenter-rsize, xCenter+rsize], im=[yCenter - rsize*5, yCenter + rsize*5], title='Riemann($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT), daxis=True, N = 25)
elif selection == 6:
    # Symmetric version, zoomed into 0.5 + 50j
    rsize = 2; xCenter = 0.5; yCenter = 50; plot_domain(domaincol_c, lambda z:rm.RiemannSymmetric(z), re=[xCenter-rsize, xCenter+rsize], im=[yCenter - rsize*5, yCenter + rsize*5], title='RiemannSymmetric($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT), daxis=True, N = 25)
elif selection == 7:
    # Gamma(s/2) function
    rsize = 30; plot_domain(domaincol_c, lambda z:gamma(z/2), re=[xCenter-rsize, xCenter+rsize], im=[yCenter - rsize*5, yCenter + rsize*5], title='gamma($z/2$)', daxis=True, N = 3)
elif selection == 8:
    # Gamma(s/2) function, square
    rsize = 125; plot_domain(domaincol_c, lambda z:gamma(z/2), re=[xCenter-rsize, xCenter+rsize], im=[yCenter - rsize, yCenter + rsize], title='gamma($z/2$)', daxis=True, N = 2)
elif selection == 9:
    # Gamma(s/2) function, zoomed into 0.5 + 50j
    rsize = 2; xCenter = 0.5; yCenter = 50; plot_domain(domaincol_c, lambda z:gamma(z/2), re=[xCenter-rsize, xCenter+rsize], im=[yCenter - rsize*5, yCenter + rsize*5], title='gamma($z/2$)', daxis=True, N = 25)


print('Completed domain coloring plot in ' + str(time.time() - t1) + ' seconds')

fig2.canvas.mpl_connect('key_press_event', on_keypress)

print("Press x to exit (focus must be on Riemann window)")
while not qFlag:
    plt.pause(0.05)
