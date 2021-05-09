#
# Code is copied from here: https://nbviewer.jupyter.org/github/empet/Math/blob/master/DomainColoring.ipynb
#
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import hsv_to_rgb
from scipy.special import zeta

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


def domaincol_c(w, s):  # Classical domain coloring
    # w is the  array of values f(z)
    # s is the constant saturation

    H = Hcomplex(w)
    S = s * np.ones(H.shape)
    modul = np.absolute(w)
    V = (1.0 - 1.0 / (1 + modul ** 2)) ** 0.2
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


#plt.figure()
plt.rcParams['figure.figsize'] = 12, 10
ab = [-2, 2]
cd = [-2, 2]
plt.subplot(1, 2, 1)
f = lambda z: (z**3 - 1)/z
plot_domain(domaincol_c, f, re=ab, im= cd, title='$f(z)=(z^3-1)/z$', daxis=True)
plt.subplot(1,2,2)
plot_domain(domaincol_c, lambda z:z, re=[-7, 7], im=[-7, 7], title='$z$', daxis=True)
plt.tight_layout()

plt.figure()
plot_domain(domaincol_c, lambda z:1/(1 - 2**(1-z)), re=[-2, 2], im=[-2, 2], title='$z$', daxis=True)

while True:
    plt.pause(0.001)

