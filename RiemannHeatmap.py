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
import ButtonManager as bm

rm.RIEMANN_ITER_LIMIT = 80
USE_BUTTON_PANEL = True

class settings:
    pass

settings.SCALE = .98       # Plot size as a function of screen height
settings.MESH_DENSITY = 0.25  # Set # of mesh points as a function of "standard"
settings.oversample = False
settings.REUSE_FIGURE = False
settings.phase_only = False
settings.top_only = False
settings.last_selection = 0  # Save last graph selection so we can recalculate


print('Done importing libraries')


# Computes hue corresponding to complex number z. Return value is in range 0 - 1
def Hcomplex(z):
    h = np.angle(z) / (2 * np.pi) + 1  # Add 1 to avoid negative values.
    return np.mod(h, 1)  # mod returns remainder. This is equivalent to taking decimal portion


def g(x):
    return (1 - 1 / (1 + x ** 2)) ** 0.2


def fmt(n: np.float64):
    if n >= 500000:
        return '{:1.2f}'.format(n / 1000000) + 'M'

    if n >= 1000:
        return '{:1.2f}'.format(n / 1000) + 'k'

    return str(n)


# Evaluates complex function f at nodes of grid determined by re, im, and N
#   re=(a, b) and im=(c, d), are pairs of real numbers that define limits of rectangular region.
#   N is a real number specifying grid points per unit interval.
def eval_grid(f, re, im, n):
    l = re[1] - re[0]
    h = im[1] - im[0]
    res_l = n * l  # horizontal resolution
    res_h = n * h  # vertical resolution
    x = np.linspace(re[0], re[1], int(res_l))
    y = np.linspace(im[0], im[1], int(res_h))
    x, y = np.meshgrid(x, y)
    z = x + 1j * y

    rm.callCount = 0
    print('Mesh has ' + str(int(res_l)) + " x " + str(int(res_h)) + " = " + fmt(np.size(z)) + ' points')
    return f(z), z


# This is the color-generating function that is passed to plot_domain().
# Maps complex number to HSV.
# Hue is based on phase angle,
# Saturation is constant,
# Intensity is based on magnitude, 0 = black, 1 = 0.5, infinity is white.
def color_to_HSV(w, s):  # Classical domain coloring
    # w is the  array of values f(z)
    # s is the constant saturation

    global settings

    H = Hcomplex(w)  # Determine hue
    # S = s * np.ones(H.shape)  # Saturation = 1.0 always

    mag = np.absolute(w)  #
    mag = np.where(mag < 1e-323, 1e-323, mag)  # To avoid divide by zero errors, impose a min magnitude

    if settings.phase_only:
        v = 1
    else:
        elas = 0.1  # Elasticity ... higher numbers give faster transition from dark to light
        intensity_range = 0.97  # Diff between max, min intensity. Should be slightly <1.0, otherwise low values go black

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
    S = 1 - spokes
    # Intensity becomes 1 when spokes == 1
    v = np.where(spokes <= 1, v + (1 - v) * spokes, v)

    # V = np.ones(H.shape)
    # the points mapped to infinity are colored with white; hsv_to_rgb(0, 0, 1)=(1, 1, 1)=white

    HSV = np.dstack((H, S, v))
    RGB = hsv_to_rgb(HSV)
    return RGB



# Creates figure, then calls plot_domain with standard options
def plot_domain2(f, re=(-1, 1), im=(-1, 1), title=''):  # Number of points per unit interval)
    global settings

    if settings.REUSE_FIGURE:
        plt.figure(settings.last_figure.number)
        # Clear flag so that we only reuse once. If we want to
        # reuse again, caller must set flag again
        settings.REUSE_FIGURE = False
    else:
        # Create new figure, and bind to keypress handler
        aspect = abs(re[1] - re[0]) / (im[1] - im[0])
        fig = bmgr.make_plot_fig(settings.SCALE * aspect, settings.SCALE)
        fig.canvas.mpl_connect('key_press_event', on_keypress)
        settings.last_figure = fig
        plt.pause(.001)


    t1 = time.time()

    density = bmgr.screen_y_pixels / abs(im[1] - im[0]) * settings.MESH_DENSITY
    mesh_size = plot_domain(color_to_HSV, f, re, im, title, 1, density, True)

    # Report time delay
    delay = time.time() - t1
    print('Completed domain coloring plot in ' + str(delay) + ' seconds')
    rate = mesh_size / delay
    print('Rate: ' + fmt(rate) + ' points per second')


def plot_domain(color_func, f, re=(-1, 1), im=(-1, 1), title='',
                s=0.9,  # Saturation
                density=200,  # Number of points per unit interval
                show_axis=None):  # Whether to show axis

    # Evaluate function f on a grid of points
    w, mesh = eval_grid(f, re, im, density)
    domc = color_func(w, s)
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


qFlag = False
next_flag = False


def on_keypress(event):
    global qFlag, next_flag
    if event.key == "x":
        qFlag = True
        print("User pressed x on keyboard")
    if event.key == "n":
        next_flag = True


def RiemannPartial(s, partial_sum):
    r_sum = 1

    if partial_sum <= 1:
        return 1

    for x in range(1, partial_sum):
        if np.mod(x, 2) == 0:
            # Note that even though index is EVEN, the denominator is ODD,
            # since it is 1/((x+1)^s)
            r_sum = r_sum + np.power(x + 1, -s)
        else:
            # ODD index for EVEN denominators
            r_sum = r_sum - np.power(x + 1, -s)
    return r_sum


def make_plot(_selection, screen_y):

    x_center = 0
    y_center = 0

    if (0 < _selection <= 6) or (_selection == 17):
        rm.precompute_coeffs()

    print('Please wait while computing heatmap... (this may take a minute or two)')

    y_max = 30
    if settings.top_only:
        y_min = 0
    else:
        y_min = -30

    # Number of mesh dots per unit interval ... roughly one per pixel
    density = screen_y / (y_max-y_min)

    if _selection == 0:
        mesh_size = 10
        plot_domain2(lambda z: z, re=[-mesh_size, mesh_size], im=[-mesh_size, mesh_size], title='$z$')
    elif _selection == 1:
        # Standard version, critical strip centered at (0,0)
        mesh_size = 30
        plot_domain2(lambda z: rm.Riemann(z),
                     re=[x_center, x_center + 2],
                     im=[y_min, y_max],
                     title='Riemann($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT))
    elif _selection == 2:
        # Standard version, square
        mesh_size = 40
        plot_domain2(lambda z: rm.Riemann(z),
                     re=[x_center - mesh_size, x_center + mesh_size],
                     im=[y_min, y_max],
                     title='Riemann($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT))
    elif _selection == 3:
        # Standard, zoomed into 0.5 + 50j
        mesh_size = 10
        x_center = 0.5
        y_center = 50
        plot_domain2(lambda z: rm.Riemann(z),
                     re=[x_center - mesh_size/5, x_center + mesh_size/5],
                     im=[y_center - mesh_size, y_center + mesh_size],
                     title='Riemann($z$), iter = ' + str(rm.RIEMANN_ITER_LIMIT))
    elif _selection == 4:
        # Symmetric version, critical strip
        mesh_size = 30
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
        mesh_size = 30
        plot_domain2(lambda z: gamma(z / 2),
                     re=[x_center - mesh_size, x_center + mesh_size],
                     im=[y_min, y_max],
                     title='gamma($z/2$)')
    elif _selection == 8:
        # Gamma(s/2) function, square
        mesh_size = 4
        plot_domain2(lambda z: gamma(z),
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
        plot_domain2(lambda z: np.sin(z / 2),
                     re=[x_center - mesh_size, x_center + mesh_size],
                     im=[y_center - mesh_size, y_center + mesh_size],
                     title='sin($z$)')
    elif _selection == 11:
        # cosine function, square
        mesh_size = 20
        plot_domain2(lambda z: np.cos(z / 2),
                     re=[x_center - mesh_size, x_center + mesh_size],
                     im=[y_center - mesh_size, y_center + mesh_size],
                     title='cos($z$)')
    elif _selection == 12:
        # complex exponential function, square
        mesh_size = 40
        plot_domain2(lambda z: np.power(3, z),
                     re=[x_center - mesh_size, x_center + mesh_size],
                     im=[y_center - mesh_size, y_center + mesh_size],
                     title='3^$z$')
    elif _selection == 13:
        # complex power function, square
        mesh_size = 40
        plot_domain2(lambda z: np.power(z, 3 - 8j),
                     re=[x_center - mesh_size, x_center + mesh_size],
                     im=[y_center - mesh_size, y_center + mesh_size],
                     title='$z$^(3-8j)')
    elif _selection == 14:
        # pi ^ s/2
        mesh_size = 40
        plot_domain2(lambda z: np.power(np.pi, z / 2),
                     re=[x_center - mesh_size, x_center + mesh_size],
                     im=[y_center - mesh_size, y_center + mesh_size],
                     title='$pi$^(z/2)')
    elif _selection == 15:
        # Partial summation of Dirichlet eta function (alternating Riemann)
        # This sum converges for Re(s) > 0
        mesh_size = 30
        global next_flag, qFlag, REUSE_FIGURE

        partial_sum = 2
        for x in range(1, 100):
            plot_domain2(lambda z: RiemannPartial(z, partial_sum),
                         re=[x_center, x_center + 2],
                         im=[y_min, y_max],
                         title='Riemann partial sum ' + str(partial_sum))

            settings.REUSE_FIGURE = True

            partial_sum = partial_sum * 2
            plt.pause(0.2)

            if next_flag:
                break

        settings.REUSE_FIGURE = False

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
        plot_domain2(lambda z: rm.Riemann(z, do_eta=True),
                     re=[x_center - 1, x_center + 2],
                     im=[y_min, y_max],
                     title='Eta($s$)')

    print()  # Riemann prints updates periodically - add newline in case Riemann did not


def make_fig_plot(event, id):
    if id < 0:
        id = settings.last_selection
        settings.REUSE_FIGURE = True
    else:
        settings.last_selection = id
    make_plot(id, bmgr.screen_y_pixels * settings.SCALE * settings.MESH_DENSITY)
    plt.pause(.001)


def do_quit(event):
    global qFlag
    qFlag = True


def do_submit(text, _id):
    print("Object id: " + str(_id) + ", entered: " + text)


def slider_update(val):
    #    print("slider: " + str(val))
    settings.MESH_DENSITY = val


def do_checkbox(label):
    checkbox_list

    index = checkbox_list.index(label)
    if index == 0:
        settings.top_only = not settings.top_only
    elif index == 1:
        settings.phase_only = not settings.phase_only


def do_slider(val, id):
    slider_val[id] = val


button_list = [None] * 20
slider_list = [None] * 20
text_box_list = [None] * 20
slider_val = [0] * 20


def AddIdButton(text, _id):
    button_list[_id] = bmgr.add_id_button(text, _id)
    button_list[_id].on_clicked(lambda x: make_fig_plot(x, button_list[_id].id))


def AddIdSlider2(text, _id):
    slider_list[_id] = bmgr.add_id_slider(text, _id)
    slider_list[_id].on_changed(lambda x: do_slider(x, slider_list[_id].id))


def AddIdTextBox(text, _id):
    text_box_list[_id] = bmgr.add_textbox(text, _id)
    text_box_list[_id].on_submit(lambda x: do_submit(x, text_box_list[_id].id))


mpl.use('TkAgg')  # Generally I prefer this. It is faster, and seems more polisehd than Qt5Agg
#mpl.use('Qt5Agg')  # Need install PyQt5. Temporary window won't close in MacOS. Windows resize when moved.


bmgr = None

# New button-based GUI selection method
plot_list = {
    0: "Pinwheel",
    1: "Riemann, strip",
    4: "Symmetric Riemann, strip",
    8: "Gamma",
    10: "Sine",
    11: "Cosine",
    12: "Exponential",
    13: "Power (spiral)",
    15: "Riemann partial sum",
    16: "1 - 2 ^ (1-s)",
    17: "Dirichlet eta"
}

checkbox_list = ["Im>0 only", "Phase only"]

bmgr = bm.ButtonManager(16,2)

#bmgr.increment_row()
#bmgr.increment_row()
cb1 = bmgr.add_checkbox(checkbox_list)
cb1.on_clicked(do_checkbox)

for k in plot_list:
    AddIdButton(plot_list[k], k)

b2 = bmgr.add_standard_button("Quit")
b2.on_clicked(do_quit)

#
# Start second column of widgets
#
bmgr.reset_button_coord(1)

b3 = bmgr.add_id_slider("mesh density",
                        _id=0, valmin=0, valmax=1,
                        valinit = settings.MESH_DENSITY)
b3.on_changed(slider_update)

b4 = AddIdButton("Recalculate", -1)
bmgr.increment_row()

for k in plot_list:
    AddIdSlider2("", k)

# Button figure also responds to key press - "x" to exit
bmgr.canvas2.mpl_connect('key_press_event', on_keypress)

plt.pause(0.01)

while not qFlag:
    bmgr.canvas2.flush_events()
    #        plt.pause(0.001) # This causes new plots to permanently overlay older ones
    time.sleep(0.025)  # This allows plots to be moved around more freely


if False:

    # Old text-based selection method
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
        next_flag = False
        while True:
            if qFlag:
                break
            if next_flag:
                break
            plt.pause(0.05)

        if qFlag:
            break
