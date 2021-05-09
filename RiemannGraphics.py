import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.widgets import Button
import numpy as np
from enum import Enum
from matplotlib.collections import LineCollection
import matplotlib.pylab as pl

import global_vars as gv
import RiemannMath as rm

# Button dimensions
BUTTON_WIDTH = 0.1
BUTTON_HEIGHT = 0.03     # Above about 0.1, buttons disappear - presumably they collide with graphs?
BUTTON_SPACING_X = BUTTON_WIDTH + 0.005
BUTTON_X = 0.05
BUTTON_Y = 0.95

# Text coordinates, horizontal
TEXT_INPUT_MATRIX_X = 0.001

# Text coordinates, vertical
TEXT_Y_COORD = 0.001  # Height of input matrix text above horizontal axis
TEXT_MATRIX_ROW_Y_SPACING = 0.08  # Vertical space between rows of input matrix
TEXT_Y_SPACING = 0.35  # Vertical space between input and output matrices

FONT_FAMILY = 'monospace'

INPUT_VECTOR_COLOR = 'r'
OUTPUT_VECTOR_COLOR = (0, 0, 0)  # Dark purple

AXIS_SCALES = [0.3, 1, 2.5]
flagAxisScaleIndex = 1      # Starting scale. Can be changed if user zooms in or out

INPUT_RANGES = [30, 100, 1000]
flagInputRangeIndex = 0

ANIMATE_SPEEDS = [0.0, 0.03, 0.1, 0.3, 1, 3, 10, 30]
flagAnimateSpeedIndex = 0
flagSavedAnimateSpeedIndex = 0  # When we pause, index is saved here so we can unpause and resume original speed

RESOLUTION_BOOST = 5    # Extra resolution when animation slows down
quitflag = False  # When true, program will exit
flagRecalc = False  # When true, will redraw existing graph items
flagRedrawAxes = False  # When true, will redraw existing graph axes
flagChangeMatrix = False  # When true, will update matrix and redraw
flagRestartAnimation = False
flagShowArrows = True
flagSideDots = False
flagCenterAtTarget = False

flagX = 0.5  # Indicates mouse x position
flagY = 0  # Indicates mouse y position
flagBidirectional = False
whichRowToAdjust = 0

localPlotSize = 0

flagMouseDown = False

arrowOutput = []

global bg1
global bg_blank

class PatchStyle(Enum):
    ARROW = 1
    CIRCLE = 2
    ARROW_CIRCLE = 3


def GetAxisScale():
    return AXIS_SCALES[flagAxisScaleIndex]


def get_animation_speed():
    return ANIMATE_SPEEDS[flagAnimateSpeedIndex]


# Set plot axis scaling to zoom in/out
def set_zoom():
    scaleTmp = GetAxisScale()
    pltHeight = localPlotSize * scaleTmp
    pltWidth = pltHeight * 1.5
    if scaleTmp > 1:
        plt.xlim([-pltWidth * .75, pltWidth * 1.25])
        plt.ylim([-pltHeight, pltHeight])
    elif scaleTmp <= 0.1 or flagCenterAtTarget:
        plt.xlim([-pltWidth, pltWidth])
        plt.ylim([-pltHeight, pltHeight])
    elif scaleTmp <= 0.3:
        plt.xlim([-pltWidth * 0.1, pltWidth * 1.9])
        plt.ylim([-pltHeight, pltHeight])
    else:
        plt.xlim([-pltWidth * .667, pltWidth * 1.333])
        plt.ylim([-pltHeight, pltHeight])



def Save_background(canvas):
    global bg1, flagAxisScaleIndex

    saveIndex = flagAxisScaleIndex
    bg1 = []
    # Save background bitmaps for each scale
    for flagAxisScaleIndex in range(0,len(AXIS_SCALES)):
        set_zoom()
        canvas.draw()
        bg1.append(canvas.copy_from_bbox(gv.ax1.bbox))

    # Restore original scale
    flagAxisScaleIndex = saveIndex
    set_zoom()
    canvas.draw()


def Restore_background(canvas):
    if flagCenterAtTarget:
        canvas.restore_region(bg_blank)
    else:
        canvas.restore_region(bg1[flagAxisScaleIndex])  # Restores static elements and erases background


def CheckArrows():
    if gv.localShowArrows != flagShowArrows:
        # Just switched on or off. Hide or show arrows.
        gv.localShowArrows = flagShowArrows
        setArrowVisible(gv.localShowArrows)


def MakeLine(PLOT_SIZE, INPUT_SCALE, RESOLUTION, LINE_X):
    V_LINE_HEIGHT = PLOT_SIZE * INPUT_SCALE
    V_LINE_SEGMENTS = int(RESOLUTION * V_LINE_HEIGHT)

    # Rainbow colors
    colors = pl.cm.jet(np.linspace(0, 1, V_LINE_SEGMENTS))

    max = LINE_X + 1j * V_LINE_HEIGHT

    if flagBidirectional:
        min = LINE_X - 1j * V_LINE_HEIGHT
    else:
        min = LINE_X + 0j

    lns = np.linspace(min, max, V_LINE_SEGMENTS)

    gv.line_min = min
    gv.line_max = max
    gv.line_length = np.abs(max - min)

    return lns, colors


def GetPoint(position):
    return gv.line_min + (gv.line_max - gv.line_min) * position / gv.line_length


# function to turn arrow visibility on and off
def setArrowVisible(b):
    for i in range(0, len(arrowOutput)):
        arrowOutput[i].set_visible(b)


def multiline(xs, ys, lc=None, color=None, **kwargs):
    """Plot lines with different colorings

    Parameters
    ----------
    xs : iterable container of x coordinates
    ys : iterable container of y coordinates
    color : iterable container of numbers mapped to colormap
    ax (optional): Axes to plot on.
    kwargs (optional): passed to LineCollection

    Notes:
        len(xs) == len(ys) == len(c) is the number of line segments
        len(xs[i]) == len(ys[i]) is the number of points for each line (indexed by i)

    Returns
    -------
    lc : LineCollection instance.
    """

    # create LineCollection
    sz = len(xs) - 1

    segments = np.zeros((sz, 2, 2))
    # Each line segment has start and end point. End of each segment is equivalent to start of next
    segments[:,0,0] = xs[0:sz]
    segments[:,0,1] = ys[0:sz]
    segments[:,1,0] = xs[1:sz+1]
    segments[:,1,1] = ys[1:sz+1]

    if lc is None:
        lc = LineCollection(segments, **kwargs)
        if color is None:
            # Initialize color scheme
            lc.set_array(np.arange(sz))
    else:
        # Update segments of existing line
        lc.set_segments(segments)
        if color is None:
            # Initialize color scheme
            lc.set_array(np.arange(sz))

    # set coloring of line segments
    #    Note: I get an error if I pass c as a list here... not sure why.
    if color is not None:
        color = color[0:sz, ]
        lc.set_color(np.asarray(color))

    return lc


# Mouse button press. Use this to start moving vector1 or vector2 in top-right plot
def on_mouse_press(event):
    if event.inaxes != gv.ax1:
        return
    global flagX, flagY, flagMouseDown, whichRowToAdjust, flagChangeMatrix
    if event.button == mpl.backend_bases.MouseButton.LEFT:
        flagMouseDown = True
        flagChangeMatrix = True
        whichRowToAdjust = 0
    if event.button == mpl.backend_bases.MouseButton.RIGHT:
        flagMouseDown = True
        flagChangeMatrix = True
        whichRowToAdjust = 1
    flagX = event.xdata
    flagY = event.ydata


def on_mouse_release(event):
    if event.inaxes != gv.ax1:
        return
    global flagMouseDown, flagChangeMatrix
    flagMouseDown = False
    flagChangeMatrix = False


def on_mouse_move(event):
    if event.inaxes != gv.ax1:
        return
    global flagMouseDown, flagRecalc, flagX, flagY, flagChangeMatrix
    if flagMouseDown:
        flagX = event.xdata
        flagY = event.ydata
        flagChangeMatrix = True


def do_1_vs_2(event=None):
    global flagBidirectional, flagRecalc
    flagBidirectional = not flagBidirectional  # Toggle between 1 and 2 rows
    flagRecalc = True


def do_quit(event=None):
    global quitflag
    quitflag = True


# Pause and unpause animation
def do_pause_unpause_animate(event=None):
    global flagAnimateSpeedIndex, flagSavedAnimateSpeedIndex

    if flagAnimateSpeedIndex == 0:
        flagAnimateSpeedIndex = flagSavedAnimateSpeedIndex if (flagSavedAnimateSpeedIndex > 0) else 1
    else:
        flagSavedAnimateSpeedIndex = flagAnimateSpeedIndex
        flagAnimateSpeedIndex = 0


def get_next_speed(s, inc = 1):
    s = s + inc
    if s >= len(ANIMATE_SPEEDS):
        s = len(ANIMATE_SPEEDS) - 1
    if s <= 0:
        s = 1
    return s


def do_speed_up(event=None):
    global flagAnimateSpeedIndex, flagSavedAnimateSpeedIndex

    if flagAnimateSpeedIndex == 0:
        flagSavedAnimateSpeedIndex = get_next_speed(flagSavedAnimateSpeedIndex)
    else:
        flagAnimateSpeedIndex = get_next_speed(flagAnimateSpeedIndex)


def do_slow_down(event=None):
    global flagAnimateSpeedIndex, flagSavedAnimateSpeedIndex

    if flagAnimateSpeedIndex == 0:
        flagSavedAnimateSpeedIndex = get_next_speed(flagSavedAnimateSpeedIndex, -1)
    else:
        flagAnimateSpeedIndex = get_next_speed(flagAnimateSpeedIndex, -1)


def do_reset_animation(event=None):
    global flagRestartAnimation, flagRecalc, flagChangeMatrix
    flagRestartAnimation = True
    flagChangeMatrix = True


def do_expand_input(event=None):
    global flagRecalc, flagInputRangeIndex, flagRedrawAxes
    flagRecalc = True
    flagInputRangeIndex = flagInputRangeIndex + 1
    if flagInputRangeIndex >= len(INPUT_RANGES):
        flagInputRangeIndex = 0
    flagRedrawAxes = True


def do_toggle_arrows(event=None):
    global flagRecalc, flagShowArrows
    flagShowArrows = not flagShowArrows


def do_zoom_out(event=None):
    global flagRedrawAxes, flagAxisScaleIndex, AXIS_SCALES
    flagAxisScaleIndex = flagAxisScaleIndex + 1
    if flagAxisScaleIndex >= len(AXIS_SCALES):
        flagAxisScaleIndex = 0
    flagRedrawAxes = True


def do_toggle_side_dots(event=None):
    global flagSideDots
    flagSideDots = not flagSideDots


def do_toggle_center(event=None):
    global flagCenterAtTarget, flagRedrawAxes
    flagCenterAtTarget = not flagCenterAtTarget
    flagRedrawAxes = True


# Keyboard press
def on_keypress(event):
    global flagChangeMatrix, whichRowToAdjust, flagX, flagChangeMatrix
    if event.key == "x":
        do_quit()
    elif event.key == " ":
        do_pause_unpause_animate()
    elif event.key == "r":
        do_reset_animation()
    elif event.key == "z":
        do_zoom_out()
    elif event.key == "v":
        do_speed_up()
    elif event.key == "up":
        do_speed_up()
    elif event.key == "down":
        do_slow_down()
    elif event.key == "left":
        flagX = flagX - 0.001
        flagChangeMatrix = True
        whichRowToAdjust = 0
    elif event.key == "right":
        flagX = flagX + 0.001
        flagChangeMatrix = True
        whichRowToAdjust = 0


# Format a single floating point number to have 3 digits after decimal
def fmt(n):
    return '{:+2.3f}'.format(n)


# Format a 2D vector as text inside brackets
def fmt_complex(r):
    return fmt(r.real) + ', ' + fmt(r.imag) + 'j'


class GraphicsObjects:
    def __init__(self, PLOT_SIZE, showbuttons = True):
        global bg_blank

        print("Available backends:")
        print(mpl.rcsetup.all_backends)
        mpl.use('TkAgg')   # Qt5Agg might also be available, but it is MUCH slower. Force TkAgg, which plots much faster
        self.backend = mpl.get_backend()
        print("Matplotlib backend is: " + self.backend) # Returns Qt5Agg after installing Qt5 ... if you don't have Qt5, I think it returns TkAgg something

        # Create figure
        fig = plt.figure() #figsize=(MAIN_WINDOW_SIZE, MAIN_WINDOW_SIZE))

        window = plt.get_current_fig_manager().window
        dpi = fig.dpi

        if self.backend == "Qt5Agg":
            # Need a hack to get screen size. Temporarily make a full-screen window, get its size, then later set "real" size
            window.showMaximized()  # Make window fullscreen
            plt.pause(.001)  # Draw items to screen so we can get size
            screen_x, screen_y = fig.get_size_inches() * fig.dpi  # size in pixels
        else:
            # window.state('zoomed')  # Make window fullscreen, for TkAgg
            screen_x, screen_y = window.wm_maxsize() # Get full scren monitor coordinates for TkAgg. Doesn't work under Qt5Agg

        screen_y = screen_y - 50  # Subtract a small amount or else the toolbar at bottom will mess things up.

        # Create window with aspect ratio 1.5:1.0
        fig.set_size_inches(screen_y * 1.5 / dpi, screen_y / dpi)  # Convert from pixels to inches by dividing by dpi

        if self.backend == "Qt5Agg":
            plt.pause(.001)  # Draw to screen so that new window size takes effect
#        fig.tight_layout() # This doesn't work, somehow

        canvas = fig.canvas
        canvas.mpl_connect('button_press_event', on_mouse_press)
        canvas.mpl_connect('button_release_event', on_mouse_release)
        canvas.mpl_connect('key_press_event', on_keypress)
        canvas.mpl_connect('motion_notify_event', on_mouse_move)
        self.fig = fig
        self.canvas = canvas
        self.ax = plt.gca()

        # Make smaller margins
        plt.subplots_adjust(left=0.04, right=0.99, top=0.99, bottom=0.05)

        canvas.draw()
        bg_blank = canvas.copy_from_bbox(self.ax.bbox)

        self.Add_axes(PLOT_SIZE)

        self.textObj = TextObjects(self.canvas, self.fig)

        if showbuttons:
            self.make_buttons()

        # Need this or else current axis will be the last button drawn, instead of main plot
        plt.sca(self.ax)


    # Create initial graphics, then save backgrounds
    def Add_axes(self, PLOT_SIZE):
        global bg1, localPlotSize

        plt.sca(self.ax)
        plt.show(block=False)

        localPlotSize = PLOT_SIZE
        set_zoom()  # Need to do this after sca, or else it won't take effect

        self.ax.axhline(y=0, color='gray')
        self.ax.axvline(x=0, color='gray')

        # Make smaller margins
        #    gv.ax1.margins(x=0.01) # This doesn't work. Use subplots_adjust() instead


    def make_buttons(self):

        # Create one new axis for each new button
        ax_animate = plt.axes([BUTTON_X, BUTTON_Y, BUTTON_WIDTH, BUTTON_HEIGHT])
        self.b_animate = Button(ax_animate, 'Toggle animate')
        self.b_animate.on_clicked(do_pause_unpause_animate)

        ax_1_vs_2 = plt.axes([BUTTON_X + BUTTON_SPACING_X, BUTTON_Y, BUTTON_WIDTH, BUTTON_HEIGHT])
        self.b_1_vs_2 = Button(ax_1_vs_2, 'Bidirectional')
        self.b_1_vs_2.on_clicked(do_1_vs_2)

        ax_output = plt.axes([BUTTON_X + BUTTON_SPACING_X * 2, BUTTON_Y, BUTTON_WIDTH, BUTTON_HEIGHT])
        self.b_output = Button(ax_output, 'Expand input')
        self.b_output.on_clicked(do_expand_input)

        ax_reset = plt.axes([BUTTON_X + BUTTON_SPACING_X * 3, BUTTON_Y, BUTTON_WIDTH, BUTTON_HEIGHT])
        self.b_reset = Button(ax_reset, 'Reset animation (r)')
        self.b_reset.on_clicked(do_reset_animation)

        ax_show_arrows = plt.axes([BUTTON_X + BUTTON_SPACING_X * 4, BUTTON_Y, BUTTON_WIDTH, BUTTON_HEIGHT])
        self.b_show_arrows = Button(ax_show_arrows, 'Toggle arrows')
        self.b_show_arrows.on_clicked(do_toggle_arrows)

        ax_button_zoom = plt.axes([BUTTON_X + BUTTON_SPACING_X * 5, BUTTON_Y, BUTTON_WIDTH, BUTTON_HEIGHT])
        self.b_zoom = Button(ax_button_zoom, 'Zoom out (z)')
        self.b_zoom.on_clicked(do_zoom_out)

        ax_button_dots = plt.axes([BUTTON_X + BUTTON_SPACING_X * 6, BUTTON_Y, BUTTON_WIDTH, BUTTON_HEIGHT])
        self.b_dots = Button(ax_button_dots, 'Toggle center')
        self.b_dots.on_clicked(do_toggle_center)

        ax_quit = plt.axes([BUTTON_X + BUTTON_SPACING_X * 7, BUTTON_Y, BUTTON_WIDTH, BUTTON_HEIGHT])
        self.b_quit = Button(ax_quit, 'Quit (x)')
        self.b_quit.on_clicked(do_quit)

    # Create input and output rainbow lines
    def InitAnimation(self):

        if gv.pltLineCollection is not None:
            gv.pltLineCollection.remove()
        # pre-generate entire line
        th2b = rm.Riemann(gv.th1)
        if flagAxisScaleIndex > 0:
            # Make lighter colors when scaling is largest
            colors2 = [1, 1, 1, 1] - ([1, 1, 1, 1] - gv.linecolors) / 4
            line2 = multiline(np.real(th2b), np.imag(th2b), color=colors2, linewidths=0.5)
        else:
            line2 = multiline(np.real(th2b), np.imag(th2b), color=gv.linecolors, linewidths=0.5)
        gv.pltLineCollection = gv.ax1.add_collection(line2)

        # Hide any items that are not static
        setArrowVisible(False)
        self.textObj.text_clear()

        # Save background
        Save_background(self.canvas)

        setArrowVisible(True)

        return gv.pltLineCollection, th2b, line2


class TextObjects:
    def __init__(self, _canvas, _figure):
        self.canvas = _canvas
        self.fig = _figure
        # Create text plot in top left
        self.ax_text = plt.gca()

        self.textObjArrayRow1 = self.ax_text.text(
            TEXT_INPUT_MATRIX_X,
            TEXT_Y_COORD,
            '',
            color='b',
            family=FONT_FAMILY)

        self.canvas.draw()  # Need this so that text will render to screen, before we capture background
        self.background = self.canvas.copy_from_bbox(self.ax_text.bbox)

    def update_input_complex(self, vector):
        self.textObjArrayRow1.set_text(fmt_complex(vector))

        # "draw_artist" draws new objects efficiently, without updating entire plot
        self.ax_text.draw_artist(self.textObjArrayRow1)

    def update_input_real(self, vector):
        self.textObjArrayRow1.set_text("x = " + fmt(vector))

        # "draw_artist" draws new objects efficiently, without updating entire plot
        self.ax_text.draw_artist(self.textObjArrayRow1)

    def text_clear(self):
        self.textObjArrayRow1.set_text("")