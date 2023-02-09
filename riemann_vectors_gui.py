import matplotlib.collections
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.widgets import Button
import numpy as np
from enum import Enum
from matplotlib.collections import LineCollection
import matplotlib.pylab as pl

import riemann_math as rm

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

AXIS_SCALES: list[float] = [0.3, 1, 2.5]

INPUT_RANGES: list[float] = [30, 100, 1000]

ANIMATE_SPEEDS: list[float] = [0.0, 0.03, 0.1, 0.3, 1, 3, 10, 30]
RESOLUTION_BOOST = 5    # Extra resolution when animation slows down


class PatchStyle(Enum):
    ARROW = 1
    CIRCLE = 2
    ARROW_CIRCLE = 3


def multiline(xs, ys, lc=None, color=None, **kwargs) -> matplotlib.collections.LineCollection:
    """Plot lines with different colorings

    Parameters
    ----------
    xs : iterable container of x coordinates
    ys : iterable container of y coordinates
    lc : line collection
    color : iterable container of numbers mapped to colormap
    ax (optional): Axes to plot on.
    kwargs (optional): passed to LineCollection

    Returns
    -------
    Notes:
        len(xs) == len(ys) == len(c) is the number of line segments
        len(xs[i]) == len(ys[i]) is the number of points for each line (indexed by i)

    """

    # create LineCollection
    sz = len(xs) - 1

    segments = np.zeros((sz, 2, 2))
    # Each line segment has start and end point. End of each segment is equivalent to start of next
    segments[:, 0, 0] = xs[0:sz]
    segments[:, 0, 1] = ys[0:sz]
    segments[:, 1, 0] = xs[1:sz+1]
    segments[:, 1, 1] = ys[1:sz+1]

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


def get_next_speed(s, inc=1):
    s = s + inc
    if s >= len(ANIMATE_SPEEDS):
        s = len(ANIMATE_SPEEDS) - 1
    if s <= 0:
        s = 1
    return s


# Format a single floating point number to have 3 digits after decimal
def fmt(n):
    return '{:+2.3f}'.format(n)


# Format a 2D vector as text inside brackets
def fmt_complex(r):
    return fmt(r.real) + ', ' + fmt(r.imag) + 'j'


class GraphicsObjects:

    pltLineCollection = None
    th1 = None
    line_colors: np.ndarray = None
    localShowArrows = None
    line_length = None

    flagAxisScaleIndex = 1  # Starting scale. Can be changed if user zooms in or out
    flagInputRangeIndex = 0
    flagAnimateSpeedIndex = 0
    flagSavedAnimateSpeedIndex = 0  # When we pause, index is saved here so we can unpause and resume original speed

    quit_flag = False  # When true, program will exit
    flagRecalculate = False  # When true, will redraw existing graph items
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

    flagMouseDown = False

    arrowOutput = []

    line_min = None
    line_max = None

    bg1 = None

    def __init__(self, plot_domain, show_buttons=True):
        self.b_reset = None
        self.b_show_arrows = None
        self.b_zoom = None
        self.b_dots = None
        self.b_quit = None
        self.b_output = None
        self.b_1_vs_2 = None
        self.b_animate = None

        self.localPlotDomain = plot_domain

        print("Available backends:")
        print(mpl.rcsetup.all_backends)
        mpl.use('TkAgg')   # Qt5Agg might also be available, but it is MUCH slower. Force TkAgg, which plots much faster
        self.backend = mpl.get_backend()
        print("Matplotlib backend is: " + self.backend)  # Returns Qt5Agg after installing Qt5, otherwise TkAgg

        # Create figure
        fig = plt.figure()  # fig_size=(MAIN_WINDOW_SIZE, MAIN_WINDOW_SIZE))

        window = plt.get_current_fig_manager().window
        dpi = fig.dpi

        if self.backend == "Qt5Agg":
            # Hack to get screen size. Make full-screen temporary window, get size, then later set "real" size
            window.showMaximized()  # Make window full screen
            plt.pause(.001)  # Draw items to screen so we can get size
            screen_x, screen_y = fig.get_size_inches() * fig.dpi  # size in pixels
        else:
            # window.state('zoomed')  # Make window full screen, for TkAgg
            screen_x, screen_y = window.wm_maxsize()  # Get full screen coordinates for TkAgg. Doesn't work with Qt5Agg
            # Shrink 5% to avoid hitting task bar
            screen_y *= 0.95

        screen_y = screen_y - 50  # Subtract a small amount or else the toolbar at bottom will mess things up.

        # Create window with aspect ratio 1.5:1.0
        fig.set_size_inches(screen_y * 1.5 / dpi, screen_y / dpi)  # Convert from pixels to inches by dividing by dpi

        if self.backend == "Qt5Agg":
            plt.pause(.001)  # Draw to screen so that new window size takes effect
#        fig.tight_layout() # This doesn't work, somehow

        canvas = fig.canvas
        canvas.mpl_connect('button_press_event', self.on_mouse_press)
        canvas.mpl_connect('button_release_event', self.on_mouse_release)
        canvas.mpl_connect('key_press_event', self.on_keypress)
        canvas.mpl_connect('motion_notify_event', self.on_mouse_move)
        self.fig = fig
        self.canvas = canvas
        self.ax = plt.gca()

        # Make smaller margins
        plt.subplots_adjust(left=0.04, right=0.99, top=0.99, bottom=0.05)

        canvas.draw()
        self.bg_blank = canvas.copy_from_bbox(self.ax.bbox)

        self.Add_axes(plot_domain)

        self.textObj = TextObjects(self.canvas, self.fig)

        if show_buttons:
            self.make_buttons()

        # Need this or else current axis will be the last button drawn, instead of main plot
        plt.sca(self.ax)

    # Create initial graphics, then save backgrounds
    def Add_axes(self, plot_domain):

        self.localPlotDomain = plot_domain

        plt.sca(self.ax)
        plt.show(block=False)

        self.set_zoom()  # Need to do this after sca, or else it won't take effect

        self.ax.axhline(y=0, color='gray')
        self.ax.axvline(x=0, color='gray')

        # Make smaller margins
        #    self.gObjects.ax.margins(x=0.01) # This doesn't work. Use subplots_adjust() instead

    def make_buttons(self):

        # Create one new axis for each new button
        ax_animate = plt.axes([BUTTON_X, BUTTON_Y, BUTTON_WIDTH, BUTTON_HEIGHT])
        self.b_animate = Button(ax_animate, 'Toggle animate')
        self.b_animate.on_clicked(self.do_pause_unpause_animate)

        ax_1_vs_2 = plt.axes([BUTTON_X + BUTTON_SPACING_X, BUTTON_Y, BUTTON_WIDTH, BUTTON_HEIGHT])
        self.b_1_vs_2 = Button(ax_1_vs_2, 'Bidirectional')
        self.b_1_vs_2.on_clicked(self.do_1_vs_2)

        ax_output = plt.axes([BUTTON_X + BUTTON_SPACING_X * 2, BUTTON_Y, BUTTON_WIDTH, BUTTON_HEIGHT])
        self.b_output = Button(ax_output, 'Expand input')
        self.b_output.on_clicked(self.do_expand_input)

        ax_reset = plt.axes([BUTTON_X + BUTTON_SPACING_X * 3, BUTTON_Y, BUTTON_WIDTH, BUTTON_HEIGHT])
        self.b_reset = Button(ax_reset, 'Reset animation (r)')
        self.b_reset.on_clicked(self.do_reset_animation)

        ax_show_arrows = plt.axes([BUTTON_X + BUTTON_SPACING_X * 4, BUTTON_Y, BUTTON_WIDTH, BUTTON_HEIGHT])
        self.b_show_arrows = Button(ax_show_arrows, 'Toggle arrows')
        self.b_show_arrows.on_clicked(self.do_toggle_arrows)

        ax_button_zoom = plt.axes([BUTTON_X + BUTTON_SPACING_X * 5, BUTTON_Y, BUTTON_WIDTH, BUTTON_HEIGHT])
        self.b_zoom = Button(ax_button_zoom, 'Zoom out (z)')
        self.b_zoom.on_clicked(self.do_zoom_out)

        ax_button_dots = plt.axes([BUTTON_X + BUTTON_SPACING_X * 6, BUTTON_Y, BUTTON_WIDTH, BUTTON_HEIGHT])
        self.b_dots = Button(ax_button_dots, 'Toggle center')
        self.b_dots.on_clicked(self.do_toggle_center)

        ax_quit = plt.axes([BUTTON_X + BUTTON_SPACING_X * 7, BUTTON_Y, BUTTON_WIDTH, BUTTON_HEIGHT])
        self.b_quit = Button(ax_quit, 'Quit (x)')
        self.b_quit.on_clicked(self.do_quit)

    # Create input and output rainbow lines
    def InitAnimation(self):

        if self.pltLineCollection is not None:
            self.pltLineCollection.remove()

        # pre-generate entire line
        th2b = rm.riemann(self.th1)
        if self.flagAxisScaleIndex > 0:
            # Make lighter colors when scaling is largest
            p1: np.ndarray = ([1.0, 1.0, 1.0, 1.0] - self.line_colors) / 4
            colors2 = -p1 + [1.0, 1.0, 1.0, 1.0]
            line2 = multiline(np.real(th2b), np.imag(th2b), color=colors2, linewidths=0.5)
        else:
            line2 = multiline(np.real(th2b), np.imag(th2b), color=self.line_colors, linewidths=0.5)
        self.pltLineCollection = self.ax.add_collection(line2)

        # Hide any items that are not static
        self.setArrowVisible(False)
        self.textObj.text_clear()

        # Save background
        self.Save_background(self.canvas)

        self.setArrowVisible(True)

        return self.pltLineCollection, th2b, line2

    def Save_background(self, canvas):

        saveIndex = self.flagAxisScaleIndex
        self.bg1 = []
        # Save background bitmaps for each scale
        for self.flagAxisScaleIndex in range(0, len(AXIS_SCALES)):
            self.set_zoom()
            canvas.draw()
            self.bg1.append(canvas.copy_from_bbox(self.ax.bbox))

        # Restore original scale
        self.flagAxisScaleIndex = saveIndex
        self.set_zoom()
        canvas.draw()

    def CheckArrows(self):
        if self.localShowArrows != self.flagShowArrows:
            # Just switched on or off. Hide or show arrows.
            self.localShowArrows = self.flagShowArrows
            self.setArrowVisible(self.localShowArrows)

    def MakeLine1(self, plot_domain, input_scale, RESOLUTION, LINE_X) -> tuple[np.ndarray, np.ndarray]:
        V_LINE_HEIGHT = plot_domain * input_scale
        V_LINE_SEGMENTS = int(RESOLUTION * V_LINE_HEIGHT)

        # Rainbow colors
        colors = pl.cm.jet(np.linspace(0, 1, V_LINE_SEGMENTS))

        max_x = LINE_X + 1j * V_LINE_HEIGHT

        if self.flagBidirectional:
            min_x = LINE_X - 1j * V_LINE_HEIGHT
        else:
            min_x = LINE_X + 0j

        lns = np.linspace(min_x, max_x, V_LINE_SEGMENTS)

        self.line_min = min_x
        self.line_max = max_x
        self.line_length = np.abs(max_x - min_x)

        return lns, colors

    def GetPoint(self, position):
        return self.line_min + (self.line_max - self.line_min) * position / self.line_length

    # Mouse button press. Use this to start moving vector1 or vector2 in top-right plot
    def on_mouse_press(self, event):
        if event.inaxes != self.ax:
            return
        if event.button == mpl.backend_bases.MouseButton.LEFT:
            self.flagMouseDown = True
            self.flagChangeMatrix = True
            self.whichRowToAdjust = 0
        if event.button == mpl.backend_bases.MouseButton.RIGHT:
            self.flagMouseDown = True
            self.flagChangeMatrix = True
            self.whichRowToAdjust = 1
        self.flagX = event.xdata
        self.flagY = event.ydata

    def on_mouse_release(self, event):
        if event.inaxes != self.ax:
            return
        self.flagMouseDown = False
        self.flagChangeMatrix = False

    def on_mouse_move(self, event):
        if event.inaxes != self.ax:
            return
        if self.flagMouseDown:
            self.flagX = event.xdata
            self.flagY = event.ydata
            self.flagChangeMatrix = True

    # Set plot axis scaling to zoom in/out
    def set_zoom(self):
        scaleTmp = self.GetAxisScale()
        pltHeight = self.localPlotDomain * scaleTmp
        pltWidth = pltHeight * 1.5
        if scaleTmp > 1:
            plt.xlim([-pltWidth * .75, pltWidth * 1.25])
            plt.ylim([-pltHeight, pltHeight])
        elif scaleTmp <= 0.1 or self.flagCenterAtTarget:
            plt.xlim([-pltWidth, pltWidth])
            plt.ylim([-pltHeight, pltHeight])
        elif scaleTmp <= 0.3:
            plt.xlim([-pltWidth * 0.1, pltWidth * 1.9])
            plt.ylim([-pltHeight, pltHeight])
        else:
            plt.xlim([-pltWidth * .667, pltWidth * 1.333])
            plt.ylim([-pltHeight, pltHeight])

    def Restore_background(self, canvas):
        if self.flagCenterAtTarget:
            canvas.restore_region(self.bg_blank)
        else:
            canvas.restore_region(self.bg1[self.flagAxisScaleIndex])  # Restores static elements and erases background

    def do_1_vs_2(self, _event=None):
        self.flagBidirectional = not self.flagBidirectional  # Toggle between 1 and 2 rows
        self.flagRecalculate = True

    # Pause and unpause animation
    def do_pause_unpause_animate(self, _event=None):
        if self.flagAnimateSpeedIndex == 0:
            self.flagAnimateSpeedIndex = self.flagSavedAnimateSpeedIndex if (self.flagSavedAnimateSpeedIndex > 0) else 1
        else:
            self.flagSavedAnimateSpeedIndex = self.flagAnimateSpeedIndex
            self.flagAnimateSpeedIndex = 0

    def do_speed_up(self, _event=None):
        if self.flagAnimateSpeedIndex == 0:
            self.flagSavedAnimateSpeedIndex = get_next_speed(self.flagSavedAnimateSpeedIndex)
        else:
            self.flagAnimateSpeedIndex = get_next_speed(self.flagAnimateSpeedIndex)

    def do_slow_down(self, _event=None):
        if self.flagAnimateSpeedIndex == 0:
            self.flagSavedAnimateSpeedIndex = get_next_speed(self.flagSavedAnimateSpeedIndex, -1)
        else:
            self.flagAnimateSpeedIndex = get_next_speed(self.flagAnimateSpeedIndex, -1)

    def do_reset_animation(self, _event=None):
        self.flagRestartAnimation = True
        self.flagChangeMatrix = True

    def do_expand_input(self, _event=None):
        self.flagRecalculate = True
        self.flagInputRangeIndex = self.flagInputRangeIndex + 1
        if self.flagInputRangeIndex >= len(INPUT_RANGES):
            self.flagInputRangeIndex = 0
        self.flagRedrawAxes = True

    def do_toggle_arrows(self, _event=None):
        self.flagShowArrows = not self.flagShowArrows

    def do_zoom_out(self, _event=None):
        self.flagAxisScaleIndex = self.flagAxisScaleIndex + 1
        if self.flagAxisScaleIndex >= len(AXIS_SCALES):
            self.flagAxisScaleIndex = 0
        self.flagRedrawAxes = True

    def do_toggle_side_dots(self, _event=None):
        self.flagSideDots = not self.flagSideDots

    def do_toggle_center(self, _event=None):
        self.flagCenterAtTarget = not self.flagCenterAtTarget
        self.flagRedrawAxes = True

    def GetAxisScale(self):
        return AXIS_SCALES[self.flagAxisScaleIndex]

    def get_animation_speed(self):
        return ANIMATE_SPEEDS[self.flagAnimateSpeedIndex]

    # function to turn arrow visibility on and off
    def setArrowVisible(self, b):
        for i in range(0, len(self.arrowOutput)):
            self.arrowOutput[i].set_visible(b)

    def do_quit(self, _event=None):
        self.quit_flag = True

    # Keyboard press
    def on_keypress(self, event):
        if event.key == "x":
            self.do_quit()
        elif event.key == " ":
            self.do_pause_unpause_animate()
        elif event.key == "r":
            self.do_reset_animation()
        elif event.key == "z":
            self.do_zoom_out()
        elif event.key == "v":
            self.do_speed_up()
        elif event.key == "up":
            self.do_speed_up()
        elif event.key == "down":
            self.do_slow_down()
        elif event.key == "left":
            self.flagX = self.flagX - 0.001
            self.flagChangeMatrix = True
            self.whichRowToAdjust = 0
        elif event.key == "right":
            self.flagX = self.flagX + 0.001
            self.flagChangeMatrix = True
            self.whichRowToAdjust = 0


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
