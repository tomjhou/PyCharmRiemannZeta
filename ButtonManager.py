
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.widgets as widgets
from matplotlib.widgets import Button


class ButtonManager:

    def __init__(self, num_buttons=5, _is_vertical = True):

        self.USE_VERTICAL_BUTTON_PANEL = _is_vertical   # Stack buttons vertically instead of horizontally
        self.BUTTON_GAP = 0.01                  # Space (as fraction of screen) between buttons

        # Button dimensions are all expressed as fraction of window size
        if self.USE_VERTICAL_BUTTON_PANEL:
            self.BUTTON_Y_COORD = 0.9
            self.BUTTON_HEIGHT = 0.8 / num_buttons
            self.BUTTON_WIDTH = 0.75
            self.BUTTON_X_START = 0.125
        else:
            self.BUTTON_Y_COORD = 0.1  # 0.95
            self.BUTTON_HEIGHT = 0.8  # 0.03  # Above about 0.1, buttons disappear - presumably they collide with graphs?
            self.BUTTON_WIDTH = 0.9 / num_buttons
            self.BUTTON_X_START = 0.05

        self.buttonX = self.BUTTON_X_START
        self.buttonY = self.BUTTON_Y_COORD

        # Create figure 1 for main plots
        self.fig1 = None
        self.canvas1 = None

        # Create figure 2 with buttons
        self.fig2 = plt.figure(2)
        self.dpi = self.fig2.dpi

        self.backend = mpl.get_backend()

        window = plt.get_current_fig_manager().window

        if self.backend == "Qt5Agg":
            # Hack to get screen size. Temporarily make a full-screen window, get size, then later set "real" size
            window.showMaximized()  # Make fullscreen
            plt.pause(.001)  # Draw items to screen so we can get size
            screen_x, screen_y = self.fig2.get_size_inches() * self.fig2.dpi  # size in pixels
        elif self.backend == "TkAgg":
            screen_x, screen_y = window.wm_maxsize()  # This works for TkAgg, but not Qt5Agg
        else:
            print("Unsupported backend " + self.backend)
            screen_y = 1024

        self.screen_y = screen_y

        # Buttons on plot window cause annoying flicker whenever mouse moves over button
        # (even if not clicked). Solve this by putting buttons on their own window

        screen_y_adj = int(self.screen_y * .95)  # Reduce height about 5% so we don't overlap windows taskbar

        if self.USE_VERTICAL_BUTTON_PANEL:
            self.fig2.set_size_inches(screen_y_adj / self.dpi / 5, screen_y_adj / self.dpi / 2)
        else:
            self.fig2.set_size_inches(screen_y_adj / self.dpi / 2, screen_y_adj / self.dpi / 15)

        self.canvas2 = self.fig2.canvas
        # Put button window at top left of screen
        self.move_window(self.canvas2, 25, 25)

    def make_plot_fig(self, xsize = 1.0, ysize = 1.0):
        # size parameters are as fraction of screen HEIGHT
#        mpl.rcParams['toolbar'] = 'toolbar2'

        self.fig1 = plt.figure()
        self.set_fig1_size(self.fig1, xsize, ysize)
        self.canvas1 = self.fig1.canvas

        # Make smaller margins
        plt.tight_layout()
        # Make even smaller margins
        plt.subplots_adjust(left=0.05, right=0.99, top=0.95, bottom=0.05)

        screen_y_adj = int(self.screen_y * .95)  # Reduce height about 5% so we don't overlap windows taskbar

        if self.USE_VERTICAL_BUTTON_PANEL:
            # Move plot window to the right to avoid overlapping buttons
            self.move_window(self.canvas1, screen_y_adj * .3, 0)
        else:
            # Move plot window to the right to avoid overlapping buttons
            self.move_window(self.canvas1, screen_y_adj * .6, 0)

        return self.fig1

    def set_fig1_size(self, fig, xsize = 1.0, ysize = 1.0):
        # size units are in fractions of screen HEIGHT.
        screen_y_adj = int(self.screen_y * .95 * ysize)  # Reduce height about 5% so we don't overlap windows taskbar

        # Make large square window for main plots
        fig.set_size_inches(screen_y_adj * xsize / self.dpi, screen_y_adj * ysize / self.dpi)

    def add_id_button(self, text, id):

        return Button2(self.next_button_axis(), text, id)

    def add_standard_button(self, text):

        ax = Button(self.next_button_axis(), text)
        return ax

    def add_checkbox(self, text):

        ax = widgets.CheckButtons(self.next_button_axis(), text)
        return ax

    def add_textbox(self, label):

        ax = widgets.TextBox(self.next_button_axis(), label)
        return ax

    def add_slider(self, label):

        ax = widgets.Slider(self.next_button_axis(), label, 0, 1)
        return ax

    def next_button_axis(self):
        # Generate axes for the next button in series (either horizontal or vertical row)
        ax = plt.axes([self.buttonX, self.buttonY, self.BUTTON_WIDTH, self.BUTTON_HEIGHT])

        # Increment coordinates in preparation for next call
        if self.USE_VERTICAL_BUTTON_PANEL:
            self.buttonY = self.buttonY - self.BUTTON_HEIGHT - self.BUTTON_GAP
        else:
            self.buttonX = self.buttonX + self.BUTTON_WIDTH + self.BUTTON_GAP

        return ax

    def move_window(self, canvas, x, y): # x and y are in pixels

        if self.backend == "Qt5Agg":
            geom = canvas.manager.window.geometry()
            x1,y1,dx,dy = geom.getRect()
            canvas.manager.window.setGeometry(x , y + 50, dx, dy)
        elif self.backend == "TkAgg":
            canvas.manager.window.wm_geometry("+%d+%d" % (x, y))
        else:
            print("Unsupported backend " + self.backend)


class Button2(Button):
    #
    # Implements Button with ID value (useful when there is a long list of buttons)
    #

    def __init__(self, ax, label, id):
        """
        """
        Button.__init__(self, ax, label, color='0.85', hovercolor='0.95')

        self.id = id
