
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.widgets as widgets
from matplotlib.widgets import Button


class ButtonManager:

    def __init__(self, num_rows=5, num_cols=2):

        self.BUTTON_Y_GAP = 0.01                  # Space (as fraction of screen) between buttons

        # Button dimensions are all expressed as fraction of window size
        self.BUTTON_Y_START = 0.9
        self.BUTTON_HEIGHT = 0.8 / num_rows
        self.BUTTON_WIDTH = 0.4
        self.BUTTON_X_GAP = 0.05
        self.BUTTON_X_START = (1 - self.BUTTON_WIDTH * num_cols - self.BUTTON_X_GAP * (num_cols-1)) / 2

        self.buttonX = self.BUTTON_X_START
        self.buttonY = self.BUTTON_Y_START

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

        self.screen_y_pixels = screen_y

        # Buttons on plot window cause annoying flicker whenever mouse moves over button
        # (even if not clicked). Solve this by putting buttons on their own window

        # menu_height = self.screen_y / 2
        self.menu_height_pixels = 700
        self.menu_width_pixels = 700 * 0.67

        self.fig2.set_size_inches(self.menu_width_pixels / self.dpi, self.menu_height_pixels / self.dpi)

        self.canvas2 = self.fig2.canvas
        # Put button window at top left of screen
        self.move_window(self.canvas2, 25, 25)

    def make_plot_fig(self, x_size = 1.0, y_size = 1.0):
        # size parameters are as fraction of screen HEIGHT
#        mpl.rcParams['toolbar'] = 'toolbar2'

        self.fig1 = plt.figure()
        self.set_fig1_size(self.fig1, x_size, y_size)
        self.canvas1 = self.fig1.canvas

        # Make smaller margins
        plt.tight_layout()
        # Make even smaller margins
        plt.subplots_adjust(left=0.05, right=0.99, top=0.95, bottom=0.05)

        screen_y_adj = int(self.screen_y_pixels * .95)  # Reduce height about 5% so we don't overlap windows taskbar

        # Move plot window to the right to avoid overlapping buttons
        self.move_window(self.canvas1, self.menu_width_pixels + 25, 0)

        return self.fig1

    def set_fig1_size(self, fig, xsize = 1.0, ysize = 1.0):
        # size units are in fractions of screen HEIGHT.
        screen_y_adj = int(self.screen_y_pixels * .95 * ysize)  # Reduce height about 5% so we don't overlap windows taskbar

        # Make large square window for main plots
        fig.set_size_inches(screen_y_adj * xsize / self.dpi, screen_y_adj * ysize / self.dpi)

    def reset_button_coord(self, col=0):
        # Call this whenever you start a new column of buttons
        self.buttonX = self.BUTTON_X_START + (self.BUTTON_WIDTH + self.BUTTON_X_GAP)*col
        self.buttonY = self.BUTTON_Y_START

    def add_blank(self, col=0):
        # Add blank space
        self.increment_row()

    def add_id_button(self, text, _id):

        return Button2(self.next_button_axis(), text, _id)

    def add_standard_button(self, text):

        ax = Button(self.next_button_axis(), text)
        return ax

    def add_checkbox(self, text):

        # Checkbox is 3 rows high, and width is same as height
        button_height_pixels = self.BUTTON_HEIGHT * 3 * self.menu_height_pixels
        button_width_fraction =button_height_pixels / self.menu_width_pixels
        ax = widgets.CheckButtons(self.next_button_axis(num_rows=3, button_width=button_width_fraction), text)
        return ax

    def add_textbox(self, label, _id):

        ax = widgets.TextBox(self.next_button_axis(), label, _id)
        return ax

    def add_id_slider(self, label, _id, valmin=0, valmax=1, **kwargs):

        ax = Slider2(self.next_button_axis(), label, _id, valmin, valmax, **kwargs)
        return ax

    def next_button_axis(self, num_rows=1, button_width=0):
        # num_rows indicates how many rows this axis should occupy
        # button_width is width of button as fraction of window width. If 0, then will use self.BUTTON_WIDTH
        self.increment_row(num_rows-1)
        if button_width == 0:
            button_width = self.BUTTON_WIDTH
        ax = plt.axes([self.buttonX, self.buttonY,
                       button_width, self.BUTTON_HEIGHT * num_rows + self.BUTTON_Y_GAP * (num_rows - 1)])
        self.increment_row()
        return ax

    def increment_row(self, height=1):
        # Increment coordinates to next button in vertical column
        self.buttonY = self.buttonY - (self.BUTTON_HEIGHT + self.BUTTON_Y_GAP) * height

    def move_window(self, canvas, x_pixels, y_pixels): # x and y are in pixels

        if self.backend == "Qt5Agg":
            geom = canvas.manager.window.geometry()
            x1,y1,dx,dy = geom.getRect()
            canvas.manager.window.setGeometry(x_pixels, y_pixels + 50, dx, dy)
        elif self.backend == "TkAgg":
            canvas.manager.window.wm_geometry("+%d+%d" % (x_pixels, y_pixels))
        else:
            print("Unsupported backend " + self.backend)


class Button2(Button):
    #
    # Implements Button with ID value (useful when there is a long list of buttons)
    #

    def __init__(self, ax, label, _id):
        """
        """
        Button.__init__(self, ax, label, color='0.85', hovercolor='0.95')

        self.id = _id


class Slider2(widgets.Slider):
    #
    # Implements Button with ID value (useful when there is a long list of buttons)
    #

    def __init__(self, ax, label, _id, valmin, valmax, **kwargs):
        """
        """
        widgets.Slider.__init__(self, ax=ax, label=label, valmin=valmin, valmax=valmax, **kwargs)

        self.id = _id


class TextBox2(widgets.TextBox):
    #
    # Implements Button with ID value (useful when there is a long list of buttons)
    #

    def __init__(self, ax, label, _id):
        """
        """
        widgets.Slider.__init__(self, ax, label)

        self.id = _id
