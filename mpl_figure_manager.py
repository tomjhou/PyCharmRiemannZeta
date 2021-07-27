
import matplotlib.pyplot as plt
import matplotlib as mpl


class MplFigureManager:

    def __init__(self):

        self.backend = mpl.get_backend()

        fig_tmp = plt.figure()
        self.dpi = fig_tmp.dpi
        window = plt.get_current_fig_manager().window

        if self.backend == "Qt5Agg":
            # Hack to get screen size. Temporarily make a full-screen window, get size, then later set "real" size
            window.showMaximized()  # Make fullscreen
            plt.pause(.001)  # Draw items to screen so we can get size
            screen_x, screen_y = fig_tmp.get_size_inches() * self.dpi  # size in pixels
        elif self.backend == "TkAgg":
            screen_x, screen_y = window.wm_maxsize()  # This works for TkAgg, but not Qt5Agg
        else:
            print("Unsupported backend " + self.backend)
            screen_y = 1024

        # This causes Tkinter mainloop to exit
        plt.close()

        self.screen_y_pixels = screen_y

        self.canvas_buttons = None
        self.menu_width_pixels = 500

        # Create figure 1 for main plots
        self.fig_plot: plt.Figure  # = None
        self.canvas_plot: mpl.backends.backend_tkagg.FigureCanvasTkAgg  # = None

    def make_plot_fig(self, x_size=1.0, y_size=1.0):
        # size parameters are as fraction of screen HEIGHT

        if x_size < 0.3:
            # Avoid extremely skinny windows, as they will flicker annoyingly when
            # matplotlib attempts to show x,y text when mouse hovers over plot area.
            x_size = 0.3

        self.fig_plot = plt.figure()
        self.set_figure_size(self.fig_plot, x_size, y_size)
        self.canvas_plot = self.fig_plot.canvas
#        self.canvas_plot.draw()   # Need this or figure won't show
        plt.pause(0.001)   # Need this or figure won't show

        # Make smaller margins
        plt.tight_layout()
        # Make even smaller margins
        plt.subplots_adjust(left=0.05, right=0.99, top=0.95, bottom=0.05)

        # Move plot window to the right to avoid overlapping buttons
        self.move_window(self.canvas_plot, self.menu_width_pixels + 25, 0)

        return self.fig_plot

    def set_figure_size(self, fig, x_size=1.0, y_size=1.0):
        # size units are in fractions of screen HEIGHT.
        screen_y_adj = int(self.screen_y_pixels * .95 * y_size)  # Reduce height ~5% so we don't overlap windows taskbar

        # Make large square window for main plots
        fig.set_size_inches(screen_y_adj * x_size / self.dpi, screen_y_adj * y_size / self.dpi)

    def move_window(self, canvas, x_pixels, y_pixels):  # x and y are in pixels

        if self.backend == "Qt5Agg":
            geom = canvas.manager.window.geometry()
            x1, y1, dx, dy = geom.getRect()
            canvas.manager.window.setGeometry(x_pixels, y_pixels + 50, dx, dy)
        elif self.backend == "TkAgg":
            canvas.manager.window.wm_geometry("+%d+%d" % (x_pixels, y_pixels))
        else:
            print("Unsupported backend " + self.backend)