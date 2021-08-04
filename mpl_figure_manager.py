
import matplotlib.pyplot as plt


class FigureManager:

    def __init__(self, root, is_android):

        self.control_window = root
        self.is_android = is_android

        self.screen_x_pixels = root.winfo_screenwidth()
        self.screen_y_pixels = root.winfo_screenheight()
        self.dpi = root.winfo_fpixels('1i')
        print("Screen y resolution: " + str(self.screen_y_pixels))
        print("Screen dpi: " + str(self.dpi))
        # Declare variable that will hold most recent figure
        self.fig_plot: plt.Figure = None

    def make_plot_fig(self, aspect=1.0, y_size=0.95):
        # size parameters are as fraction of screen HEIGHT

        if aspect < 0.3:
            # Avoid extremely skinny windows, as they will flicker annoyingly when
            # matplotlib attempts to show x,y text when mouse hovers over plot area.
            aspect = 0.3

        if aspect > 16/9:
            # Avoid extremely wide windows
            aspect = 16/9

        self.fig_plot = plt.figure()

        height_pixels = int(self.screen_y_pixels * y_size)
        width_pixels = int(self.screen_y_pixels * aspect * y_size)
        x_inches = width_pixels / self.fig_plot.dpi
        y_inches = height_pixels / self.fig_plot.dpi
        print("Figure dpi is %1.2f" % self.fig_plot.dpi)
        print("Creating window with size (in inches) %1.2f x %1.2f" % (x_inches, y_inches))

        # Make large square window for main plots. For some reason, on Android window comes out much
        # smaller than expected.
        self.fig_plot.set_size_inches(x_inches, y_inches)

        # Make smaller margins
        plt.tight_layout()
        # Make even smaller margins
        plt.subplots_adjust(left=0.05, right=0.99, top=0.95, bottom=0.05)

        # Move plot window to the right to avoid overlapping buttons
        if not self.is_android:
            menu_width_pixels = self.control_window.winfo_width() + self.control_window.winfo_x()
            self.fig_plot.canvas.manager.window.wm_geometry("+%d+%d" % (menu_width_pixels + 5, 0))

        #        self.canvas_plot.draw()   # Need this or figure won't show
        plt.pause(0.001)  # Need this or figure won't show

        return self.fig_plot

