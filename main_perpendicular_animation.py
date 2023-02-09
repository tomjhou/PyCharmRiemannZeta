#
# Calculates Riemann zeta function for values along vertical line. Line is then animated as sweeping from left to
# right. In the middle of the animation, line will be values for which Re(x) = 0.5
#
# Uses Euler transform calculation of Dirichlet eta function, which converges faster than Riemann sum.
#
# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import matplotlib.pyplot as plt
import numpy as np

# My files
import riemann_vectors_gui as mg
import riemann_math as rm

OUTPUT_SCALE = 1    # Scales transformed lines so they are more easily visible
INPUT_SCALE_BASE = 30   # Size of initial line. Will scale to fit plot. Larger numbers allow input to range more widely
INPUT_SCALE_EXPANDED = 100  # Expanded input range
RIEMANN_ITER_LIMIT = 275
INPUT_LINE_RESOLUTION = 40     # Number of points per unit on input line
OUTPUT_ARROW_SCALE = 0.03
INPUT_ARROW_SCALE = 0.03
CIRCLE_PATCH_SCALE = 0.01

InputStyle = mg.PatchStyle.ARROW_CIRCLE

INPUT_SCALE = INPUT_SCALE_BASE
PLOT_DOMAIN = np.pi / 2


class PerpendicularAnimate(mg.GraphicsObjects):

    flagExpandInput = False

    def __int__(self, plot_domain):
        # Create figure graphics
        super().__init__(plot_domain, True)

    def make_graph(self):

        global INPUT_SCALE

        if self.flagExpandInput:
            INPUT_SCALE = INPUT_SCALE_EXPANDED

        canvas = self.canvas

        mg.flagBidirectional = True
        mg.flagAxisScale = 0.1
        ANIMATION_RANGE_BASE = 0.5
        ANIMATION_STEP_BASE = 0.01

        ANIMATION_RANGE = ANIMATION_RANGE_BASE * self.GetAxisScale()
        ANIMATION_STEP = ANIMATION_STEP_BASE * self.GetAxisScale()

        LINE_X = 0.5 - ANIMATION_RANGE / 2
        th1, line_colors = self.MakeLine1(PLOT_DOMAIN, INPUT_SCALE, INPUT_LINE_RESOLUTION, LINE_X)
        tmp, numArrowsPlusOne = rm.riemann(2, True)  # Last argument returns required array size
        # Create new outArray that will be used to generate arrows. Has one more element than arrow
        rm.outArray = np.zeros(numArrowsPlusOne, dtype=complex)

        # Create input rainbow line
        line1 = mg.multiline(np.real(th1), np.imag(th1), color=line_colors, linewidths=0.5)
        self.ax.add_collection(line1)

        # Create output rainbow line
        th2 = rm.riemann(th1)
        th2b = [x * OUTPUT_SCALE for x in th2]
        line2 = mg.multiline(np.real(th2b), np.imag(th2b), color=line_colors, linewidths=0.5)
        self.pltLineCollection = self.ax.add_collection(line2)

        canvas.flush_events()
        # Need to call this the first time or else objects won't draw later
        plt.pause(0.01)

        localAnimate = self.flagAnimateSpeedIndex  # Need local copy of this flag so we can detect when it changes state
        # Update text
        self.textObj.update_input_real(LINE_X)

        while self.quit_flag == 0:

            # Render to screen
            if self.flagAnimateSpeedIndex > 0:
                plt.pause(0.001)  # This is necessary to redraw line segments

            # This gives better performance after initial draw
            canvas.flush_events()

            while not self.quit_flag:

                if localAnimate != self.flagAnimateSpeedIndex:
                    # Animation has just changed.
                    localAnimate = self.flagAnimateSpeedIndex

                if self.flagRedrawAxes:
                    self.flagRedrawAxes = False
                    self.set_zoom()
                    ANIMATION_RANGE = ANIMATION_RANGE_BASE * self.GetAxisScale()
                    ANIMATION_STEP = ANIMATION_STEP_BASE * self.GetAxisScale()
                    LINE_X = 0.5 - ANIMATION_STEP * 2

                if self.flagChangeMatrix:
                    LINE_X = self.flagX

                if localAnimate or self.flagChangeMatrix:
                    mg.flagRecalculate = False
                    LINE_X = LINE_X + ANIMATION_STEP
                    if LINE_X > 0.5 + ANIMATION_RANGE / 2:
                        LINE_X = 0.5 - ANIMATION_RANGE / 2
                    INPUT_SCALE = INPUT_SCALE_EXPANDED if self.flagExpandInput else INPUT_SCALE_BASE

                    th1, line_colors = self.MakeLine1(PLOT_DOMAIN, INPUT_SCALE, INPUT_LINE_RESOLUTION, LINE_X)
                    th2 = rm.riemann(th1)
                    th2b = [x * OUTPUT_SCALE for x in th2]

                    mg.multiline(np.real(th1), np.imag(th1), line1, color=line_colors)
                    mg.multiline(np.real(th2b), np.imag(th2b), line2, color=line_colors)
                    # Update text
                    self.textObj.update_input_real(LINE_X)
                    break

                canvas.flush_events()
                plt.pause(0.01)

            if self.quit_flag:
                break

            if self.flagAnimateSpeedIndex > 0:
                plt.pause(0.001)

        plt.close()


if __name__ == '__main__':
    print('\nRiemann-Zeta demo instructions:\n\n' +
          'Left mouse button drags red input vector\n' +
          'Right mouse button drags vertical line\n' +
          '"x" to exit')

    o = PerpendicularAnimate(PLOT_DOMAIN)
    o.make_graph()
