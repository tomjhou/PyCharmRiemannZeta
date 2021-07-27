#
# Calculates Riemann zeta function for values along vertical line. Line is then animated as sweeping from left to
# right. In the middle of the animation, line will be values for which Re(x) = 0.5
#
# Uses Euler transform calculation of Dirichlet eta function, which converges faster than Riemann sum.
#
# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.pylab as pl
import numpy as np
from enum import Enum

# My files
import riemann_graphics as mg
import riemann_math as rm
import global_vars as gv
import matplotlib.lines as mlines

OUTPUT_SCALE = 1    # Scales transformed lines so they are more easily visible
INPUT_SCALE_BASE = 30       # Size of initial line. Will be scaled to fit on plot. Larger numbers allow input to range more widely
INPUT_SCALE_EXPANDED = 100  # Expanded input range
RIEMANN_ITER_LIMIT = 275
INPUT_LINE_RESOLUTION = 40     # Number of points per unit on input line
PLOT_DOMAIN = np.pi / 2
OUTPUT_ARROW_SCALE = 0.03
INPUT_ARROW_SCALE = 0.03
CIRCLE_PATCH_SCALE = 0.01
INPUT_SCALE = INPUT_SCALE_EXPANDED if mg.flagExpandInput else INPUT_SCALE_BASE

InputStyle = mg.PatchStyle.ARROW_CIRCLE

if __name__ == '__main__':
    print('\nRiemann-Zeta demo instructions:\n\n' +
          'Left mouse button drags red input vector\n' +
          'Right mouse button drags vertical line\n' +
          '"x" to exit')

# Create figure graphics
gObjects = mg.GraphicsObjects()
gv.ax1 = gObjects.ax
canvas = gObjects.canvas

mg.flagBidirectional = True
mg.flagAxisScale = 0.1
ANIMATION_RANGE_BASE = 0.5
ANIMATION_STEP_BASE = 0.01

# Create static objects, and save background bitmaps so we don't have to redraw them over and over.
mg.Add_axes(canvas, PLOT_DOMAIN)

ANIMATION_RANGE = ANIMATION_RANGE_BASE * mg.GetAxisScale()
ANIMATION_STEP = ANIMATION_STEP_BASE * mg.GetAxisScale()

LINE_X = 0.5 - ANIMATION_RANGE / 2
th1, line_colors = mg.MakeLine(PLOT_DOMAIN, INPUT_SCALE, INPUT_LINE_RESOLUTION, LINE_X)
tmp, numArrowsPlusOne = rm.Riemann(2, True) # Last argument returns required array size
# Create new outArray that will be used to generate arrows. Has one more element than arrow
rm.outArray = np.zeros(numArrowsPlusOne, dtype=complex)

# Create input rainbow line
line1 = mg.multiline(np.real(th1), np.imag(th1), color=line_colors, linewidths=0.5)
gv.ax1.add_collection(line1)

# Create output rainbow line
th2 = rm.Riemann(th1)
th2b = [x * OUTPUT_SCALE for x in th2]
line2 = mg.multiline(np.real(th2b), np.imag(th2b), color=line_colors, linewidths=0.5)
pltLineCollection = gv.ax1.add_collection(line2)

canvas.flush_events()

# Need to call this the first time or else objects won't draw later
plt.pause(0.01)

localAnimate = mg.flagAnimate  # Need local copy of this flag so we can detect when it changes state
# Update text
gObjects.textObj.update_input_real(LINE_X)

while mg.quit_flag == 0:

    # Render to screen
    if not mg.flagAnimate:
        plt.pause(0.001)  # This is necessary to redraw line segments

    # This gives better performance after initial draw
    canvas.flush_events()

    while not mg.quit_flag:

        if localAnimate != mg.flagAnimate:
            # Animation has just changed.
            localAnimate = mg.flagAnimate

        if mg.flagRedrawAxes:
            mg.flagRedrawAxes = False
            mg.set_zoom()
            ANIMATION_RANGE = ANIMATION_RANGE_BASE * mg.GetAxisScale()
            ANIMATION_STEP = ANIMATION_STEP_BASE * mg.GetAxisScale()
            LINE_X = 0.5 - ANIMATION_STEP * 2

        if mg.flagChangeMatrix:
            LINE_X = mg.flagX

        if localAnimate or mg.flagChangeMatrix:
            mg.flagRecalculate = False
            LINE_X = LINE_X + ANIMATION_STEP
            if LINE_X > 0.5 + ANIMATION_RANGE / 2:
                LINE_X = 0.5 - ANIMATION_RANGE / 2
            INPUT_SCALE = INPUT_SCALE_EXPANDED if mg.flagExpandInput else INPUT_SCALE_BASE

            th1, line_colors = mg.MakeLine(PLOT_DOMAIN, INPUT_SCALE, INPUT_LINE_RESOLUTION, LINE_X)
            th2 = rm.Riemann(th1)
            th2b = [x * OUTPUT_SCALE for x in th2]

            mg.multiline(np.real(th1), np.imag(th1), line1, color=line_colors)
            mg.multiline(np.real(th2b), np.imag(th2b), line2, color=line_colors)
            # Update text
            gObjects.textObj.update_input_real(LINE_X)
            break

        canvas.flush_events()
        plt.pause(0.01)

    if mg.quit_flag:
        break

    if mg.flagAnimate:
        plt.pause(0.001)

plt.close()
