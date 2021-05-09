#
# Pre-calculates Riemann Zeta transform on vertical line, then allows user to browse along pre-calculated line
#
# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.pylab as pl
#from matplotlib.lines import Line2D
import numpy as np

# My files
#import MatrixDemoMath as mm
import RiemannGraphics as mg
import global_vars as gv
import matplotlib.lines as mlines

LINE_X = mg.flagX
OUTPUT_SCALE = 1    # Scales transformed lines so they are more easily visible
INPUT_SCALE_BASE = 10       # Size of initial line. Will be scaled to fit on plot. Larger numbers allow input to range more widely
INPUT_SCALE_EXPANDED = 30  # Expanded input range
ARROW_SCALE = 1
RIEMANN_ITER_LIMIT = 500
RESOLUTION = 12     # Number of points per unit
PLOT_SIZE = np.pi * 1.5

INPUT_SCALE = INPUT_SCALE_EXPANDED if mg.flagExpandInput else INPUT_SCALE_BASE

inputVal = LINE_X+0j
axisLimit = PLOT_SIZE  # Coordinate limits for x-y plots

global bg1

# Create initial x-y, text, and bar plots, then save backgrounds
def create_initial_graphics():
    global bg1

#    plt.subplot(111)
    plt.sca(gv.ax1)

    plt.xlim([-axisLimit*0.5, axisLimit*2.5])
    plt.ylim([-axisLimit, axisLimit])
    plt.show(block=False)

    # Need this so that background bitmaps will be up to date
    canvas.draw()

    # Save background bitmaps
    bg1 = canvas.copy_from_bbox(gv.ax1.bbox)


def Riemann(v):
    if np.size(v) > 1:
        return [Riemann(x)[0] for x in v]

    cumSum = 0 + 0j
    outArray[0] = cumSum
    plotNum = 1
    for i in range(1, RIEMANN_ITER_LIMIT):
        cumSum = cumSum + (1 / (i ** v))
        if i < 30:
            outArray[plotNum] = cumSum
            plotNum = plotNum + 1
        elif i < 300:
            if np.mod(i,10) == 0:
                outArray[plotNum] = cumSum
                plotNum = plotNum + 1
        elif np.mod(i,100) == 0:
            outArray[plotNum] = cumSum
            plotNum = plotNum + 1

#    if len(outArray) > plotNum:
#        outArray = outArray[:plotNum]

    return cumSum, plotNum

def MakeLine():
    V_LINE_HEIGHT = PLOT_SIZE * INPUT_SCALE
    V_LINE_SEGMENTS = int(RESOLUTION * V_LINE_HEIGHT)

    # Rainbow colors
    colors = pl.cm.jet(np.linspace(0, 1, V_LINE_SEGMENTS))

    if mg.flagBidirectional:
        lns = np.linspace(LINE_X - 1j * V_LINE_HEIGHT, LINE_X + 1j * V_LINE_HEIGHT, V_LINE_SEGMENTS)
    else:
        lns = np.linspace(LINE_X + 0j, LINE_X + 1j * V_LINE_HEIGHT, V_LINE_SEGMENTS)

    return lns, colors

if __name__ == '__main__':
    print('\nReimann-Zeta demo instructions:\n\n' +
          'Left mouse button drags red input vector\n' +
          'Right mouse button drags vertical line\n' +
          '"x" to exit')

# Create figure
gObjects = mg.GraphicsObjects()
gv.ax1 = gObjects.ax

canvas = gObjects.canvas

# Create static objects, and save background bitmaps so we don't have to redraw them over and over.
# This greatly speeds up animation
create_initial_graphics()

currentStep = 0
cycles = 0

# Create red "input" arrow
arrowInput = mpatches.FancyArrowPatch((0, 0), (0, 0), color=mg.INPUT_VECTOR_COLOR, mutation_scale=ARROW_SCALE)
gv.ax1.add_patch(arrowInput)

LINE_X = mg.flagX
th1, linecolors = MakeLine()
th1b = [x / INPUT_SCALE for x in th1]
outArray = np.zeros(9999, dtype=complex) # Make temporary outArray that is too large
tmp, numArrowsPlusOne = Riemann(1.0) # Use single value to get number of arrows

# Create new outArray that will be used to generate arrows. Has one more element than arrow
outArray = np.zeros(numArrowsPlusOne, dtype=complex)

# Create black "output" arrows
arrowOutput = [mpatches.FancyArrowPatch((0, 0), (0, 0), color=mg.OUTPUT_VECTOR_COLOR, mutation_scale=ARROW_SCALE)]
gv.ax1.add_patch(arrowOutput[0])
for i in range(1, numArrowsPlusOne - 1):
    arrowOutput.append(mpatches.FancyArrowPatch((0, 0), (0, 0), color=mg.OUTPUT_VECTOR_COLOR, mutation_scale=ARROW_SCALE))
    gv.ax1.add_patch(arrowOutput[i])


# Create input and output lines
line1 = mg.multiline(np.real(th1b), np.imag(th1b), color=linecolors, linewidths = 0.5)
gv.ax1.add_collection(line1)

#th2 = Riemann(th1) # If we pass array to Riemann, we no longer get # of arrows
#th2b = [x * OUTPUT_SCALE for x in th2]

th2b = Riemann(th1[0:2]) # First two points
line2 = mg.multiline(np.real(th2b), np.imag(th2b), color=linecolors, linewidths = 0.5)
gv.ax1.add_collection(line2)
currentPoint = 2
for x in range(2,len(th1)):
    t = th1[x]
    th2b.append(Riemann(t)[0])
    line2 = mg.multiline(np.real(th2b), np.imag(th2b), line2, color=linecolors, linewidths = 0.5)
    if np.mod(x,20) == 0:
        canvas.flush_events()
        plt.pause(0.01)

canvas.flush_events()
plt.pause(0.01)

# Need to call this the first time or else objects won't draw later
plt.pause(0.01)
localRedrawAxes = False
localAnimate = mg.flagAnimate  # Need local copy of this flag so we can detect when it changes state

while mg.quitflag == 0:

    Riemann(inputVal)

    canvas.restore_region(bg1)  # Restores static elements and erases background
    # Draw red arrow for input value
    arrowInput.set_positions((0, 0), tuple([inputVal.real / INPUT_SCALE, inputVal.imag / INPUT_SCALE]))
    gv.ax1.draw_artist(arrowInput)

    # Draw head-to-tail purple arrows
    for i in range(0, len(arrowOutput)):
        s1 = outArray[i] * OUTPUT_SCALE
        s2 = outArray[i+1] * OUTPUT_SCALE
        arrowOutput[i].set_positions(tuple([s1.real, s1.imag]), tuple([s2.real, s2.imag]))
        gv.ax1.draw_artist(arrowOutput[i])

    # Update text
    gObjects.textObj.update_input(inputVal)

    # Render to screen
    if not mg.flagAnimate:
#        canvas.blit(gv.ax1.bbox)    # This will only show the draw_artist() items and background, i.e. head-to-tail arrows
        plt.pause(0.01)  # This is necessary to redraw line segments

    # This gives better performance after initial draw
    if localRedrawAxes:
        plt.pause(0.01)
        localRedrawAxes = False
    else:
        canvas.flush_events()

    while not mg.quitflag:
        if localAnimate and not mg.flagAnimate:
            # Just switched off. Turn off animation, but run one more loop
            localAnimate = False
            break
        if not localAnimate and mg.flagAnimate:
            # Just switched on
            localAnimate = True

        if mg.flagChangeMatrix:
            mg.flagChangeMatrix = False

            if mg.flagRestartAnimation:
                inputVal = 1.1+0j
                LINE_X = 1.1
                mg.flagRestartAnimation = False

            elif mg.flagX is not None and mg.flagY is not None:  # Because x, y won't be valid if mouse went out of bounds

                if mg.whichRowToAdjust == 0:
                    # Left mouse moves single input point
                    inputVal = complex(LINE_X, mg.flagY * INPUT_SCALE)
                    break
                else:
                    # Right mouse moves entire line
                    LINE_X = mg.flagX * INPUT_SCALE
                    inputVal = complex(LINE_X, mg.flagY * INPUT_SCALE)

                    mg.flagRecalc = True

        if mg.flagRecalc:
            mg.flagRecalc = False
            INPUT_SCALE = INPUT_SCALE_EXPANDED if mg.flagExpandInput else INPUT_SCALE_BASE

            th1, linecolors = MakeLine()
            th2 = Riemann(th1)
            th1b = [x / INPUT_SCALE for x in th1]
            th2b = [x * OUTPUT_SCALE for x in th2]

            mg.multiline(np.real(th1b), np.imag(th1b), line1, color=linecolors)
            mg.multiline(np.real(th2b), np.imag(th2b), line2, color=linecolors)
            break

        if mg.flagRedrawAxes:
            mg.flagRedrawAxes = False
            localRedrawAxes = True
            break

        if localAnimate:
            break

        canvas.flush_events()

    if mg.quitflag:
        break

    if mg.flagAnimate:
        plt.pause(0.1)
        inputVal = inputVal + 0.1j
        if np.imag(inputVal) > np.imag(th1[len(th1)-1]):
            inputVal = complex(np.real(inputVal), 0)

plt.close()
