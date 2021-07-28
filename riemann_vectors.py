#
# Calculates Riemann function for values along vertical line in complex plane. Animates result in complex number plane.
#
# Calculation uses Euler's transformation to evaluate Dirichlet eta function, then multiples by 1/(1-2^(1-s))
# to get Riemann zeta. This converges much faster. Calculation is further sped up by precomputing the weighted
# binomial coefficient sum used in Euler transform.
#

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.pylab as pl
import numpy as np
import time

# My files
import riemann_graphics as mg
import riemann_math as rm
import global_vars as gv
import matplotlib.lines as mlines

# Change the default cursor to any valid TK cursor
# To hide it, you'd use the string "none" (or possibly "no" on windows)
old_cursor = []
# If true, generates a nice smooth demo for background screensaver. Very computationally expensive
DEMO_MODE = False
rm.RIEMANN_ITER_LIMIT = 150  # Use higher value to get more precision


def wait_cursor(v=True):
    global old_cursor
    if v:
        old_cursor = gObjects.canvas.toolbar._lastCursor
        gObjects.canvas.toolbar.set_cursor(1)
#        tkagg.cursord[cursors.POINTER] = 'coffee_mug'
    else:
        gObjects.canvas.toolbar.set_cursor(old_cursor)

#        tkagg.cursord[cursors.POINTER] = 'arrow'


if DEMO_MODE:
    # Generate demo that we can use to make background screensaver
    rm.RIEMANN_ITER_LIMIT = 600
    mg.flagAxisScaleIndex = len(mg.AXIS_SCALES) - 1       # Zoomed out
    mg.flagInputRangeIndex = 2      # Max input range
else:
    # Standard mode, less computationally expensive
    rm.RIEMANN_ITER_LIMIT = 150
    mg.flagAxisScaleIndex = len(mg.AXIS_SCALES) - 1       # Zoomed out
    mg.flagInputRangeIndex = 1      # Normal (reduced) input range

rm.precompute_coeffs()

INPUT_LINE_RESOLUTION = 40     # Number of points per unit on input line

PLOT_SIZE = np.pi / 2
OUTPUT_ARROW_SCALE = 0.03
INPUT_ARROW_SCALE = 0.03
CIRCLE_PATCH_SCALE = 0.007
INPUT_SCALE = mg.INPUT_RANGES[mg.flagInputRangeIndex]

InputStyle = mg.PatchStyle.ARROW_CIRCLE

print('\nRiemann-Zeta demo instructions:\n\n' +
      'Left mouse button drags red input vector\n' +
      'Right mouse button drags vertical line\n' +
      'Space bar to pause/unpause animation\n' +
      'Up/down arrows to change animation velocity\n' +
      '"r" to reset animation\n' +
      '"z" to zoom out\n' +
      '"x" to exit')

# Create figure graphics
gObjects = mg.GraphicsObjects(PLOT_SIZE, not DEMO_MODE)
gv.ax1 = gObjects.ax
canvas = gObjects.canvas

LINE_X = mg.flagX


def MakeLine():
    gv.th1, gv.linecolors = mg.MakeLine(PLOT_SIZE, INPUT_SCALE, INPUT_LINE_RESOLUTION, LINE_X)
    th1_scaled = [x / INPUT_SCALE for x in gv.th1]
    # Create input rainbow line
    line1 = mg.multiline(np.real(th1_scaled), np.imag(th1_scaled), color=gv.linecolors, linewidths=0.5)

    return th1_scaled, line1


th1b, line1 = MakeLine()
plt.pause(.001)  # Draws items to screen

tmp, numArrowsPlusOne = rm.riemann(2, True)  # Second argument returns required array size

# Create red "input" arrow
if InputStyle == mg.PatchStyle.ARROW or InputStyle == mg.PatchStyle.ARROW_CIRCLE:
    inputArrowPatch = []  # mpatches.Arrow(0, 0, 0, 0, color=mg.INPUT_VECTOR_COLOR, width=INPUT_PATCH_SCALE)

inputCirclePatch: mpl.patches.Circle
inputCirclePatch2: mpl.patches.Circle
inputCirclePatch3: mpl.patches.Circle

# Create optional side patches
if InputStyle == mg.PatchStyle.CIRCLE or InputStyle == mg.PatchStyle.ARROW_CIRCLE:

    scaleTmp = mg.GetAxisScale()

    inputCirclePatch = mpatches.Circle((0, 0), color=mg.INPUT_VECTOR_COLOR, radius=CIRCLE_PATCH_SCALE * scaleTmp)
    inputCirclePatch.set_edgecolor('none')
    gv.ax1.add_patch(inputCirclePatch)

    inputCirclePatch2 = mpatches.Circle((0, 0), color=mg.INPUT_VECTOR_COLOR, radius=CIRCLE_PATCH_SCALE * scaleTmp)
    inputCirclePatch2.set_edgecolor('none')
    gv.ax1.add_patch(inputCirclePatch2)
    inputCirclePatch3 = mpatches.Circle((0, 0), color=mg.INPUT_VECTOR_COLOR, radius=CIRCLE_PATCH_SCALE * scaleTmp)
    inputCirclePatch3.set_edgecolor('none')
    gv.ax1.add_patch(inputCirclePatch3)

    outputCirclePatch = mpatches.Circle((0, 0), color=(0, 0, 0), radius=CIRCLE_PATCH_SCALE * scaleTmp)
    outputCirclePatch.set_edgecolor('none')
    gv.ax1.add_patch(outputCirclePatch)

    outputCirclePatch2 = mpatches.Circle((0, 0), color=(0, 0, 0), radius=CIRCLE_PATCH_SCALE * scaleTmp)
    outputCirclePatch3 = mpatches.Circle((0, 0), color=(0, 0, 0), radius=CIRCLE_PATCH_SCALE * scaleTmp)
    gv.ax1.add_patch(outputCirclePatch2)
    gv.ax1.add_patch(outputCirclePatch3)

    # Side patches are invisible unless selected later
    inputCirclePatch2.set_visible(False)
    inputCirclePatch3.set_visible(False)
    outputCirclePatch2.set_visible(False)
    outputCirclePatch3.set_visible(False)

# Create new outArray that will be used to generate arrows. Has one more element than arrow
rm.outArray = np.zeros(numArrowsPlusOne, dtype=complex)

# Empty output arrow collection
outputArrowCollection = []
gv.pltLineCollection = None
lineList = []

# Create static objects that do not change
gv.pltLineCollection, th2b, line2 = gObjects.InitAnimation()
currentPoint = 0.0
currentCycle = 0

# Need to call this the first time or else objects won't draw later
plt.pause(0.01)

localAnimateSpeed = mg.get_animation_speed()  # Need local copy of this flag so we can detect when it changes state
gv.localShowArrows = mg.flagShowArrows
localShowSideDots = False

DesiredStepsPerSecond_base = 1

inputVal = gv.th1[0]
outputVal = rm.riemann(inputVal)
mg.arrowOutput = []

wait_cursor(True)
# Create dynamic objects that will change every step
for x in range(0, len(gv.th1)):
    lineTmp = gv.ax1.add_line(mpl.lines.Line2D([0, 0], [0, 0], color=gv.linecolors[x]))
    lineList.append(lineTmp)

arrowColors = pl.cm.jet(np.linspace(0, 4, numArrowsPlusOne - 1))

# Create arrow segments
for i in range(1, numArrowsPlusOne):
    patchTmp = gv.ax1.add_line(mpl.lines.Line2D([0, 0], [0, 0],
                                                color=arrowColors[i-1],
                                                linewidth=5 / np.log2(i + 2)))
    mg.arrowOutput.append(patchTmp)

wait_cursor(False)

plt.pause(0.001)
t1 = time.time()
max_speed = 3
while mg.quit_flag == 0:

    if DEMO_MODE:
        DesiredStepsPerSecond_actual = DesiredStepsPerSecond_base / np.abs(rm.EtaToZetaScale(inputVal))
    else:
        DesiredStepsPerSecond_actual = \
            DesiredStepsPerSecond_base * mg.get_animation_speed() / np.abs(rm.EtaToZetaScale(inputVal))

    inputVal = mg.GetPoint(currentPoint)
    outputVal = rm.riemann(inputVal)

    mg.CheckArrows()

    # Draw head-to-tail purple arrows
    if gv.localShowArrows:
        if mg.GetAxisScale() > 1:
            # Reduce # of arrows when zoomed out, as they will be invisible
            localNumArrows = int(numArrowsPlusOne - 1)
        elif mg.GetAxisScale() == 1:
            localNumArrows = int(numArrowsPlusOne - 1)
        else:
            localNumArrows = numArrowsPlusOne - 1
        for i in range(0, localNumArrows):
            if mg.flagCenterAtTarget:
                s1 = outputVal - rm.outArray[i]
                s2 = outputVal - rm.outArray[i + 1]
                mg.arrowOutput[i].set_data([s1.real, s2.real], [s1.imag, s2.imag])
            else:
                s1 = rm.outArray[i]
                s2 = rm.outArray[i+1]
                mg.arrowOutput[i].set_data([s1.real, s2.real], [s1.imag, s2.imag])

    if mg.flagSideDots != localShowSideDots:
        localShowSideDots = mg.flagSideDots
        inputCirclePatch2.set_visible(localShowSideDots)
        inputCirclePatch3.set_visible(localShowSideDots)
        outputCirclePatch2.set_visible(localShowSideDots)
        outputCirclePatch3.set_visible(localShowSideDots)

    # Draw input/output patches
    if mg.flagCenterAtTarget:
        inputCirclePatch.set_radius(0)
        outputCirclePatch.set_radius(0)
        if inputArrowPatch:
            inputArrowPatch.remove()
            inputArrowPatch = []
    else:
        if InputStyle == mg.PatchStyle.ARROW or InputStyle == mg.PatchStyle.ARROW_CIRCLE:
            # Arrow patch seems to have no way to update, so just remove and replace
            if inputArrowPatch:
                inputArrowPatch.remove()
            inputArrowPatch = mpatches.Arrow(0, 0, inputVal.real / INPUT_SCALE, inputVal.imag / INPUT_SCALE,
                                             color=mg.INPUT_VECTOR_COLOR, width=INPUT_ARROW_SCALE * mg.GetAxisScale())
            gv.ax1.add_patch(inputArrowPatch)
        if (InputStyle == mg.PatchStyle.CIRCLE) or (InputStyle == mg.PatchStyle.ARROW_CIRCLE):
            inputCirclePatch.center = inputVal.real / INPUT_SCALE, inputVal.imag / INPUT_SCALE

            tmpColorIndex = int(currentPoint * len(gv.linecolors) / gv.line_length)
            tmpColor = gv.linecolors[tmpColorIndex]

            inputCirclePatch.set_facecolor(tmpColor)
            inputCirclePatch.set_radius(CIRCLE_PATCH_SCALE * mg.GetAxisScale())
            outputCirclePatch.center = outputVal.real, outputVal.imag
            outputCirclePatch.set_facecolor(tmpColor)
            outputCirclePatch.set_radius(CIRCLE_PATCH_SCALE * mg.GetAxisScale())

            if localShowSideDots:
                inputVal2 = inputVal - 0.2
                inputVal3 = inputVal + 0.2
                inputCirclePatch2.center = inputVal2.real / INPUT_SCALE, inputVal2.imag / INPUT_SCALE
                inputCirclePatch3.center = inputVal3.real / INPUT_SCALE, inputVal3.imag / INPUT_SCALE
                inputCirclePatch2.set_radius(CIRCLE_PATCH_SCALE * mg.GetAxisScale() / 2)
                inputCirclePatch3.set_radius(CIRCLE_PATCH_SCALE * mg.GetAxisScale() / 2)

                outputVal2 = rm.riemann(inputVal2)
                outputVal3 = rm.riemann(inputVal3)

                outputCirclePatch2.center = outputVal2.real, outputVal2.imag
                outputCirclePatch3.center = outputVal3.real, outputVal3.imag
                outputCirclePatch2.set_radius(CIRCLE_PATCH_SCALE * mg.GetAxisScale() / 2)
                outputCirclePatch3.set_radius(CIRCLE_PATCH_SCALE * mg.GetAxisScale() / 2)

    mg.Restore_background(canvas)

    # Update text
    gObjects.textObj.update_input_complex(inputVal)
#    for x in range(0,currentPoint):
#        gv.ax1.draw_artist(lineList[x])
    if gv.localShowArrows:
        for i in range(0, localNumArrows):
            gv.ax1.draw_artist(mg.arrowOutput[i])
    for p in gv.ax1.patches:
        # Draw input/output circles last, so they will be on top
        gv.ax1.draw_artist(p)

    # This draws the restored canvas background, and draw_artist() items, but not patch collections.
    #    gv.ax1.draw_artist(gv.ax1.patch)
    canvas.blit(gv.ax1.bbox)

    # Render to screen
    if mg.get_animation_speed() == 0:
        plt.pause(0.001)  # This is necessary to redraw line segments

    # This gives better performance after initial draw
    canvas.flush_events()

    while not mg.quit_flag:

        mg.CheckArrows()
        if gv.localShowArrows != mg.flagShowArrows:
            # Just switched on or off. Hide or show arrows.
            localShowArrows = mg.flagShowArrows
            mg.setArrowVisible(gv.localShowArrows)

        if localAnimateSpeed != mg.get_animation_speed():
            # Animation has just changed.
            if localAnimateSpeed == 0:
                # If currently paused, update t1 or else we'll get a large jump
                t1 = 0
            localAnimateSpeed = mg.get_animation_speed()
            DesiredStepsPerSecond_actual = DesiredStepsPerSecond_base * localAnimateSpeed

        if mg.flagChangeMatrix:
            mg.flagChangeMatrix = False

            if mg.flagRestartAnimation:
                LINE_X = mg.flagX
                inputVal = LINE_X+0j
                currentCycle = 0
                currentPoint = 0
                mg.flagRestartAnimation = False

            elif mg.flagX is not None and mg.flagY is not None:  # Because x, y won't be valid if mouse out of bounds

                if mg.whichRowToAdjust == 0:
                    # Left mouse moves single input point
                    inputVal = complex(LINE_X, mg.flagY * INPUT_SCALE)
                    break
                else:
                    # Right mouse moves entire line
                    LINE_X = mg.flagX * INPUT_SCALE
                    inputVal = complex(LINE_X, mg.flagY * INPUT_SCALE)

                    mg.flagRecalculate = True

        if mg.flagRecalculate:
            mg.flagRecalculate = False
            INPUT_SCALE = mg.INPUT_RANGES[mg.flagInputRangeIndex]

            mg.Restore_background(canvas)
            th1b, line1 = MakeLine()
            gv.pltLineCollection, th2b, line2 = gObjects.InitAnimation()

#            mg.multiline(np.real(th1b), np.imag(th1b), line1, color=gv.linecolors)
#            mg.multiline(np.real(th2b), np.imag(th2b), line2, color=gv.linecolors)
            break

        if mg.flagRedrawAxes:
            mg.flagRedrawAxes = False
            mg.set_zoom()
            canvas.draw()
            break

        if localAnimateSpeed > 0:
            break

        canvas.flush_events()

    if mg.quit_flag:
        break

    if localAnimateSpeed > 0:

        if t1 == 0:
            # Just resumed animation, don't increment
            t1 = time.time()
            continue

        t2 = time.time()
        tdiff = t2 - t1
        t1 = t2

        currentPoint = currentPoint + DesiredStepsPerSecond_actual * tdiff
        if currentPoint >= gv.line_length:
            currentPoint = 0.0
            currentCycle = currentCycle + 1

plt.close()
