#
# Calculates Riemann function for values along vertical line in complex plane. Animates result in complex number plane.
#
# Calculation uses Euler's transformation to evaluate Dirichlet eta function, then multiples by 1/(1-2^(1-s))
# to get Riemann zeta. This converges much faster. Calculation is further sped up by precomputing the weighted
# binomial coefficient sum used in Euler transform.
#

import matplotlib as mpl
import matplotlib.collections
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.pylab as pl
import numpy as np
import time
import typing

# My files
import riemann_vectors_gui as mg
import riemann_math as rm
import matplotlib.lines as mlines


class RiemannVectors(mg.GraphicsObjects):
    # Change the default cursor to any valid TK cursor
    # To hide it, you'd use the string "none" (or possibly "no" on windows)
    old_cursor = None
    # If true, generates a nice smooth demo for background screensaver. Very computationally expensive
    DEMO_MODE = False
    rm.RIEMANN_ITER_LIMIT = 150  # Use higher value to get more precision

    def __init__(self):

        if self.DEMO_MODE:
            # Generate demo that we can use to make background screensaver
            rm.RIEMANN_ITER_LIMIT = 600
            self.flagAxisScaleIndex = len(mg.AXIS_SCALES) - 1  # Zoomed out
            self.flagInputRangeIndex = 2  # Max input range
        else:
            # Standard mode, less computationally expensive
            rm.RIEMANN_ITER_LIMIT = 150
            self.flagAxisScaleIndex = len(mg.AXIS_SCALES) - 1  # Zoomed out
            self.flagInputRangeIndex = 1  # Normal (reduced) input range

        rm.precompute_coeffs()

        # Create figure graphics
        self.PLOT_SIZE = np.pi / 2

        super().__init__(self.PLOT_SIZE, not self.DEMO_MODE)

        self.INPUT_LINE_RESOLUTION = 40  # Number of points per unit on input line

        self.OUTPUT_ARROW_SCALE = 0.03
        self.INPUT_ARROW_SCALE = 0.03
        self.CIRCLE_PATCH_SCALE: float = 0.007
        self.INPUT_SCALE = mg.INPUT_RANGES[self.flagInputRangeIndex]

        self.InputStyle = mg.PatchStyle.ARROW_CIRCLE

        print('\nRiemann-Zeta demo instructions:\n\n' +
              'Left mouse button drags red input vector\n' +
              'Right mouse button drags vertical line\n' +
              'Space bar to pause/unpause animation\n' +
              'Up/down arrows to change animation velocity\n' +
              '"r" to reset animation\n' +
              '"z" to zoom out\n' +
              '"x" to exit')

        self.LINE_X = self.flagX

    def wait_cursor(self, v=True):
        """
        Change cursor to hourglass so user knows that computations are happening.
        Parameters
        ----------
        v

        Returns
        -------

        """
        if v:
            try:
                # PyCharm complains that _last_cursor is protected element
                self.old_cursor = self.canvas.toolbar._last_cursor
                self.canvas.toolbar.set_cursor(1)
            except AttributeError as e:
                print("Error setting hourglass cursor: " + str(e))
        else:
            if self.old_cursor is not None:
                self.canvas.toolbar.set_cursor(self.old_cursor)

    def MakeLine2(self) -> [list, matplotlib.collections.LineCollection]:
        self.th1, self.line_colors = super().MakeLine1(self.PLOT_SIZE, self.INPUT_SCALE,
                                                      self.INPUT_LINE_RESOLUTION, self.LINE_X)
        th1_scaled = [x2 / self.INPUT_SCALE for x2 in self.th1]
        # Create input rainbow line
        line = mg.multiline(np.real(th1_scaled), np.imag(th1_scaled), color=self.line_colors, linewidths=0.5)

        return th1_scaled, line

    def make_graph(self):

        th1b, line1 = self.MakeLine2()
        plt.pause(.001)  # Draws items to screen

        tmp, numArrowsPlusOne = rm.riemann(2, True)  # Second argument returns required array size

        localNumArrows = 0
        inputArrowPatch = None

        # Create red "input" arrow
        if self.InputStyle == mg.PatchStyle.ARROW or self.InputStyle == mg.PatchStyle.ARROW_CIRCLE:
            inputArrowPatch = None  # mpatches.Arrow(0, 0, 0, 0, color=mg.INPUT_VECTOR_COLOR, width=INPUT_PATCH_SCALE)

        inputCirclePatch: typing.Optional[mpl.patches.Circle] = None
        inputCirclePatch2: typing.Optional[mpl.patches.Circle] = None
        inputCirclePatch3: typing.Optional[mpl.patches.Circle] = None
        outputCirclePatch: typing.Optional[mpl.patches.Circle] = None
        outputCirclePatch2: typing.Optional[mpl.patches.Circle] = None
        outputCirclePatch3: typing.Optional[mpl.patches.Circle] = None

        # Create optional side patches
        if self.InputStyle == mg.PatchStyle.CIRCLE or self.InputStyle == mg.PatchStyle.ARROW_CIRCLE:
            scaleTmp = self.GetAxisScale()

            inputCirclePatch = mpatches.Circle((0, 0), color=mg.INPUT_VECTOR_COLOR,
                                               radius=self.CIRCLE_PATCH_SCALE * scaleTmp)
            inputCirclePatch.set_edgecolor('none')
            self.ax.add_patch(inputCirclePatch)

            inputCirclePatch2 = mpatches.Circle((0, 0), color=mg.INPUT_VECTOR_COLOR,
                                                radius=self.CIRCLE_PATCH_SCALE * scaleTmp)
            inputCirclePatch2.set_edgecolor('none')
            self.ax.add_patch(inputCirclePatch2)
            inputCirclePatch3 = mpatches.Circle((0, 0), color=mg.INPUT_VECTOR_COLOR,
                                                radius=self.CIRCLE_PATCH_SCALE * scaleTmp)
            inputCirclePatch3.set_edgecolor('none')
            self.ax.add_patch(inputCirclePatch3)

            outputCirclePatch = mpatches.Circle((0, 0), color=(0, 0, 0), radius=self.CIRCLE_PATCH_SCALE * scaleTmp)
            outputCirclePatch.set_edgecolor('none')
            self.ax.add_patch(outputCirclePatch)

            outputCirclePatch2 = mpatches.Circle((0, 0), color=(0, 0, 0), radius=self.CIRCLE_PATCH_SCALE * scaleTmp)
            outputCirclePatch3 = mpatches.Circle((0, 0), color=(0, 0, 0), radius=self.CIRCLE_PATCH_SCALE * scaleTmp)
            self.ax.add_patch(outputCirclePatch2)
            self.ax.add_patch(outputCirclePatch3)

            # Side patches are invisible unless selected later
            inputCirclePatch2.set_visible(False)
            inputCirclePatch3.set_visible(False)
            outputCirclePatch2.set_visible(False)
            outputCirclePatch3.set_visible(False)

        # Create new outArray that will be used to generate arrows. Has one more element than arrow
        rm.outArray = np.zeros(numArrowsPlusOne, dtype=complex)

        # Empty output arrow collection
        outputArrowCollection = []
        self.pltLineCollection = None
        lineList = []

        # Create static objects that do not change
        self.pltLineCollection, th2b, line2 = self.InitAnimation()
        currentPoint = 0.0
        currentCycle = 0

        # Need to call this the first time or else objects won't draw later
        plt.pause(0.01)

        localAnimateSpeed = self.get_animation_speed()  # Need local copy of flag to detect when it changes state
        self.localShowArrows = self.flagShowArrows
        localShowSideDots = False

        DesiredStepsPerSecond_base = 1

        inputVal = self.th1[0]
        outputVal = rm.riemann(inputVal)
        mg.arrowOutput = []

        self.wait_cursor(True)
        # Create dynamic objects that will change every step
        for x in range(0, len(self.th1)):
            lineTmp = self.ax.add_line(mpl.lines.Line2D([0, 0], [0, 0], color=self.line_colors[x]))
            lineList.append(lineTmp)

        arrowColors = pl.cm.jet(np.linspace(0, 4, numArrowsPlusOne - 1))

        # Create arrow segments
        for i in range(1, numArrowsPlusOne):
            patchTmp = self.ax.add_line(mpl.lines.Line2D([0, 0], [0, 0],
                                                                  color=arrowColors[i - 1],
                                                                  linewidth=5 / np.log2(i + 2)))
            self.arrowOutput.append(patchTmp)

        self.wait_cursor(False)

        plt.pause(0.001)
        t1 = time.time()
        max_speed = 3
        while self.quit_flag == 0:

            if self.DEMO_MODE:
                DesiredStepsPerSecond_actual = DesiredStepsPerSecond_base / np.abs(rm.eta_zeta_scale(inputVal))
            else:
                DesiredStepsPerSecond_actual = \
                    DesiredStepsPerSecond_base * self.get_animation_speed() / np.abs(rm.eta_zeta_scale(inputVal))

            inputVal = self.GetPoint(currentPoint)
            outputVal = rm.riemann(inputVal)

            self.CheckArrows()

            # Draw head-to-tail purple arrows
            if self.localShowArrows:
                if self.GetAxisScale() > 1:
                    # Reduce # of arrows when zoomed out, as they will be invisible
                    localNumArrows = int(numArrowsPlusOne - 1)
                elif self.GetAxisScale() == 1:
                    localNumArrows = int(numArrowsPlusOne - 1)
                else:
                    localNumArrows = numArrowsPlusOne - 1
                for i in range(0, localNumArrows):
                    if self.flagCenterAtTarget:
                        s1 = outputVal - rm.outArray[i]
                        s2 = outputVal - rm.outArray[i + 1]
                        self.arrowOutput[i].set_data([s1.real, s2.real], [s1.imag, s2.imag])
                    else:
                        s1 = rm.outArray[i]
                        s2 = rm.outArray[i + 1]
                        self.arrowOutput[i].set_data([s1.real, s2.real], [s1.imag, s2.imag])

            if self.flagSideDots != localShowSideDots:
                localShowSideDots = self.flagSideDots
                inputCirclePatch2.set_visible(localShowSideDots)
                inputCirclePatch3.set_visible(localShowSideDots)
                outputCirclePatch2.set_visible(localShowSideDots)
                outputCirclePatch3.set_visible(localShowSideDots)

            # Draw input/output patches
            if self.flagCenterAtTarget:
                inputCirclePatch.set_radius(0)
                outputCirclePatch.set_radius(0)
                if inputArrowPatch is not None:
                    inputArrowPatch.remove()
                    inputArrowPatch = None
            else:
                if self.InputStyle == mg.PatchStyle.ARROW or self.InputStyle == mg.PatchStyle.ARROW_CIRCLE:
                    # Arrow patch seems to have no way to update, so just remove and replace
                    if inputArrowPatch is not None:
                        inputArrowPatch.remove()
                    inputArrowPatch = mpatches.Arrow(0, 0, inputVal.real / self.INPUT_SCALE,
                                                     inputVal.imag / self.INPUT_SCALE,
                                                     color=mg.INPUT_VECTOR_COLOR,
                                                     width=self.INPUT_ARROW_SCALE * self.GetAxisScale())
                    self.ax.add_patch(inputArrowPatch)
                if (self.InputStyle == mg.PatchStyle.CIRCLE) or (self.InputStyle == mg.PatchStyle.ARROW_CIRCLE):
                    inputCirclePatch.center = inputVal.real / self.INPUT_SCALE, inputVal.imag / self.INPUT_SCALE

                    tmpColorIndex = int(currentPoint * len(self.line_colors) / self.line_length)
                    tmpColor = self.line_colors[tmpColorIndex]

                    inputCirclePatch.set_facecolor(tmpColor)
                    inputCirclePatch.set_radius(self.CIRCLE_PATCH_SCALE * self.GetAxisScale())
                    outputCirclePatch.center = outputVal.real, outputVal.imag
                    outputCirclePatch.set_facecolor(tmpColor)
                    outputCirclePatch.set_radius(self.CIRCLE_PATCH_SCALE * self.GetAxisScale())

                    if localShowSideDots:
                        inputVal2 = inputVal - 0.2
                        inputVal3 = inputVal + 0.2
                        inputCirclePatch2.center = inputVal2.real / self.INPUT_SCALE, inputVal2.imag / self.INPUT_SCALE
                        inputCirclePatch3.center = inputVal3.real / self.INPUT_SCALE, inputVal3.imag / self.INPUT_SCALE
                        inputCirclePatch2.set_radius(self.CIRCLE_PATCH_SCALE * self.GetAxisScale() / 2)
                        inputCirclePatch3.set_radius(self.CIRCLE_PATCH_SCALE * self.GetAxisScale() / 2)

                        outputVal2 = rm.riemann(inputVal2)
                        outputVal3 = rm.riemann(inputVal3)

                        outputCirclePatch2.center = outputVal2.real, outputVal2.imag
                        outputCirclePatch3.center = outputVal3.real, outputVal3.imag
                        outputCirclePatch2.set_radius(self.CIRCLE_PATCH_SCALE * self.GetAxisScale() / 2)
                        outputCirclePatch3.set_radius(self.CIRCLE_PATCH_SCALE * self.GetAxisScale() / 2)

            self.Restore_background(self.canvas)

            # Update text
            self.textObj.update_input_complex(inputVal)
            #    for x in range(0,currentPoint):
            #        self.ax.draw_artist(lineList[x])
            if self.localShowArrows:
                for i in range(0, localNumArrows):
                    self.ax.draw_artist(self.arrowOutput[i])
            for p in self.ax.patches:
                # Draw input/output circles last, so they will be on top
                self.ax.draw_artist(p)

            # This draws the restored canvas background, and draw_artist() items, but not patch collections.
            #    self.ax.draw_artist(self.ax.patch)
            self.canvas.blit(self.ax.bbox)

            # Render to screen
            if self.get_animation_speed() == 0:
                plt.pause(0.001)  # This is necessary to redraw line segments

            # This gives better performance after initial draw
            self.canvas.flush_events()

            while not self.quit_flag:

                self.CheckArrows()
                if self.localShowArrows != self.flagShowArrows:
                    # Just switched on or off. Hide or show arrows.
                    localShowArrows = self.flagShowArrows
                    self.setArrowVisible(self.localShowArrows)

                if localAnimateSpeed != self.get_animation_speed():
                    # Animation has just changed.
                    if localAnimateSpeed == 0:
                        # If currently paused, update t1 or else we'll get a large jump
                        t1 = 0
                    localAnimateSpeed = self.get_animation_speed()
                    DesiredStepsPerSecond_actual = DesiredStepsPerSecond_base * localAnimateSpeed

                if self.flagChangeMatrix:
                    self.flagChangeMatrix = False

                    if self.flagRestartAnimation:
                        LINE_X = self.flagX
                        inputVal = LINE_X + 0j
                        currentCycle = 0
                        currentPoint = 0
                        self.flagRestartAnimation = False

                    elif self.flagX is not None and self.flagY is not None:  # x, y invalid if mouse out of bounds

                        if self.whichRowToAdjust == 0:
                            # Left mouse moves single input point
                            inputVal = complex(self.LINE_X, self.flagY * self.INPUT_SCALE)
                            break
                        else:
                            # Right mouse moves entire line
                            self.LINE_X = self.flagX * self.INPUT_SCALE
                            inputVal = complex(self.LINE_X, self.flagY * self.INPUT_SCALE)

                            self.flagRecalculate = True

                if self.flagRecalculate:
                    self.flagRecalculate = False
                    INPUT_SCALE = mg.INPUT_RANGES[self.flagInputRangeIndex]

                    self.Restore_background(self.canvas)
                    th1b, line1 = self.MakeLine2()
                    self.pltLineCollection, th2b, line2 = self.InitAnimation()
                    break

                if self.flagRedrawAxes:
                    self.flagRedrawAxes = False
                    self.set_zoom()
                    self.canvas.draw()
                    break

                if localAnimateSpeed > 0:
                    break

                self.canvas.flush_events()

            if self.quit_flag:
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
                if currentPoint >= self.line_length:
                    currentPoint = 0.0
                    currentCycle = currentCycle + 1

        plt.close()
