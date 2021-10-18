import tkinter as tk
from tkinter import ttk
from functools import partial, partialmethod
from os import environ

import riemann_heatmap as rh


class WinHeatMap:
    def __init__(self, win, is_android):

        # global rh.fig_mgr

        self.win = win
        self.is_android = is_android

        win.title(string="Complex functions")

        #
        # Controllers at top of window
        #
        # These include sliders, buttons, and checkboxes
        #
        frame_top_controls = tk.Frame(win, highlightbackground="black",
                                      highlightthickness=1, relief="flat", borderwidth=5)
        frame_top_controls.pack(side=tk.TOP, fill=tk.X, padx=8, pady=5)

        self.label1 = ttk.Label(frame_top_controls, text="Resolution =")
        self.label1.pack()

        self.slider1 = ttk.Scale(frame_top_controls, from_=0.05, to=1, command=self.do_density_slider)
        self.slider1.bind("<ButtonRelease-1>", self.do_density_slider_button_release)
        self.slider1.pack(fill=tk.X)

        self.labelA = ttk.Label(frame_top_controls, text="Parameter A")
        self.labelA.pack()

        # slider for parameter A
        self.sliderA = ttk.Scale(frame_top_controls, from_=-10, to=10, command=self.do_parameterA)
        self.sliderA.bind("<ButtonRelease-1>", self.do_sliderA_button_release)
        self.sliderA.pack(fill=tk.X)

        self.labelB = ttk.Label(frame_top_controls, text="Parameter B")
        self.labelB.pack()

        # slider for parameter B
        self.sliderB = ttk.Scale(frame_top_controls, from_=-10, to=10, command=self.do_parameterB)
        self.sliderB.bind("<ButtonRelease-1>", self.do_sliderB_button_release)
        self.sliderB.pack(fill=tk.X)

        self.label2 = ttk.Label(frame_top_controls, text="Riemann iterations =")
        self.label2.pack()

        #
        # Frame for iterations
        #
        # This contains both a scale/slider, and a spinbox
        #
        frame_iter = ttk.Frame(frame_top_controls)
        frame_iter.pack(fill=tk.X)

        # slider for iterations
        self.slider_iter = ttk.Scale(frame_iter, from_=20, to=200, command=self.update_iter_label)
        self.slider_iter.bind("<ButtonRelease-1>", self.do_slider_iter_button_release)
        self.slider_iter.pack(side=tk.LEFT, expand=1, fill=tk.X, pady=5)

        self.spin_iter = ttk.Spinbox(frame_iter, from_=20, to=2000, command=self.do_iter_spin, width=7)
        self.spin_iter.bind('<Return>', self.do_iter_spin_event)
        self.spin_iter.bind('<FocusOut>', self.do_iter_spin_event)
        self.spin_iter.pack(side=tk.LEFT)  # , fill=tk.X)

        self.progress = ttk.Progressbar(frame_top_controls, orient="horizontal", length=100, mode="determinate")
        self.progress.pack(fill=tk.X)

        # Frame for plot range spinners and accompanying label
        frame_spinner = ttk.Frame(frame_top_controls)
        frame_spinner.pack(fill=tk.X, pady=5)

        # spinboxes for plot range
        #
        # I initially used a textvariable to control spinner, and put this variable into riemann_heatmap.settings.
        # I thought this would simplify code since there would then be no need to create a separate state variable
        # that I would have to manually update. However, this has two disadvantages:
        #   1. Lots of other code (well, just make_plot) has to know that this is a tkinter.StringVar, and not a
        #      normal number.
        #   2. It is hard upon entering callback whether value was just changed, versus whether user has simply
        #      changed focus without changing anything. This is important to avoid extraneous auto-recalculations.
        #      To track this, we need to create a separate variable to hold "old" state, and manually update it.
        #      This obviates all the work we saved by using the textvariable to begin with.
        #

        # Plot x range
        self.label_spinner2 = ttk.Label(frame_spinner, text="Plot x ± range: ")
        self.label_spinner2.pack(side=tk.LEFT)

        # spinbox for plot x range
        # rh.settings.plot_range_x = tk.StringVar(win)
        self.spin_x_range = ttk.Spinbox(frame_spinner, from_=1, to=100, command=self.do_spin_x_range, width=5)
        self.spin_x_range.pack(side=tk.LEFT)

        # Need to bind keys or else value doesn't update
        self.spin_x_range.bind('<Return>', self.do_spin_x_range_event)
        self.spin_x_range.bind('<FocusOut>', self.do_spin_x_range_event)

        self.label_spinner = ttk.Label(frame_spinner, text="  Plot y ± range: ")
        self.label_spinner.pack(side=tk.LEFT)
        self.spin_y_range = ttk.Spinbox(frame_spinner, from_=1, to=100, command=self.do_spin_y_range, width=5)
        self.spin_y_range.pack(side=tk.LEFT)

        # Need to bind keys or else value doesn't update
        self.spin_y_range.bind('<Return>', self.do_spin_y_range_event)
        self.spin_y_range.bind('<FocusOut>', self.do_spin_y_range_event)
        #            self.spin1.bind('<FocusIn>', do_spin1_event)

        self.label_spinner = ttk.Label(frame_spinner, text="  Plot y center: ")
        self.label_spinner.pack(side=tk.LEFT)
        self.spin_y_center = ttk.Spinbox(frame_spinner, from_=1, to=100, command=self.do_spin_y_center, width=5)
        self.spin_y_center.pack(side=tk.LEFT)

        # Need to bind keys or else value doesn't update
        self.spin_y_center.bind('<Return>', self.do_spin_y_center_event)
        self.spin_y_center.bind('<FocusOut>', self.do_spin_y_center_event)
        #            self.spin1.bind('<FocusIn>', do_spin1_event)

        #
        # Row of check boxes controlling critical strip plotting, positive y plotting,
        # and auto-recalculation
        #
        frame_checks_plot = tk.Frame(frame_top_controls)
        frame_checks_plot.pack(fill=tk.X, pady=5)

        # Narrow x only
        self.var_critical_strip = tk.IntVar(win)
        self.checkbox_critical_strip = ttk.Checkbutton(frame_checks_plot,
                                                       text="0 < x < 1 only (critical strip)",
                                                       command=self.do_critical_strip,
                                                       variable=self.var_critical_strip)
        self.checkbox_critical_strip.pack(side=tk.LEFT)

        # Aspect ratio checkbox
        rh.settings.keep_square = tk.IntVar(win)
        self.checkbox_keep_square = ttk.Checkbutton(frame_checks_plot,
                                                    text="1:1 aspect ratio",
                                                    command=self.do_square,
                                                    variable=rh.settings.keep_square)
        self.checkbox_keep_square.pack(side=tk.LEFT)

        #
        # Row of buttons
        #
        self.frame_buttons = tk.Frame(frame_top_controls)
        self.frame_buttons.pack()  # fill=tk.X)
        self.button_recalculate = ttk.Button(self.frame_buttons, text="Recalculate", command=self.recalculate)
        self.button_recalculate.pack(side=tk.LEFT)
        self.button_cancel = ttk.Button(self.frame_buttons, text="Cancel", command=self.do_cancel)
        self.button_cancel.pack(side=tk.LEFT)

        #
        # Row of check boxes
        #
        frame_checks = tk.Frame(frame_top_controls)
        frame_checks.pack()  # fill=tk.X)

        self.var_phase = tk.IntVar(win)
        self.checkbox_phase = ttk.Checkbutton(frame_checks,
                                              text="Phase only",
                                              command=self.do_checkbox_phase_only,
                                              variable=self.var_phase)
        self.checkbox_phase.pack(side=tk.LEFT)

        self.var_magnitude_only = tk.IntVar(win)
        self.checkbox_magnitude_only = ttk.Checkbutton(frame_checks,
                                                       text="Mag only",
                                                       command=self.do_checkbox_magnitude_only,
                                                       variable=self.var_magnitude_only)
        self.checkbox_magnitude_only.pack(side=tk.LEFT)

        self.var_oversample = tk.IntVar(win)
        self.checkbox_oversample = ttk.Checkbutton(frame_checks,
                                                   text="4x oversample",
                                                   command=self.do_checkbox_oversample,
                                                   variable=self.var_oversample)
        self.checkbox_oversample.pack(side=tk.LEFT)

        self.var_auto_recalculate = tk.IntVar(win)
        self.checkbox_auto_recalculate = ttk.Checkbutton(frame_checks,
                                                         text="Auto recalc",
                                                         command=self.do_checkbox_auto_recalculate,
                                                         variable=self.var_auto_recalculate)
        self.checkbox_auto_recalculate.pack(side=tk.LEFT)

        #
        #  Grid of Buttons and controls for various graph types
        #
        frame2 = tk.Frame(win, highlightbackground="black", highlightthickness=1, relief="flat", borderwidth=5)
        frame2.pack(side=tk.TOP, padx=8, pady=5)

        row_num = 0
        col_num = 0
        self.button_ax_list = {}
        for k in rh.plot_list:
            self.button_ax_list[k] = ttk.Button(frame2, text=rh.plot_list[k],
                                                # For some reason, lambda doesn't work, but partial does
                                                command=partial(self.do_button, k))
            self.button_ax_list[k].grid(column=col_num, row=row_num, sticky=tk.E + tk.W)

            col_num += 1
            if col_num == 2:
                col_num = 0
                row_num += 1

        #                tmp = ttk.Label(frame2, text="controls for " + str(k))
        #                tmp.grid(column=1, row=row_num)
        #                row_num = row_num + 1

        # Bottom controls
        frame_bottom_controls = tk.Frame(win)
        frame_bottom_controls.pack(side=tk.TOP)

        self.button_q = ttk.Button(frame_bottom_controls, text="Close", command=self.do_quit)
        self.button_q.pack()

        rh.update_progress_callback = self.update_slider
        rh.rm.computation_progress_callback = rh.update_computation_status

    def set_initial_values(self):

        rh.make_figure_manager(self.win, self.is_android)

        # These must go outside constructor, as they will trigger callback
        # which needs access to mainApp object
        self.slider1.set(rh.settings.MESH_DENSITY)

        # Set iterations in spinbox and slider
        self.slider_iter.set(rh.rm.RIEMANN_ITER_LIMIT)
        self.spin_iter.set(rh.rm.RIEMANN_ITER_LIMIT)

        # Checkboxes for general plot behavior
        self.var_phase.set(rh.settings.phase_only)
        self.var_magnitude_only.set(rh.settings.magnitude_only)
        self.var_oversample.set(rh.settings.oversample)
        self.var_auto_recalculate.set(rh.settings.auto_recalculate)

        # Checkboxes for plot domain range
        self.var_critical_strip.set(rh.settings.critical_strip)
        rh.settings.keep_square.set(1)

        self.sliderA.set(rh.settings.parameterA)
        self.sliderB.set(rh.settings.parameterB)

        self.spin_y_range.set(rh.settings.plot_range_y)
        self.spin_x_range.set(rh.settings.plot_range_x)
        self.spin_y_center.set(rh.settings.plot_y_center)

    def update_slider(self, percent):
        #        print("  \r" + '{:1.2f}'.format(percent) + "%", end="")
        self.progress['value'] = percent
        self.win.update()

    # partialmethod() does not work in combination with tkinter. So this is not used.
    # Instead, there is a function defined outside this class that handles button press.
    def do_button(self, wid):
        print()
        print("User chose plot: " + rh.plot_list[wid])
        rh.settings.last_selection = wid
        rh.make_plot(_selection=wid)

    def do_square(self):
        # Because we have a backing variable tracking this, we don't need to update any other variable
        val = rh.settings.keep_square.get()
        print("Keep square val: " + str(val))

        if val == True:
            # Ensure that y-range is same as x-range
            rh.settings.plot_range_y = rh.settings.plot_range_x
            self.spin_y_range.set(rh.settings.plot_range_y)
            # Turn off critical strip option, which is incompatible with this one.
            rh.settings.critical_strip = not val
            self.var_critical_strip.set(not val)

    def do_quit(self):

        # Need this or else control will not return to mainloop
        rh.rm.quit_computation_flag = True

        # This closes heatmap window once computations return
        self.win.destroy()

    def do_density_slider_button_release(self, _event):
        val = self.slider1.get()
        self.do_density_slider(val)
        if rh.settings.auto_recalculate:
            self.recalculate()

    def do_density_slider(self, val):
        # Note: val is a string when called during mouse move, and a float when
        # called during button release
        rh.settings.MESH_DENSITY = float(val)
        if __name__ == "__main__":
            # In test mode, rh.fig_mgr.screen_y_pixels will not exist, so substitute 1000
            val = float(val) * 1000
        else:
            val = float(val) * rh.fig_mgr.screen_y_pixels
        st = str(int(float(val)))
        self.label1['text'] = "Resolution = " + st + " x " + st

    def do_sliderA_button_release(self, _event):
        self.do_parameterA(self.sliderA.get())
        if rh.settings.auto_recalculate:
            self.recalculate()

    def do_parameterA(self, val):
        rh.settings.parameterA = float(val)
        self.labelA['text'] = "Parameter A = " + str(val)

    def do_sliderB_button_release(self, _event):
        self.do_parameterB(self.sliderB.get())
        if rh.settings.auto_recalculate:
            self.recalculate()

    def do_parameterB(self, val):
        rh.settings.parameterB = float(val)
        self.labelB['text'] = "Parameter B = " + str(val)

    def do_slider_iter_button_release(self, _event):
        # User has let go of button. Now we can really update things

        # Convert string to int, then back to string
        val_str = self.slider_iter.get()
        val_int = int(float(val_str))
        val_str = str(val_int)

        # Set label text and spinbox
        self.update_iter_label(val_str)
        self.spin_iter.set(val_str)
        rh.rm.RIEMANN_ITER_LIMIT = val_int
        rh.rm.precompute_coeffs()
        if rh.settings.auto_recalculate:
            self.recalculate()

    def update_iter_label(self, val):
        # Callback associated with slider/scale object
        #
        # We come here whenever slider is updated, either by GUI input or by code
        # This updates text label, but does NOT update anything else, for two reasons:
        # 1. When user moves slider, it generates lots of calls for every intermediate value.
        # 2. When value is set by outside code (e.g. slider_iter.set()), we don't want to update value,
        #    e.g. if spinner has selected value > 200
        val_int = int(float(val))
        val_str = str(val_int)

        if val_int == self.slider_iter['to']:
            # If we have reached max value, then show ">="
            self.label2['text'] = "Riemann iter >= " + val_str
        else:
            self.label2['text'] = "Riemann iter = " + val_str

    def do_iter_spin(self):
        # Callback to handle user clicking spin arrows
        #
        # calling spin_iter.set() does NOT seem to come here, whereas calling slider_iter.set() does seem to
        # trigger the callback
        #
        val_int = int(self.spin_iter.get())
        self.slider_iter.set(str(val_int))
        rh.rm.RIEMANN_ITER_LIMIT = val_int
        # Don't recalculate here, as user might click spinner arrows multiple times in rapid succession

    def do_iter_spin_event(self, _event):
        # Handle user text entry
        val_int = int(self.spin_iter.get())
        if rh.rm.RIEMANN_ITER_LIMIT != val_int:
            # Handle only if value has changed
            self.slider_iter.set(str(val_int))
            rh.rm.RIEMANN_ITER_LIMIT = val_int
            if rh.settings.auto_recalculate:
                self.recalculate()

    def recalculate(self):
        # Handle user pressing Recalculate button, or programmatic recalculations if auto-calculate is checked
        wid = rh.settings.last_selection
        if wid < 0:
            # No previous selection exists, e.g. we have just launched program
            return
        rh.settings.REUSE_FIGURE = True
        rh.make_plot(wid)

    def do_remap_color(self):
        # Handle user pressing Recalculate button, or programmatic recalculations if auto-calculate is checked
        if rh.settings.last_selection < 0:
            # No previous selection exists, e.g. we have just launched program
            return
        rh.map_color()

    def do_cancel(self):
        rh.rm.quit_computation_flag = True

    def do_spin_y_range(self):
        # User has clicked on spinner arrows
        # This is the y-range spinner
        val = self.spin_y_range.get()
        float_val = float(val)
        rh.settings.plot_range_y = float_val
        if rh.settings.keep_square.get():
            self.spin_x_range.set(val)
            rh.settings.plot_range_x = float_val

    # Attempting to make a generic updater ... this does not work because range1, range2 are floats, which are immutable
    # and so changes in here won't propagate out
    def update_range(self, spin1, spin2, range1, range2):
        val = spin1.get()
        float_val = float(val)
        if range1 != float_val:
            range1 = float_val
            # Only update if it has changed.
            if rh.settings.keep_square.get():
                range2 = float_val
                spin2.set(val)
            if rh.settings.auto_recalculate:
                self.recalculate()

    def do_spin_y_range_event(self, _event):
        # User has moved focus out of spinner, or hit enter key
        # This is the y-range spinner
        val = self.spin_y_range.get()  # Get value from GUI. Value is text.
        float_val = float(val)
        if rh.settings.plot_range_y != float_val:
            # If changed from before, then store new value in rh.settings,
            # then handle possible updating of x-range (to keep square)
            # and possible auto-recalculation.
            rh.settings.plot_range_y = float_val
            if rh.settings.keep_square.get():
                rh.settings.plot_range_x = float_val
                self.spin_x_range.set(val)
            if rh.settings.auto_recalculate:
                self.recalculate()

    def do_spin_x_range(self):
        # User has clicked on spinner arrows
        # This is the x-range spinner
        val = self.spin_x_range.get()
        float_val = float(val)
        rh.settings.plot_range_x = float_val
        if rh.settings.keep_square.get():
            self.spin_y_range.set(val)
            rh.settings.plot_range_y = float_val

    def do_spin_x_range_event(self, _event):
        # User has moved focus out of spinner, or hit enter key
        # This is the x-range spinner
        val = self.spin_x_range.get()
        float_val = float(val)
        if rh.settings.plot_range_x != float_val:
            rh.settings.plot_range_x = float_val
            # Only update if it has changed.
            if rh.settings.keep_square.get():
                rh.settings.plot_range_y = float_val
                self.spin_y_range.set(val)
            if rh.settings.auto_recalculate:
                self.recalculate()

    def do_spin_y_center(self):
        # User has clicked on spinner arrows
        # This is the y-center spinner
        val = self.spin_y_center.get()
        float_val = float(val)
        rh.settings.plot_y_center = float_val

    def do_spin_y_center_event(self, _event):
        # User has moved focus out of spinner, or hit enter key
        # This is the y-center spinner
        val = self.spin_y_center.get()
        float_val = float(val)
        if rh.settings.plot_y_center != float_val:
            rh.settings.plot_y_center = float_val
            # Only update if it has changed.
            if rh.settings.auto_recalculate:
                self.recalculate()

    def do_checkbox_phase_only(self):
        # Handle user checking-unchecking this box
        rh.settings.phase_only = self.var_phase.get()  # not rh.Settings.phase_only
        if rh.settings.phase_only:
            # If checking, then uncheck the other one ... should we use a radio button instead?
            rh.settings.magnitude_only = 0
            self.var_magnitude_only.set(0)
        self.do_remap_color()

    def do_checkbox_magnitude_only(self):
        # Handle user checking-unchecking this box
        rh.settings.magnitude_only = self.var_magnitude_only.get()  # not rh.Settings.phase_only
        if rh.settings.magnitude_only:
            # If checking, then uncheck the other
            rh.settings.phase_only = 0
            self.var_phase.set(0)
        self.do_remap_color()

    def do_checkbox_oversample(self):
        # Handle user checking-unchecking this box
        rh.settings.oversample = self.var_oversample.get()  # not rh.Settings.oversample
        if rh.settings.auto_recalculate:
            self.recalculate()

    def do_checkbox_auto_recalculate(self):
        # Handle user checking-unchecking this box
        rh.settings.auto_recalculate = self.var_auto_recalculate.get()  # not rh.Settings.auto_recalculate

    def do_critical_strip(self):
        rh.settings.critical_strip = self.var_critical_strip.get()
        if rh.settings.critical_strip == True:
            # Turn off 1:1 aspect ratio, which is incompatible with this option
            rh.settings.keep_square.set(False)

            # Set y-range to 8
            self.spin_y_range.set("8")
            rh.settings.plot_range_y = 8.0

            if rh.settings.plot_y_center == 0:
                # If y-center is 0, then set to 20 to get more interesting region.
                self.spin_y_center.set("20")
                rh.settings.plot_y_center = 20.0
        if rh.settings.auto_recalculate:
            self.recalculate()

    def make_heatmap_gui(self):
        # User lower value to speed up calculations
        # Must set this value before creating TopLevel window, as that will
        # use it to set slider/scale
        rh.rm.RIEMANN_ITER_LIMIT = 100

        self.set_initial_values()


def do_main():
    root = tk.Tk()

    if 'ANDROID_BOOTLOGO' in environ:
        is_android = True

    win = WinHeatMap(root, is_android)
    win.make_heatmap_gui()
    root.mainloop()


# Make this topmost. This is needed if there is text entry beforehand, which sends focus away from Tkinter objects
# root.lift()
# root.attributes('-topmost', True)
# root.after_idle(root.attributes, '-topmost', False)
# root.mainloop()
if __name__ == "__main__":
    do_main()
else:
    print(__name__)
