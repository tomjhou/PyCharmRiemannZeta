import tkinter as tk
from tkinter import ttk
from functools import partial, partialmethod

import riemann_heatmap as rh


def do_square():
    # diagnostic function only. Because we have a backing variable tracking this, we don't need to handle anything
    print("Keep square val: " + str(rh.settings.keep_square.get()))


def do_button(wid):
    print("do_button " + str(wid))
    rh.settings.last_selection = wid
    rh.make_plot(_selection=wid)


class WinHeatMap:
    def __init__(self, win):

        # global rh.fig_mgr

        self.win = win

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
        self.label_spinner2 = ttk.Label(frame_spinner, text="Plot x range: ")
        self.label_spinner2.pack(side=tk.LEFT)

        # spinbox for plot x range
        # rh.settings.plot_range_x = tk.StringVar(win)
        self.spin2 = ttk.Spinbox(frame_spinner, from_=1, to=100, command=self.do_spin2, width=5)
        self.spin2.pack(side=tk.LEFT)

        # Need to bind keys or else value doesn't update
        self.spin2.bind('<Return>', self.do_spin2_event)
        self.spin2.bind('<FocusOut>', self.do_spin2_event)

        self.label_spinner = ttk.Label(frame_spinner, text="  Plot y range: ")
        self.label_spinner.pack(side=tk.LEFT)
        self.spin1 = ttk.Spinbox(frame_spinner, from_=1, to=100, command=self.do_spin1, width=5)
        self.spin1.pack(side=tk.LEFT)

        # Need to bind keys or else value doesn't update
        self.spin1.bind('<Return>', self.do_spin1_event)
        self.spin1.bind('<FocusOut>', self.do_spin1_event)
        #            self.spin1.bind('<FocusIn>', do_spin1_event)

        self.label_spinner = ttk.Label(frame_spinner, text="  Plot y start: ")
        self.label_spinner.pack(side=tk.LEFT)
        self.spin3 = ttk.Spinbox(frame_spinner, from_=1, to=100, command=self.do_spin3, width=5)
        self.spin3.pack(side=tk.LEFT)

        # Need to bind keys or else value doesn't update
        self.spin3.bind('<Return>', self.do_spin3_event)
        self.spin3.bind('<FocusOut>', self.do_spin3_event)
        #            self.spin1.bind('<FocusIn>', do_spin1_event)

        # Aspect ratio checkbox
        rh.settings.keep_square = tk.IntVar(win)
        self.checkbox_keep_square = ttk.Checkbutton(frame_spinner,
                                                    text="1:1 aspect ratio",
                                                    command=do_square,
                                                    variable=rh.settings.keep_square)
        self.checkbox_keep_square.pack(side=tk.LEFT)

        #
        # Row of check boxes controlling critical strip plotting, positive y plotting,
        # and auto-recalculation
        #
        frame_checks_plot = tk.Frame(frame_top_controls)
        frame_checks_plot.pack(fill=tk.X, pady=5)

        # Narrow x only
        self.var_critical_strip = tk.IntVar(win)
        self.checkbox_critical_strip = ttk.Checkbutton(frame_checks_plot,
                                                       text="-2 < x < 2 only (critical strip)",
                                                       command=self.do_critical_strip,
                                                       variable=self.var_critical_strip)
        self.checkbox_critical_strip.pack(side=tk.LEFT)

        # Positive y only
        self.var_top_only = tk.IntVar(win)
        self.checkbox_top_only = ttk.Checkbutton(frame_checks_plot,
                                                 text="Positive y only",
                                                 command=self.do_positive_y_only,
                                                 variable=self.var_top_only)
        self.checkbox_top_only.pack(side=tk.LEFT)

        #
        # Row of buttons
        #
        self.frame_buttons = tk.Frame(frame_top_controls)
        self.frame_buttons.pack()  # fill=tk.X)
        self.button_recalculate = ttk.Button(self.frame_buttons, text="Recalculate", command=self.do_recalculate)
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
#                                                command=self.do_button2)
                                                command=partial(do_button, k))
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

        rh.make_figure_manager(self.win)

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
        self.var_top_only.set(rh.settings.top_only)
        self.var_critical_strip.set(rh.settings.critical_strip)
        rh.settings.keep_square.set(1)

        self.sliderA.set(rh.settings.parameterA)
        self.sliderB.set(rh.settings.parameterB)

        self.spin1.set(rh.settings.plot_range_y)
        self.spin2.set(rh.settings.plot_range_x)
        self.spin3.set(rh.settings.plot_y_start)

    def update_slider(self, percent):
        #        print("  \r" + '{:1.2f}'.format(percent) + "%", end="")
        self.progress['value'] = percent
        self.win.update()

    # partialmethod() does not work in combination with tkinter. So this is not used.
    # Instead, there is a function defined outside this class that handles button press.
    def do_button(self, wid):
        print("do_button " + str(wid))
        rh.settings.last_selection = wid
        rh.make_plot(_selection=wid)

    def do_quit(self):

        # Need this or else control will not return to mainloop
        rh.rm.quit_computation_flag = True

        # This closes heatmap window once computations return
        self.win.destroy()

    def do_density_slider_button_release(self, _event):
        val = self.slider1.get()
        self.do_density_slider(val)
        if rh.settings.auto_recalculate:
            self.do_recalculate()

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
            self.do_recalculate()

    def do_parameterA(self, val):
        rh.settings.parameterA = float(val)
        self.labelA['text'] = "Parameter A = " + str(val)

    def do_sliderB_button_release(self, _event):
        self.do_parameterB(self.sliderB.get())
        if rh.settings.auto_recalculate:
            self.do_recalculate()

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
            self.do_recalculate()

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
                self.do_recalculate()

    def do_recalculate(self):
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

    def do_spin1(self):
        # User has clicked on spinner arrows
        val = self.spin1.get()
        float_val = float(val)
        rh.settings.plot_range_y = float_val
        if rh.settings.keep_square.get():
            self.spin2.set(val)
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
                self.do_recalculate()

    def do_spin1_event(self, _event):
        # User has moved focus out of spinner, or hit enter key
        val = self.spin1.get()
        float_val = float(val)
        if rh.settings.plot_range_y != float_val:
            rh.settings.plot_range_y = float_val
            # Only update if it has changed.
            if rh.settings.keep_square.get():
                rh.settings.plot_range_x = float_val
                self.spin2.set(val)
            if rh.settings.auto_recalculate:
                self.do_recalculate()

    def do_spin2(self):
        # User has clicked on spinner arrows
        val = self.spin2.get()
        float_val = float(val)
        rh.settings.plot_range_x = float_val
        if rh.settings.keep_square.get():
            self.spin1.set(val)
            rh.settings.plot_range_y = float_val

    def do_spin2_event(self, _event):
        # User has moved focus out of spinner, or hit enter key
        val = self.spin2.get()
        float_val = float(val)
        if rh.settings.plot_range_x != float_val:
            rh.settings.plot_range_x = float_val
            # Only update if it has changed.
            if rh.settings.keep_square.get():
                rh.settings.plot_range_y = float_val
                self.spin1.set(val)
            if rh.settings.auto_recalculate:
                self.do_recalculate()

    def do_spin3(self):
        # User has clicked on spinner arrows
        val = self.spin3.get()
        float_val = float(val)
        rh.settings.plot_y_start = float_val

    def do_spin3_event(self, _event):
        # User has moved focus out of spinner, or hit enter key
        val = self.spin3.get()
        float_val = float(val)
        if rh.settings.plot_y_start != float_val:
            rh.settings.plot_y_start = float_val
            # Only update if it has changed.
            if rh.settings.auto_recalculate:
                self.do_recalculate()

    def do_checkbox_phase_only(self):
        # Handle user checking-unchecking this box
        rh.settings.phase_only = self.var_phase.get()  # not rh.Settings.phase_only
        if rh.settings.phase_only:
            # If checking, then uncheck the other one ... should we use a radio button instead?
            rh.settings.magnitude_only = 0
            self.var_magnitude_only.set(0)
        if rh.settings.auto_recalculate:
            self.do_remap_color()

    def do_checkbox_magnitude_only(self):
        # Handle user checking-unchecking this box
        rh.settings.magnitude_only = self.var_magnitude_only.get()  # not rh.Settings.phase_only
        if rh.settings.magnitude_only:
            # If checking, then uncheck the other
            rh.settings.phase_only = 0
            self.var_phase.set(0)
        if rh.settings.auto_recalculate:
            self.do_remap_color()

    def do_checkbox_oversample(self):
        # Handle user checking-unchecking this box
        rh.settings.oversample = self.var_oversample.get()  # not rh.Settings.oversample
        if rh.settings.auto_recalculate:
            self.do_recalculate()

    def do_checkbox_auto_recalculate(self):
        # Handle user checking-unchecking this box
        rh.settings.auto_recalculate = self.var_auto_recalculate.get()  # not rh.Settings.auto_recalculate
        if rh.settings.auto_recalculate:
            self.do_recalculate()

    def do_positive_y_only(self):
        rh.settings.top_only = self.var_top_only.get()
        if rh.settings.auto_recalculate:
            self.do_recalculate()

    def do_critical_strip(self):
        rh.settings.critical_strip = self.var_critical_strip.get()
        if rh.settings.auto_recalculate:
            self.do_recalculate()

    def make_heatmap_gui(self):
        # User lower value to speed up calculations
        # Must set this value before creating TopLevel window, as that will
        # use it to set slider/scale
        rh.rm.RIEMANN_ITER_LIMIT = 30

        self.set_initial_values()


def do_main():
    root = tk.Tk()
    win = WinHeatMap(root)
    win.make_heatmap_gui()
    # Place just below earlier window
    root.geometry("+5+250")
    root.mainloop()

# Make this topmost. This is needed if there is text entry beforehand, which sends focus away from Tkinter objects
# root.lift()
# root.attributes('-topmost', True)
# root.after_idle(root.attributes, '-topmost', False)
# root.mainloop()

if __name__ == "__main__":
    do_main()

