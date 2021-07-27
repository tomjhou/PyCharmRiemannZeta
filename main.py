import tkinter as tk
from tkinter import ttk
import matplotlib as mpl
import mpl_figure_manager as mfm


quit_flag = False
mpl.use('TKAgg')


def do_vectors():
    import riemann_vectors


def do_heatmap():

    import riemann_heatmap as rh
    from functools import partial

    def do_button(wid):
        rh.Settings.last_selection = wid
        rh.make_plot(_selection=wid)

    def do_quit():
        global quit_flag

        # Need this or else control will not return to mainloop
        rh.rm.quit_computation_flag = True
        quit_flag = True

        # This quits mainApp, but only after control returns to mainloop
        win.quit()

    def do_density_slider(val):
        rh.Settings.MESH_DENSITY = float(val)
        res = float(val) * rh.fig_mgr.screen_y_pixels
        st = str(int(res))
        win_heatmap.label1['text'] = "Resolution = " + st + " x " + st

    def do_iter_slider(val):
        rh.rm.RIEMANN_ITER_LIMIT = int(float(val))
        win_heatmap.label2['text'] = "Riemann iter = " + str(rh.rm.RIEMANN_ITER_LIMIT)
        rh.rm.precompute_coeffs()

    def do_recalculate():
        wid = rh.Settings.last_selection
        rh.Settings.REUSE_FIGURE = True
        rh.make_plot(wid)

    def do_cancel():
        rh.rm.quit_computation_flag = True

    def update_slider(percent, update_only=False):

        if update_only:
            root.update()
            return

        #        print("  \r" + '{:1.2f}'.format(percent) + "%", end="")
        win_heatmap.progress['value'] = percent
        root.update()

    #        mainApp.update_idletasks()

    def do_phase():
        rh.Settings.phase_only = not rh.Settings.phase_only

    def do_oversample():
        rh.Settings.oversample = not rh.Settings.oversample

    class WinHeatMap:
        def __init__(self, win):
            self.win = win

            win.title(string="Complex functions")

            #
            # Controllers at top of window
            #
            # These includes sliders, buttons, and checkboxes
            #
            frame_top_controls = tk.Frame(win, highlightbackground="black",
                                          highlightthickness=1, relief="flat", borderwidth=5)
            frame_top_controls.pack(side=tk.TOP, fill=tk.X, padx=8, pady=5)

            self.label1 = ttk.Label(frame_top_controls, text="Resolution =")
            self.label1.pack()

            self.slider1 = ttk.Scale(frame_top_controls, from_=0.1, to=1, command=do_density_slider)
            self.slider1.pack(fill=tk.X)

            self.progress = ttk.Progressbar(frame_top_controls, orient="horizontal", length=100, mode="determinate")
            self.progress.pack(fill=tk.X)

            self.label2 = ttk.Label(frame_top_controls, text="Riemann iters =")
            self.label2.pack()

            self.slider2 = ttk.Scale(frame_top_controls, from_=20, to=100, command=do_iter_slider)
            self.slider2.pack(fill=tk.X)

            #
            # Row of buttons
            #
            self.frame_buttons = tk.Frame(frame_top_controls)
            self.frame_buttons.pack(fill=tk.X)
            self.button_recalculate = ttk.Button(self.frame_buttons, text="Recalculate", command=do_recalculate)
            self.button_recalculate.pack(side=tk.LEFT)
            self.button_cancel = ttk.Button(self.frame_buttons, text="Cancel", command=do_cancel)
            self.button_cancel.pack(side=tk.LEFT)

            #
            # Row of check boxes
            #
            self.frame_checks = tk.Frame(frame_top_controls)
            self.frame_checks.pack(fill=tk.X)

            self.var_phase = tk.IntVar(win)
            self.checkbox_phase = ttk.Checkbutton(self.frame_checks,
                                                  text="Phase only",
                                                  command=do_phase,
                                                  variable=self.var_phase)
            self.checkbox_phase.pack(side=tk.LEFT)

            self.var_oversample = tk.IntVar(win)
            self.checkbox_oversample = ttk.Checkbutton(self.frame_checks,
                                                       text="4x oversample",
                                                       command=do_oversample,
                                                       variable=self.var_oversample)
            self.checkbox_oversample.pack(side=tk.LEFT)

            #
            #  Grid of Buttons and controls for various graph types
            #
            frame2 = tk.Frame(win, highlightbackground="black", highlightthickness=1, relief="flat", borderwidth=5)
            frame2.pack(side=tk.TOP, padx=8, pady=5)

            row_num = 0
            self.button_ax_list = {}
            for k in rh.plot_list:
                self.button_ax_list[k] = ttk.Button(frame2, text=rh.plot_list[k],
                                                    # For some reason, lambda doesn't work, but partial does
                                                    command=partial(do_button, k))
                self.button_ax_list[k].grid(column=0, row=row_num, sticky=tk.E + tk.W)

                tmp = ttk.Label(frame2, text="controls for " + str(k))
                tmp.grid(column=1, row=row_num)

                row_num = row_num + 1

            # Bottom controls
            frame_bottom_controls = tk.Frame(win)
            frame_bottom_controls.pack(side=tk.TOP)

            self.button_q = ttk.Button(frame_bottom_controls, text="Quit", command=lambda: do_quit())
            self.button_q.pack()

            rh.update_progress_callback = update_slider
            rh.rm.computation_progress_callback = rh.update_computation_status

        def set_initial_values(self):

            # These must go outside constructor, as they will trigger callback
            # which needs access to mainApp object
            self.slider1.set(rh.Settings.MESH_DENSITY)
            self.slider2.set(rh.rm.RIEMANN_ITER_LIMIT)
            self.var_phase.set(0)
            self.var_oversample.set(0)

    # Note that the following will close a temporary figure, causing tk.mainloop to quit.
    rh.fig_mgr = mfm.MplFigureManager()

    win = tk.Toplevel(root)
    win_heatmap = WinHeatMap(win)
    win_heatmap.set_initial_values()

    # User lower value to speed up calculations
    rh.rm.RIEMANN_ITER_LIMIT = 30

    # Make this topmost. This is needed if there is text entry beforehand, which sends focus away from Tkinter objects
    # root.lift()
    # root.attributes('-topmost', True)
    # root.after_idle(root.attributes, '-topmost', False)
    # root.mainloop()


def do_exit():
    global quit_flag
    quit_flag = True
    root.quit()


root = tk.Tk()
root.title(string="Choose")
frame1 = tk.Frame(root, highlightbackground="black", highlightthickness=1, relief="flat", borderwidth=5)
frame1.pack(side=tk.TOP, fill=tk.BOTH, padx=20, pady=20)

ttk.Label(frame1, text="Choose program").pack(fill=tk.X, pady=5)
ttk.Button(frame1, text="Riemann vectors", command=do_vectors).pack(fill=tk.X, padx=10, pady=5)
ttk.Button(frame1, text="Heatmaps", command=do_heatmap).pack(fill=tk.X, padx=10, pady=5)
ttk.Button(frame1, text="Exit", command=do_exit).pack(fill=tk.X, padx=10, pady=5)

# Because matplotlib.close() will terminate this loop, we wrap it in another loop
while not quit_flag:
    tk.mainloop()

print("Finished")
