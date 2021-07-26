import matplotlib as mpl
import time

import ButtonManager as bm


if True:
    print('\nSelect option\n' +
          '1. Riemann vectors\n' +
          '2. Riemann heatmap\n')

    ans = input("Select option: ")
else:
    ans = "2"

if ans == "1":
    import RiemannVectors
elif ans == "2":

    import tkinter as tk
    from tkinter import ttk
    import RiemannHeatmap as rh
    from functools import partial

    def do_button(wid):
#        print("Button " + str(wid))
        rh.Settings.last_selection = wid
        rh.make_plot(_selection=wid)

    def do_quit():
        # Need this or else control will not return to mainloop
        rh.rm.quit_computation_flag = True

        # This quits mainApp, but only after control returns to mainloop
        mainApp.quit()

    def do_density_slider(val):
        rh.Settings.MESH_DENSITY=float(val)
        res = float(val)*rh.bmgr.screen_y_pixels
        st = str(int(res))
        mainApp.label1['text'] = "Resolution = " + st + " x " + st

    def do_iter_slider(val):
        rh.rm.RIEMANN_ITER_LIMIT=int(float(val))
        mainApp.label2['text'] = "Riemann iter = " + str(rh.rm.RIEMANN_ITER_LIMIT)
        rh.rm.precompute_coeffs()

    def do_recalc():
        wid = rh.Settings.last_selection
        rh.Settings.REUSE_FIGURE = True
        rh.make_plot(wid)

    def do_cancel():
        rh.rm.quit_computation_flag = True

    def update_slider(percent, update_only=False):

        if update_only:
            mainApp.update()
            return

#        print("  \r" + '{:1.2f}'.format(percent) + "%", end="")
        mainApp.progress['value']=percent
        mainApp.update()
#        mainApp.update_idletasks()

    def do_phase():
        rh.Settings.phase_only = not rh.Settings.phase_only

    def do_oversample():
        rh.Settings.oversample = not rh.Settings.oversample

    rh.bmgr = bm.ButtonManager(use_matplotlib_widgets=False)

    class MainApp(tk.Tk):
        def __init__(self, *args, **kwargs):
            tk.Tk.__init__(self, *args, **kwargs)

            self.title(string="Complex functions")

            #
            # Controllers at top of window
            #
            # These includes sliders, buttons, and checkboxes
            #
            frame_top_controls = tk.Frame(self, highlightbackground="black", highlightthickness=1, relief="flat", borderwidth=5)
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
            self.button_recalc = ttk.Button(self.frame_buttons, text="Recalculate",command=do_recalc)
            self.button_recalc.pack(side=tk.LEFT)
            self.button_cancel = ttk.Button(self.frame_buttons, text="Cancel",command=do_cancel)
            self.button_cancel.pack(side=tk.LEFT)

            #
            # Row of check boxes
            #
            self.frame_checks = tk.Frame(frame_top_controls)
            self.frame_checks.pack(fill=tk.X)

            self.var_phase = tk.IntVar(self)
            self.checkbox_phase = ttk.Checkbutton(self.frame_checks,
                                                  text="Phase only",
                                                  command=do_phase,
                                                  variable=self.var_phase)
            self.checkbox_phase.pack(side=tk.LEFT)

            self.var_oversample = tk.IntVar(self)
            self.checkbox_oversample = ttk.Checkbutton(self.frame_checks,
                                                       text="4x oversample",
                                                       command=do_oversample,
                                                       variable=self.var_oversample)
            self.checkbox_oversample.pack(side=tk.LEFT)

            #
            #  Grid of Buttons and controls for various graph types
            #
            frame2 = tk.Frame(self, highlightbackground="black", highlightthickness=1, relief="flat", borderwidth=5)
            frame2.pack(side=tk.TOP, padx=8, pady=5)

            row_num = 0
            self.button_ax_list = {}
            for k in rh.plot_list:
                self.button_ax_list[k] = ttk.Button(frame2, text=rh.plot_list[k],
                                                 # For some reason, lambda doesn't work, but partial does
                                                 command=partial(do_button, k))
                self.button_ax_list[k].grid(column=0, row=row_num,sticky=tk.E+tk.W)

                tmp = ttk.Label(frame2, text="controls for " + str(k))
                tmp.grid(column=1,row=row_num)

                row_num = row_num + 1

            # Botton controls
            frame_bottom_controls = tk.Frame(self)
            frame_bottom_controls.pack(side=tk.TOP)

            self.button_q = ttk.Button(frame_bottom_controls, text="Quit", command=lambda:do_quit())
            self.button_q.pack()

            rh.update_progress_callback = update_slider
            rh.rm.computation_progress_callback = rh.update_computation_status

        def set_initial_vals(self):
            # These must go outside constructor, as they will trigger callback
            # which needs access to mainApp object
            self.slider1.set(rh.Settings.MESH_DENSITY)
            self.slider2.set(rh.rm.RIEMANN_ITER_LIMIT)
            self.var_phase.set(0)
            self.var_oversample.set(0)



    mainApp = MainApp()
    mainApp.set_initial_vals()
    mainApp.lift()
    mainApp.attributes('-topmost', True)
    mainApp.after_idle(mainApp.attributes, '-topmost', False)
    mainApp.mainloop()
