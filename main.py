import matplotlib as mpl
import time


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
        rh.settings.last_selection = wid
        rh.make_plot(_selection=wid)

    def do_quit():
        mainApp.quit()

    def do_slider(val):
        rh.settings.MESH_DENSITY=float(val)

    def do_recalc():
        wid = rh.settings.last_selection
        rh.settings.REUSE_FIGURE = True
        rh.make_plot(wid)

    def do_cancel():
        rh.rm.quit_computation_flag = True

    def update_slider(val):
        percent = val * 100 / rh.mesh_points
        print("  \r" + '{:1.2f}'.format(percent) + "%", end="")
        mainApp.progress['value']=percent
        mainApp.update()
#        mainApp.update_idletasks()

    def do_phase():
        rh.settings.phase_only = not rh.settings.phase_only

    def do_oversample():
        rh.settings.oversample = not rh.settings.oversample

    rh.bmgr = rh.bm.ButtonManager(use_matplotlib_widgets=False)

    class MainApp(tk.Tk):
        def __init__(self, *args, **kwargs):
            tk.Tk.__init__(self, *args, **kwargs)

            self.title(string="Complex functions")

            frame1 = tk.Frame(self, highlightbackground="black", highlightthickness=1, relief="flat", borderwidth=5)
            frame1.pack(side=tk.TOP, fill=tk.X, padx=8, pady=5)

            self.label1 = ttk.Label(frame1, text="Resolution =")
            self.label1.pack()

            self.slider1 = ttk.Scale(frame1, from_=0, to=1, command=do_slider)
            self.slider1.pack(fill=tk.X)
            self.slider1.set(rh.settings.MESH_DENSITY)

            self.progress = ttk.Progressbar(frame1, orient="horizontal", length=100, mode="determinate")
            self.progress.pack(fill=tk.X)

            self.frame_buttons = tk.Frame(frame1)
            self.frame_buttons.pack(fill=tk.X)
            self.button_recalc = ttk.Button(self.frame_buttons, text="Recalculate",command=do_recalc)
            self.button_recalc.pack(side=tk.LEFT)
            self.button_cancel = ttk.Button(self.frame_buttons, text="Cancel",command=do_cancel)
            self.button_cancel.pack(side=tk.LEFT)

            self.frame_checks = tk.Frame(frame1)
            self.frame_checks.pack(fill=tk.X)
            self.checkbox_phase = ttk.Checkbutton(self.frame_checks, text="Phase only", command=do_phase)
            self.checkbox_phase.pack(side=tk.LEFT)
            self.checkbox_oversample = ttk.Checkbutton(self.frame_checks, text="4x oversample", command=do_oversample)
            self.checkbox_oversample.pack(side=tk.LEFT)

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

            frame3 = tk.Frame(self)
            frame3.pack(side=tk.TOP)

            self.button_q = ttk.Button(frame3, text="Quit", command=lambda:do_quit())
            self.button_q.pack()

            rh.rm.slider_progress_callback = update_slider

    class Button2(ttk.Button):
        def __init__(self, text, _id, **kwargs):
            """
            """
            ttk.Button.__init__(self, text=text, **kwargs)

            self.id = _id

        def get_id(self):
            return self.id

    mainApp = MainApp()
    mainApp.mainloop()
