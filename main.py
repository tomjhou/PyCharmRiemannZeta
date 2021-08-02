import tkinter as tk
from tkinter import ttk
import matplotlib as mpl
from heatmap_gui import *

quit_flag = False
mpl.use('TKAgg')


def do_vectors():
    import riemann_vectors


def do_heatmap():

    win = tk.Toplevel(root)
    win_heatmap = WinHeatMap(win)
    win_heatmap.make_heatmap_gui()


def do_critical_line():
    import plot_critical_line_values


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
ttk.Button(frame1, text="Critical line plot", command=do_critical_line).pack(fill=tk.X, padx=10, pady=5)
ttk.Button(frame1, text="Exit", command=do_exit).pack(fill=tk.X, padx=10, pady=5)

root.geometry("+5+5")   # Place in very top left corner of screen

# Because matplotlib.close() will terminate this loop, we wrap it in another loop
while not quit_flag:
    tk.mainloop()
    print("Mainloop exited.", end="")
    if not quit_flag:
        print(" Reentering now.")

print("Finished")
