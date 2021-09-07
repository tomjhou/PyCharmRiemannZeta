
from functools import partial
import matplotlib as mpl
import plot_critical_line_values as pclv

mpl.use('TKAgg')

from heatmap_gui import *

quit_flag = False

# Gaps between various menu windows
PADDING_PIXELS = 5
is_android = False
if 'ANDROID_BOOTLOGO' in environ:
    is_android = True


def do_vectors():
    import riemann_vectors


def do_heatmap(root):

    win = tk.Toplevel(root)
    # Place just below earlier window

    if not is_android:
        # Move to avoid overlapping first window

        # winfo_rooty() gets y-coordinate of window contents BELOW the titlebar, whereas winfo_y(), which gets
        # coordinate of the top of the titlebar. So we use the first, in order to place new window entirely below old.
        # Noet that we also add

#        new_geom = "+%d+%d" % (PADDING_PIXELS, PADDING_PIXELS + root.winfo_rooty() + root.winfo_height())
        new_geom = "+%d+%d" % (25, 25)
        print("Creating heatmap control window at location " + new_geom)
        win.geometry(new_geom)

    win_heatmap = WinHeatMap(win, is_android)
    win_heatmap.make_heatmap_gui()


def do_critical_line():
    pclv.show_critial_line()

def do_exit(root):
    global quit_flag
    quit_flag = True
    root.quit()


def do_main():
    global quit_flag

    root = tk.Tk()
    root.title(string="Choose")

    if is_android:
        # Use simpler GUI on PyDroid, and larger buttons
        frame1 = root
        ipad = 20
    else:
        frame1 = tk.Frame(root, highlightbackground="black", highlightthickness=1, relief="flat", borderwidth=5)
        frame1.pack(side=tk.TOP, fill=tk.BOTH, padx=20, pady=20)
        ipad = 10

    ttk.Label(frame1, text="Choose program").pack(fill=tk.X, pady=5)
    ttk.Button(frame1, text="Heatmaps", command=partial(do_heatmap, root)).pack(fill=tk.X,
                                                                 padx=10, pady=5,
                                                                 ipadx=ipad, ipady=ipad)
    ttk.Button(frame1, text="Critical line plot", command=do_critical_line).pack(fill=tk.X,
                                                                                 padx=10, pady=5,
                                                                                 ipadx=ipad, ipady=ipad)
    ttk.Button(frame1, text="Riemann vectors", command=do_vectors).pack(fill=tk.X,
                                                                        padx=10, pady=5,
                                                                        ipadx=ipad, ipady=ipad)
    ttk.Button(frame1, text="Exit", command=partial(do_exit, root)).pack(fill=tk.X,
                                                          padx=10, pady=5,
                                                          ipadx=ipad, ipady=ipad)


    if is_android:
        print("Android!")
        # On PyDroid, I can't get root.geometry() to do anything. So just skip it.
    else:
        # Place in top left corner of screen
        root.geometry("+%d+%d" % (PADDING_PIXELS, PADDING_PIXELS))

    # Force window to show, otherwise winfo_geometry() will return zero
    frame1.update()
    print("Created main window with geometry " + root.winfo_geometry())

    # frame1.tkraise()

    # Clear flag in case we are running for the second time
    quit_flag = False

    # Because matplotlib.close() will terminate this loop, we wrap it in another loop
    while not quit_flag:
        tk.mainloop()
        print("Mainloop exited.")
        if not quit_flag:
            print(" Reentering now.")

    print("Finished")
    root.destroy()

#if __name__ == '__main__':   #
do_main()
