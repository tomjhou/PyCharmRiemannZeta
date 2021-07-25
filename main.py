import matplotlib as mpl
import time

print('\nSelect option\n' +
      '1. Riemann vectors\n' +
      '2. Riemann heatmap using Matplotlib buttons\n' +
      '3. Riemann heatmap using TkInter buttons\n')

ans = input("Select option: ")

if ans == "1":
    import RiemannVectors
elif ans == "2":
    import RiemannHeatmap as rh

    # Dictionaries to store widget axes.
    # Each widget needs a unique axis object used to assign callback.
    # Even though we don't use these objects after initial definition,
    # it seems we can't recycle the variable name, as the callbacks will then
    # only work for the last item assigned.
    button_ax_list = {}
    slider_ax_list = {}
    text_box_ax_list = {}

    # Dictionary to store user input from each slider
    slider_val = {}

    def AddIdButton(text, _id):
        button_ax_list[_id] = rh.bmgr.add_id_button(text, _id)
        button_ax_list[_id].on_clicked(lambda x: rh.make_fig_plot(x, button_ax_list[_id].id))


    def AddIdSlider2(text, _id, **kwargs):
        slider_ax_list[_id] = rh.bmgr.add_id_slider(text, _id, **kwargs)
        slider_ax_list[_id].on_changed(lambda x: rh.do_slider(x, slider_ax_list[_id].id))
        return slider_ax_list[_id]


    def AddIdTextBox(text, _id):
        text_box_ax_list[_id] = rh.bmgr.add_textbox(text, _id)
        text_box_ax_list[_id].on_submit(lambda x: rh.do_submit(x, text_box_ax_list[_id].id))

    rh.bmgr = rh.bm.ButtonManager(num_rows=len(rh.plot_list)+5, num_cols=2)

    cb1 = rh.bmgr.add_checkbox(rh.checkbox_list)
    cb1.on_clicked(rh.do_checkbox)

    for k in rh.plot_list:
        AddIdButton(rh.plot_list[k], k)

    b2 = rh.bmgr.add_standard_button("Quit")
    b2.on_clicked(rh.do_quit)

    #
    # Start second column of widgets
    #
    rh.bmgr.reset_button_coord(1)

    b3 = rh.bmgr.add_id_slider("mesh density",
                            _id=0, valmin=0, valmax=1,
                            valinit=rh.settings.MESH_DENSITY)
    b3.on_changed(rh.do_slider_density_update)

    AddIdButton("Recalculate", -1)
    rh.slider_progress = AddIdSlider2("% Progress", _id=-1, valmin=0, valmax=100)
    rh.rm.slider_progress_callback = rh.do_slider_progress_update
    # bmgr.increment_row()

    for k in rh.plot_list:
        AddIdSlider2("", k)

    # Button figure also responds to key press - "x" to exit
    rh.bmgr.canvas_buttons.mpl_connect('key_press_event', rh.on_keypress)

    mpl.pyplot.pause(0.01)

    while not rh.qFlag:
        rh.bmgr.canvas_buttons.flush_events()
        #        plt.pause(0.001) # This causes new plots to permanently overlay older ones
        time.sleep(0.025)  # This allows plots to be moved around more freely
