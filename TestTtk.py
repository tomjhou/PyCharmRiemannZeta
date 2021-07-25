import tkinter as tk
from tkinter import ttk


class SampleApp(tk.Tk):

    def __init__(self, *args, **kwargs):
        tk.Tk.__init__(self, *args, **kwargs)

        self.button = ttk.Button(text="start", command=self.start)
        self.button.pack()

        self.button1 = ttk.Button(text="test1", command=(lambda : self.test(1)))
        self.button1.pack()

        self.button2 = ttk.Button(text="test2", command=self.reset)
        self.button2.pack()

#        self.labelvar = tk.StringVar()
#        self.labelvar.set("text")
        self.text = ttk.Label(text="t")  # textvariable=self.labelvar)
        self.text.pack()

        self.progress = ttk.Progressbar(self, orient="horizontal",
                                        length=200, mode="determinate")
        self.progress.pack()

        self.bytes = 0
        self.maxbytes = 0

        self.text.config(text="starting")

        self.reset()

    def reset(self):
        self.progress["value"] = 0
        self.maxbytes = 50000
        self.progress["maximum"] = 50000
        self.text.config(text="starting")
        self.progress["value"]=25000

    def start(self):
        self.read_bytes()

    def test(self, wid):
#        self.labelvar.set("new " + str(wid))
        self.text.config(text="new " + str(wid))

    def read_bytes(self):
        '''simulate reading 500 bytes; update progress bar'''
        self.bytes += 500
        self.progress["value"] = self.bytes
        if self.bytes < self.maxbytes:
            # read more bytes after 100 ms
#            self.labelvar.set(str(self.bytes))
            # self.text['text'] = str(self.bytes)
            self.text.config(text=str(self.bytes))
            print(str(self.bytes), end=" ")
            self.after(100, self.read_bytes)

app = SampleApp()
app.mainloop()