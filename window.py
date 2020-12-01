from tkinter import *
import tkinter.ttk as ttk
from methods import lagrange_method, newtone_method, newtone_diff, least_square_method, cubic_spline_method, \
    koeff_a_least_square, koeff_cubic_spline
from matplotlib import pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
from PIL import ImageTk, Image

bg_color = "#F6F1E5"
bg_caution = "#F6EAEA"
bg_button = "#EBE3BF"
font = ("Arial", 11)


def px_to_inch(value):
    return .0104 * value


def graph_output(event):
    try:
        values_x = [float(x0.get()), float(x1.get()), float(x2.get()), float(x3.get()), float(x4.get())]
        values_y = [float(y0.get()), float(y1.get()), float(y2.get()), float(y3.get()), float(y4.get())]
    except (ValueError, UnboundLocalError):
        caution["text"] = "Введите данные\nкорректно"
        caution["fg"] = "red"
        return
    caution["text"] = "OK"
    caution["fg"] = "green"
    fig, ax = plt.subplots(figsize=(px_to_inch(512), px_to_inch(383.2)))
    step = .1
    num_points = 5
    if var.get() == 0:
        lagrange_x = np.arange(values_x[0], values_x[-1] + step, step)
        lagrange_y = [lagrange_method(values_x, values_y, x, num_points) for x in lagrange_x]
        ax.plot(lagrange_x, lagrange_y)
    elif var.get() == 1:
        newtone_x = np.arange(values_x[0], values_x[-1] + step, step)
        diff_storage = newtone_diff(values_x, values_y, num_points)
        newtone_y = [newtone_method(values_x, values_y, x, num_points, diff_storage) for x in newtone_x]
        ax.plot(newtone_x, newtone_y, color='g')
    elif var.get() == 5:
        spline_x = np.arange(values_x[0], values_x[-1] + step, step)
        koeff = koeff_cubic_spline(values_x, values_y, num_points)
        spline_y = [cubic_spline_method(values_x, koeff, x, num_points) for x in spline_x]
        ax.plot(spline_x, spline_y, color='y')
    else:
        a = koeff_a_least_square(values_x, values_y, num_points, var.get() - 1)
        smooth_x = np.arange(values_x[0], values_x[-1] + step, step)
        smooth_y = [least_square_method(a, x, var.get() - 1) for x in smooth_x]
        ax.plot(smooth_x, smooth_y, color='m')

    ax.plot(values_x, values_y, 'r')
    ax.scatter(values_x, values_y, color='r', s=32, marker='o')
    canvas = FigureCanvasTkAgg(fig, master=main)
    canvas.get_tk_widget().place(x=78, y=263)


main = Tk()
main.geometry('687x708+0+0')
main.title('lab2')
main.iconbitmap('logo.ico')

photo = PhotoImage(file=r"C:\Users\User\PycharmProjects\interpolation\back_main.png")
background = Label(main, image=photo)
background.pack()

label_x = Label(main, text="x", bg=bg_color, font=font)
label_y = Label(main, text="y", bg=bg_color, font=font)
label_x.place(x=222, y=18)
label_y.place(x=222, y=52)

x0 = Entry(main, width=3, justify=CENTER)
x1 = Entry(main, width=3, justify=CENTER)
x2 = Entry(main, width=3, justify=CENTER)
x3 = Entry(main, width=3, justify=CENTER)
x4 = Entry(main, width=3, justify=CENTER)
y0 = Entry(main, width=3, justify=CENTER)
y1 = Entry(main, width=3, justify=CENTER)
y2 = Entry(main, width=3, justify=CENTER)
y3 = Entry(main, width=3, justify=CENTER)
y4 = Entry(main, width=3, justify=CENTER)
x0.place(x=265, y=20)
x1.place(x=311, y=20)
x2.place(x=357, y=20)
x3.place(x=403, y=20)
x4.place(x=449, y=20)
y0.place(x=265, y=54)
y1.place(x=311, y=54)
y2.place(x=357, y=54)
y3.place(x=403, y=54)
y4.place(x=449, y=54)

caution = Label(main, text="", fg="green", bg=bg_caution)
caution.place(x=553, y=25)

interpolation = Label(main, text="Интерполяционный многочлен\nс помощью формулы", bg=bg_color, font=font)
interpolation.place(x=111, y=99)
spline = Label(main, text="Сглаживающий многочлен\nпо методу\nнаименьших квадратов", bg=bg_color, font=font)
spline.place(x=448, y=94)

var = IntVar()
var.set(0)
var_lagrange = Radiobutton(main, text="Лагранжа", variable=var, value=0, bg=bg_color, font=font)
var_newtone = Radiobutton(main, text="Ньютона", variable=var, value=1, bg=bg_color, font=font)
var_pow1 = Radiobutton(main, text="1 степени", variable=var, value=2, bg=bg_color, font=font)
var_pow2 = Radiobutton(main, text="2 степени", variable=var, value=3, bg=bg_color, font=font)
var_pow3 = Radiobutton(main, text="3 степени", variable=var, value=4, bg=bg_color, font=font)
var_spline = Radiobutton(main, text="Кубический сплайн дефекта 1", variable=var, value=5, bg=bg_color, font=font)
var_lagrange.place(x=105, y=137)
var_newtone.place(x=219, y=137)
var_pow1.place(x=489, y=150)
var_pow2.place(x=489, y=172)
var_pow3.place(x=489, y=194)
var_spline.place(x=82, y=189)

run = Button(main, text="Начать", width=15, height=1, bg=bg_button, font=font)
run.bind("<Button-1>", graph_output)
run.place(x=291, y=230)

ex = Button(main, text="Выход", width=15, height=1, bg=bg_button, font=font)
ex.bind("<Button-1>", exit)
ex.place(x=538, y=664)

main.mainloop()
