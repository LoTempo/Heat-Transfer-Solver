from tkinter import *
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import numpy as np
from tkinter import ttk

FONT_NAME = "Arial"
FONT_SIZE = 12
FONT = (FONT_NAME, FONT_SIZE)

def center_window(window):
    window.update_idletasks()
    width = window.winfo_width()
    height = window.winfo_height()
    screen_width = window.winfo_screenwidth()
    screen_height = window.winfo_screenheight()

    x = (screen_width // 2) - (width // 2)
    y = (screen_height // 2) - (height // 2)

    window.geometry(f"{width}x{height}+{x}+{y}")

root = Tk()

root.geometry('1000x700')
center_window(root)
root.configure(bg='#DCEAF2')
root.title('Heat Transfer Solver v.1.0.0')

style = ttk.Style()
style.theme_use('clam') #'alt', 'clam', 'default', 'classic'
style.configure("TLabel", font=FONT)
style.configure("TButton", font=FONT)
style.configure("TEntry", font=FONT)

frame = Frame(root, bg='#DCEAF2')
frame.place(relwidth=1, relheight=1,y=5)

canvas2 = Canvas(frame, height=600, width=600)
canvas2.grid(row=1, column=3, rowspan=18,padx=35)

line_f = None
line_s = None
T_f = None
T_s = None
fig = None
dt = None
time_text = None
label = None

def render_formula(master, formula, row, column, bg='#DCEAF2'):
    fig, ax = plt.subplots(figsize=(2, 0.3), facecolor=bg)
    ax.set_facecolor(bg)
    ax.text(0.5, 0.5, formula, fontsize=14, ha='center', va='center', transform=ax.transAxes)
    ax.axis("off")
    fig.tight_layout()

    canvas = FigureCanvasTkAgg(fig, master=master)
    canvas_widget = canvas.get_tk_widget()
    canvas_widget.grid(row=row, column=column, sticky='e', padx=1, pady=1)

def btn_click():
    global line_f, line_s, T_f, T_s, fig, dt, time_text, label

    L = float(entries["L_input"].get())
    T_ext1 = float(entries["T_ext1_input"].get())
    T_ext2 = float(entries["T_ext2_input"].get())
    T_init_f = float(entries["T_init_f_input"].get())
    K_ff = float(entries["K_ff_input"].get())
    K_ss = float(entries["K_ss_input"].get())
    K_karks = float(entries["K_karks_input"].get())
    h_sf = float(entries["h_sf_input"].get())
    h_ext1 = float(entries["h_ext1_input"].get())
    h_ext2 = float(entries["h_ext2_input"].get())
    rho_f = float(entries["rho_f_input"].get())
    rho_s = float(entries["rho_s_input"].get())
    co_f = float(entries["co_f_input"].get())
    co_s = float(entries["co_s_input"].get())
    epsilon = float(entries["epsilon_input"].get())
    u = float(entries["u_input"].get())
    N_x = int(entries["N_x_input"].get())
    N_t = int(entries["N_t_input"].get())
    dt = float(entries["dt_input"].get())
    T_in = T_ext1
    T_init_s = T_init_f
    dx = L / N_x

    T_f = np.ones((N_t, N_x)) * T_init_f
    T_s = np.ones((N_t, N_x)) * T_init_s

    T_f[:, 0] = T_in

    for i in range(1, N_t):
        for j in range(1, N_x - 1):
            T_s[i, j] = (1 / ((1 - epsilon) *rho_s*co_s )) * (
                    ((T_s[i - 1, j - 1] - 2 * T_s[i - 1, j] + T_s[i - 1, j + 1]) / dx ** 2) * K_ss + h_sf * (
                        T_f[i - 1, j] - T_s[i - 1, j])) * dt + T_s[i - 1, j]

            T_f[i, j] = (((K_ff * ((T_f[i - 1, j - 1] - 2 * T_f[i - 1, j] + T_f[i - 1, j + 1]) / dx ** 2) + h_sf * (
                        T_s[i - 1, j] - T_f[i - 1, j])) / (
                                  epsilon * rho_f*co_f)) - u * (T_f[i - 1, j] - T_f[i - 1, j - 1]) / dx) * dt + T_f[
                            i - 1, j]

        T_s[i, 0] = (h_ext1 * dx * T_ext1 + K_karks * T_s[i, 1]) / (K_karks + h_ext1 * dx)

        T_s[i, -1] = (K_karks * T_s[i, -2] + h_ext2 * dx * T_ext2) / (K_karks + h_ext2 * dx)

        T_f[i, -1] = T_f[i, -2]

    x = np.linspace(0, L, N_x)

    fig, ax = plt.subplots(figsize=(6, 4))
    fig.patch.set_facecolor('#DCEAF2')
    fig.tight_layout(pad=2)

    line_f, = ax.plot(x, T_f[N_t-1, :], linestyle='--', label=f"Жидкость")
    line_s, = ax.plot(x, T_s[N_t-1, :], linestyle='-', label=f"Твердое тело")

    ax.set_xlabel(r"$x, м$", fontsize=14)
    ax.set_ylabel(r"$T, ^{\circ}C$", fontsize=14)
    ax.set_ylim(min(T_f.min(), T_s.min()), max(T_f.max(), T_s.max()))
    ax.set_xlim(0, L)

    ax.legend()
    time_text = ax.text(0.72, 1.04, '', transform=ax.transAxes,fontsize=12)

    for widget in canvas2.winfo_children():
        widget.destroy()

    canvas = FigureCanvasTkAgg(fig, master=canvas2)
    canvas.draw()
    canvas.get_tk_widget().pack()

    label = Label(frame, text=f"Шаг по времени",font=(12), bg='#DCEAF2')
    label.grid(row=19, column=3)

    slider = ttk.Scale(frame, from_=0, to=(int(entries["N_t_input"].get())-1), orient=HORIZONTAL, length=500)
    slider.grid(row=20, column=3, padx=50, pady=5)

    slider.bind("<Motion>", lambda event: update_graph(slider.get()))

def update_graph(val):
    global line_f, line_s, T_f, T_s, fig, dt, time_text, label

    t = int(float(val))
    line_f.set_ydata(T_f[t, :])
    line_s.set_ydata(T_s[t, :])
    time_text.set_text(f'Время = {(t) * dt:.2f} с')
    fig.canvas.draw_idle()

    label.config(text=f"Шаг по времени: {t}")

fields = [
    (r"$L, м$", "L_input", "0.04"),
    (r"$T_{1}, ^{\circ}C$", "T_ext1_input", "50.0"),
    (r"$T_{2}, ^{\circ}C$", "T_ext2_input", "0.0"),
    (r"$T_{0}, ^{\circ}C$", "T_init_f_input", "0.0"),
    (r"$\lambda_{\text{эфф.ж}}, Вт/(м\cdot ^{\circ}C)$", "K_ff_input", "0.296"),
    (r"$\lambda_{\text{эфф.т}}, Вт/(м\cdot ^{\circ}C)$", "K_ss_input", "0.059"),
    (r"$\lambda_{\text{т}}, Вт/(м\cdot ^{\circ}C)$", "K_karks_input", "0.375"),
    (r"$\alpha_{\text{т.ж}}, Вт/(м^2\cdot ^{\circ}C)$", "h_sf_input", "500.0"),
    (r"$\alpha_{1}, Вт/(м^2\cdot ^{\circ}C)$", "h_ext1_input", "1000.0"),
    (r"$\alpha_{2}, Вт/(м^2\cdot ^{\circ}C)$", "h_ext2_input", "10.0"),
    (r"$\rho_{\text{ж}}, кг/м^3$", "rho_f_input", "1000.0"),
    (r"$\rho_{\text{т}}, кг/м^3$", "rho_s_input", "1412.0"),
    (r"$c_{p.ж}, Дж/(кг\cdot ^{\circ}C)$", "co_f_input", "4200.0"),
    (r"$c_{p.т}, Дж/(кг\cdot ^{\circ}C)$", "co_s_input", "800.0"),
    (r"$\phi$", "epsilon_input", "0.8"),
    (r"$u, м/с$", "u_input", "0.0002"),
    (r"$N$", "N_x_input", "50"),
    (r"$M$", "N_t_input", "5000"),
    (r"$\Delta t, с$", "dt_input", "0.01")
]

entries = {}

for i, (formula, var_name, default_value) in enumerate(fields):
    render_formula(frame, formula, row=i + 1, column=0)
    entry = ttk.Entry(frame)
    entry.grid(row=i + 1, column=1, sticky='ew', padx=1, pady=1)
    entry.insert(0, default_value)
    entries[var_name] = entry

btn = ttk.Button(frame, text='Обновить график', command=btn_click)
btn.grid(row=len(fields) + 1, column=0, columnspan=2, pady=10)

root.mainloop()
