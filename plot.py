import numpy as np
import matplotlib.pyplot as plt


def single_plot():
    data0 = np.loadtxt("u0.txt").transpose()
    r0 = data0[0]
    u0 = data0[1]
    
    data1 = np.loadtxt("u1.txt").transpose()
    r1 = data1[0]
    u1 = data1[1]
    
    data2 = np.loadtxt("u2.txt").transpose()
    r2 = data2[0]
    u2 = data2[1]
    
    data3 = np.loadtxt("u3.txt").transpose()
    r3 = data3[0]
    u3 = data3[1]

    plt.figure(figsize=(9, 5))
    plt.plot(r0, u0, linewidth=2)
    plt.plot(r1, u1, linewidth=2)
    plt.plot(r2, u2, linewidth=2)
    plt.plot(r3, u3, linewidth=2)
    plt.xlabel(r"Ось $\eta = y\sqrt{\frac{U_1}{2\nu x}}$", labelpad=-4, fontsize=12)
    plt.ylabel(r"$f\ '(\eta)$", fontsize=12)
    plt.title(r"Численное решение краевой задачи. f'(0) = 0.107852")
    plt.grid(True)
    plt.grid(which="major", linewidth=1)
    plt.grid(which="minor", linewidth=0.4)
    plt.minorticks_on()

    plt.show()

single_plot()
