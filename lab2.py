from matplotlib import pyplot as plt
import math
import numpy as np
import mpmath as mp
import sympy as sp
import pandas as pd

# f(x) = 64(x1)^2 + 126x1x2 + 64(x2)^2 - 10x1 + 30x2 + 13, x э [-1; 0]

# 1. Метод наискорейшего спуска                (-)
# 2. Метод сопряженных градиентов              (-)
# 3. Метод Ньютона                             (-)
# 4. Метод правильного симплекса               (-)
# 5. Метод циклического покоординатного спуска (-)
# 6. Метод Хука-Дживса                         (-)
# 7. Метод случайного поиска                   (-)

# Овражная функция: f(x1, x2) = (x1)^2 + a(x2)^2

a = [1, 250, 1000]
eps = [10e-4, 10e-6]

def o_f(a, x1: sp.Symbol, x2: sp.Symbol):
    return x1**2 + a*(x2**2)

def f(x1: sp.Symbol, x2: sp.Symbol):
    return 64 * (x1 ** 2) + 126 * x1 * x2 + 64 * (x2 ** 2) - 10 * x1 + 30 * x2 + 13

# Задание 1

def method_naisk_spuska(func, x, x0, eps, b=np.array([[0], [0]])):
    steps = 0

    x1 = x[0][0]
    x2 = x[1][0]

    diff1 = sp.diff(func, x1)
    diff2 = sp.diff(func, x2)

    while True:
        x10 = x0[0][0]
        x20 = x0[1][0]
        H = np.array(
            [
                [float(sp.diff(diff1, x1).subs(x1, x10).subs(x2, x20)), float(sp.diff(diff1, x2).subs(x1, x10).subs(x2, x20))],
                [float(sp.diff(diff2, x1).subs(x1, x10).subs(x2, x20)), float(sp.diff(diff2, x2).subs(x1, x10).subs(x2, x20))]
            ]
        )
        n_f = np.array([[float(diff1.subs(x1, x10).subs(x2, x20))], [float(diff2.subs(x1, x10).subs(x2, x20))]])
        tmp = 0
        for i in n_f:
            tmp += i[0] ** 2
        norm_n_f = math.sqrt(tmp)
        if norm_n_f <= eps:
            res_x = np.array(x0)
            res_f = float(func.subs(x1, x10).subs(x2, x20))
            break
        p = -n_f
        alpha = -sum([i[0] ** 2 for i in n_f])
        tmp = H.dot(p)
        tmp = sum(tmp[i][0] * p[i][0] for i in range(len(p)))
        alpha /= tmp
        tmp = alpha * n_f
        x0 = x0 + tmp
        steps += 1

    return res_x, res_f, steps

def method_sopr_grad(func, eps):
    res_x1, res_x2, res_f, steps = 0, 0, 0, 0
    return res_x1, res_x2, res_f, steps

def method_newton(func, eps):
    res_x1, res_x2, res_f, steps = 0, 0, 0, 0
    return res_x1, res_x2, res_f, steps

def method_true_simplex(func, eps):
    res_x1, res_x2, res_f, steps = 0, 0, 0, 0
    return res_x1, res_x2, res_f, steps

def method_acycle_coord_spusk(func, eps):
    res_x1, res_x2, res_f, steps = 0, 0, 0, 0
    return res_x1, res_x2, res_f, steps

def method_huka_djives(func, eps):
    res_x1, res_x2, res_f, steps = 0, 0, 0, 0
    return res_x1, res_x2, res_f, steps

def method_sluch_poiska(func, eps):
    res_x1, res_x2, res_f, steps = 0, 0, 0, 0
    return res_x1, res_x2, res_f, steps

# Задание 2

x1 = sp.symbols('x1')
x2 = sp.symbols('x2')

a1 = []
a2 = []
a3 = []
a4 = []
a5 = []
a6 = []
a7 = []
a8 = []
a9 = []

x = np.array([[x1], [x2]])
x0 = np.array([[3], [7]])

for i in a:
    for j in eps:
        f_ = o_f(i, x1, x2)
        ff = f(x1, x2)
        a1.append(i)
        a2.append(j)
        res_x1, res_f1, steps1 = method_naisk_spuska(f_, x, x0, j)
        res_x2, res_f2, steps2 = method_naisk_spuska(ff, x, x0, j)
        a3.append(steps1)
        # res_x12, res_x22, res_f2, steps2 = method_sopr_grad(f, j)
        a4.append(steps2)
        # res_x13, res_x23, res_f3, steps3 = method_newton(f, j)
        # a5.append(steps3)
        # res_x14, res_x24, res_f4, steps4 = method_true_simplex(f, j)
        # a6.append(steps4)
        # res_x15, res_x25, res_f5, steps5 = method_acycle_coord_spusk(f, j)
        # a7.append(steps5)
        # res_x16, res_x26, res_f6, steps6 = method_huka_djives(f, j)
        # a8.append(steps6)
        # res_x17, res_x27, res_f7, steps7 = method_sluch_poiska(f, j)
        # a9.append(steps7)

df = pd.DataFrame({
    'a' : a1, 'Точность' : a2,
    'Метод наискорейшего спуска' : a3,
    'Метод наискорейшего спуска 2' : a4,
    # 'Метод сопряженных градиентов' : a4, 'Метод Ньютона' : a5,
    # 'Метод правильного симплекса' : a6,
    # 'Метод циклического покоординатного спуска' : a7,
    # 'Метод Хука-Дживса' : a8,
    # 'Метод случайного поиска' : a9
})

print(df)
