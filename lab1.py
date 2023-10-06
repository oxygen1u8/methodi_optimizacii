from matplotlib import pyplot as plt
import math
import numpy as np
import sympy as sp

# f(x) = x^4 + x^2 + x + 1, x э [-1; 0]

# 1. Метод перебора            (+)
# 2. Метод поразрядного поиска (+)
# 3. Метод дихотомии           (+)
# 4. Метод золотого сечения    (?)
# 5. Метод парабол             (+)
# 6. Метод средней точки       (+)
# 7. Метод хорд                (+)
# 8. Метод Ньютона             (-)

a = -1
b = 0
eps = 0.01

def f(x):
    return x**4 + x**2 + x + 1   

def f_sym(x: sp.Symbol):
    return x**4 + x**2 + x + 1

def print_res(res_x, res_y, steps):
    print('Min y =', res_y, 'at x =', res_x, 'steps =', steps)

def method_perebora(a, b, eps):
    n = math.ceil((b - a) / eps)
    i = 0
    x_i = a + i * (b - a) / n
    y_i = f(x_i)
    res_x = 0
    res_y = 0
    steps = 0
    while i < (n + 1):
        i += 1
        x_i1 = a + i * (b - a) / n
        y_i1 = f(x_i1)
        if y_i1 >= y_i:
            break
        x_i, y_i = x_i1, y_i1
        steps += 1
    res_x, res_y = x_i, y_i
    print_res(res_x, res_y, steps)
    return res_y

def method_porazryadnogo_poiska(a, b, eps, delta):
    x_i = a
    y_i = f(x_i)
    res_x = x_i
    res_y = y_i
    steps = 0
    while True:
        x_i1 = x_i + delta
        if x_i1 > b:
            break
        y_i1 = f(x_i1)
        if y_i1 > y_i:
            if abs(delta) < eps:
                res_x = x_i
                res_y = y_i
                break
            delta = -delta / 2
        x_i, y_i = x_i1, y_i1
        steps += 1
    print_res(res_x, res_y, steps)
    return res_y

def method_dihtomia(a, b, eps, delta):    
    res_x = 0
    res_y = 0
    steps = 0
    eps_n = (b - a) / 2
    while eps_n > eps:
        x1 = (b + a - delta) / 2
        x2 = (b + a + delta) / 2
        if f(x1) <= f(x2):
            b = x2
        else:
            a = x1
        eps_n = (b - a) / 2
        steps += 1
    res_x = (a + b) / 2
    res_y = f(res_x)
    print_res(res_x, res_y, steps)
    return res_y

def method_zolotogo_secheniya(a, b, eps):
    steps = 0
    res_x = 0
    res_y = 0
    eps_n = (b - a) / 2
    t = (math.sqrt(5) - 1) / 2

    x1 = a + (3 - math.sqrt(5)) / 2 * (b - a)
    x2 = a + b - x1
    # x2 = a + (math.sqrt(5) - 1) / 2 * (b - a)
    y1 = f(x1)
    y2 = f(x2)

    while abs(eps_n) > eps:
        steps += 1
        if y1 <= y2:
            b = x2
            x2 = x1
            y2 = y1
            x1 = a + t * (b - a)
            y1 = f(x1)
        else:
            a = x1
            x1 = x2
            y1 = y2
            x2 = b - t * (b - a)
            y2 = f(x2)
        eps_n = t * (b - a)

    res_x = (a + b) / 2
    res_y = f(res_x)
    print_res(res_x, res_y, steps)
    return res_y

def method_parabol(a, b, eps, delta):
    res_x = 0
    res_y = 0
    steps = 0

    x1 = a + delta
    x2 = x1 + delta
    x3 = x2 + delta

    x_values = [x1, x2, x3]

    x_prev = 0xffffffff

    while True:
        steps += 1
    
        y1 = f(x_values[0])
        y2 = f(x_values[1])
        y3 = f(x_values[2])

        a1 = (y2 - y1) / (x_values[1] - x_values[0])
        a2 = 1 / (x_values[2] - x_values[1]) * ((y3 - y1) / (x_values[2] - x_values[0]) - (y2 - y1) / (x_values[1] - x_values[0]))

        x_ = 1/2  * (x_values[0] + x_values[1] - a1 / a2)

        if x_prev != 0xffffffff:
            if abs(x_prev - x_) <= eps:
                res_x = x_
                res_y = f(x_)
                break

        # y_ = f(x_)
        for i in range(len(x_values)):
            if x_ <= x_values[i]:
                for j in range(i - 1):
                    x_values[j] = x_values[j + 1]
                x_values[i - 1] = x_
                break
        x_prev = x_
    
    print_res(res_x, res_y, steps)
    return res_y

def method_sredney_tochki(a, b, eps):
    res_x = res_y = steps = 0
    x = sp.symbols('x')
    func = f_sym(x)

    while True:
        steps += 1

        x_ = (a + b) / 2
        f_diff = sp.diff(func, x)
        diff_value = f_diff.subs(x, x_)

        if abs(diff_value) <= eps:
            res_x = x_
            res_y = f(res_x)
            break
        
        if diff_value > 0:
            b = x_
        else:
            a = x_

    print_res(res_x, res_y, steps)
    return res_y

def method_hord(a, b, eps):
    res_x = res_y = steps = 0
    x = sp.symbols('x')
    func = f_sym(x)
    diff_f = sp.diff(func, x)

    while True:
        steps += 1
        tmp1 = float(diff_f.subs(x, a))
        tmp2 = float(diff_f.subs(x, b))
        x_ = a - tmp1 / (tmp1 - tmp2) * (a - b)
        if abs(diff_f.subs(x, x_)) <= eps:
            res_x = x_
            res_y = f(res_x)
            break
        if diff_f.subs(x, x_) > 0:
            b = x_ 
        else:
            a = x_
    
    print_res(res_x, res_y, steps)
    return res_y

x = np.linspace(a, b, 1000, endpoint=True)
y = f(x)

method_perebora(a, b, eps)
method_porazryadnogo_poiska(a, b, eps, 0.25)
method_dihtomia(a, b, eps, eps / 2)
method_zolotogo_secheniya(a, b, eps)
method_parabol(a, b, eps, 0.25)
method_sredney_tochki(a, b, eps)
method_hord(a, b, eps)

plt.title('f(x) = x^4 + x^2 + x + 1')
plt.xlabel('X')
plt.ylabel('Y')
plt.grid(visible=True)
plt.plot(x, y)
plt.show()
