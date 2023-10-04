from matplotlib import pyplot as plt
import random as rand
import math
import numpy as np

# f(x) = x^4 + x^2 + x + 1, x э [-1; 0]

# 1. Метод перебора            (+)
# 2. Метод поразрядного поиска (+)
# 3. Метод дихотомии           (+)
# 4. Метод золотого сечения    (-)
# 5. Метод парабол             (-)
# 6. Метод средней точки       (-)
# 7. Метод хорд                (-)
# 8. Метод Ньютона             (-)

a = -1
b = 0
eps = 0.01

def f(x):
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
    
    return

x = np.linspace(a, b, 1000, endpoint=True)
y = f(x)

method_perebora(a, b, eps)
method_porazryadnogo_poiska(a, b, eps, 0.25)
method_dihtomia(a, b, eps, eps / 2)

plt.title('f(x) = x^4 + x^2 + x + 1')
plt.xlabel('X')
plt.ylabel('Y')
plt.grid(visible=True)
plt.plot(x, y)
plt.show()
