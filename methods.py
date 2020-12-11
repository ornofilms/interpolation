import numpy as np
import matplotlib.pyplot as plt


def lagrange_method(x, y, t, n):
    num, denom = 1, 1
    lagrange = 0
    for i in range(n):
        for j in range(n):
            if i != j:
                num *= t - x[j]
                denom *= x[i] - x[j]
        lagrange += y[i] * num / denom
        num, denom = 1, 1
    return lagrange


def newtone_diff(x, y, n):
    f_diff = [0 for i in range(n)]
    f_diff[0] = [i for i in y]
    for i in range(n):
        print(f"f(x{i}) = {x[i]}")
    for i in range(1, n):
        f_diff[i] = []
        for j in range(n - i):
            f_diff[i].append((f_diff[i - 1][j] - f_diff[i - 1][j + 1]) / (x[j] - x[j + i]))
            if (i == 1):
                print(f"f(x{j};x{j + i}) = (f(x{j}) - f(x{j + 1})) / (x{j} - x{j + i}) = "
                      f"({f_diff[i - 1][j]} - {f_diff[i - 1][j + 1]}) / ({x[j]} - {x[j + i]}) = {f_diff[i][-1]}")
            else:
                print(f"f(x{j};x{j + i}) = (f(x{j};x{j + i - 1}) - f(x{j + 1};x{j + i})) / (x{j} - x{j + i}) = "
                      f"({f_diff[i - 1][j]} - {f_diff[i - 1][j + 1]}) / ({x[j]} - {x[j + i]}) = {f_diff[i][-1]}")
    return f_diff


def newtone_method(x, y, t, n, diff_storage):
    fx = y[0]
    product = 1
    for i in range(1, n):
        for j in range(i):
            product *= (t - x[j])
        fx += product * diff_storage[i][0]
        product = 1
    return fx


def gauss_slae(matrix, k):
    # forward stroke
    for iter in range(k):
        imax = iter
        jmax = 0
        for i in range(iter, k+1): # find maximum modulo value in all matrix
            for j in range(k+1):
                if abs(matrix[i][j]) > abs(matrix[imax][jmax]):
                    imax = i
                    jmax = j
        for j in range(k+2):      # swap row witn max value and current-stap row
            matrix[iter][j], matrix[imax][j] = matrix[imax][j], matrix[iter][j]
        c_max = matrix[iter][jmax]
        for j in range(k + 2):
            if c_max != 0:
                matrix[iter][j] /= c_max
        for i in range(iter + 1, k + 1):
            c_max = matrix[i][jmax]
            if c_max != 0:
                for j in range(k + 2):
                    matrix[i][j] -= matrix[iter][j] * c_max
    # reverse stroke
    solve = {}
    for i in reversed(range(k + 1)):
        b = matrix[i][-1]
        for j in range(k + 1):
            if j in solve.keys():
                b -= solve[j] * matrix[i][j]
        for j in range(k + 1):
            if matrix[i][j] != 0 and j not in solve.keys():
                solve[j] = b / matrix[i][j]
    return solve


def run_through_slae(a, b, c, d, n):
    y = []
    alpha = []
    beta = []
    # forward stroke
    for i in range(n):
        if i == 0:
            y.append(b[i])
            alpha.append(- c[i] / y[i])
            beta.append(d[i] / y[i])
        else:
            y.append(b[i] + a[i] * alpha[i - 1])
            beta.append((d[i] - a[i] * beta[i - 1]) / y[i])
            if i < n - 1:
                alpha.append(- c[i] / y[i])
    x = {}
    # reverse stroke
    for i in reversed(range(n)):
        if i == n - 1:
            x[i] = beta[i]
        else:
            x[i] = alpha[i] * x[i + 1] + beta[i]
    return x


def koeff_a_least_square(x, y, n, k):
    c = [0] * (2 * k + 1)
    d = [0] * (k + 1)
    for m in range(len(c)):
        for i in range(n):
            c[m] += x[i] ** m
    for j in range(len(d)):
        for i in range(n):
            d[j] += y[i] * x[i] ** j
    matrix = []
    for i in range(k + 1):  # forming a matrix from c values
        matrix.append([c[j + i] for j in range(k + 1)])
        matrix[i].append(d[i])
    return gauss_slae(matrix, k)


def least_square_method(a, t, k):
    p = 0
    for i in reversed(range(k + 1)):
        p += a[i] * t**i
    return p


def koeff_cubic_spline(x, y, n):
    koeff = {'a': {}, 'b': {}, 'c': {}, 'd': {}}  # a-d koeff of cubic spline
    a = []                                        # a-d for run-through method
    b = []
    c = []
    d = []
    for i in range(2, n):
        j = i - 2
        hi_1 = abs(x[i - 2] - x[i - 1])
        hi = abs(x[i - 1] - x[i])
        a.append(hi_1 if j != 0 else 0)
        b.append(2 * (hi_1 + hi))
        c.append(hi if i != n - 1 else 0)
        d.append(3 * ((y[i] - y[i - 1]) / hi - (y[i - 1] - y[i - 2]) / hi_1))
    c_koeff = run_through_slae(a, b, c, d, n - 2)
    koeff['c'][1] = 0
    for i in c_koeff.keys():
        koeff['c'][i + 2] = c_koeff[i]
    koeff['c'][n] = 0
    for i in range(1, n):
        h = abs(x[i - 1] - x[i])
        koeff['a'][i] = y[i - 1]
        koeff['b'][i] = (y[i] - y[i - 1]) / h - ((koeff['c'][i + 1] + 2 * koeff['c'][i]) * h) / 3
        koeff['d'][i] = (koeff['c'][i + 1] - koeff['c'][i]) / (3 * h)
    return koeff


def cubic_spline_method(x, koeff, t, n):
    t_round = round(t, 2)
    for i in range(1, n):
        if x[i - 1] <= t_round <= x[i]:
            return koeff['a'][i] + koeff['b'][i] * (t_round - x[i - 1]) + \
                   koeff['c'][i] * (t_round - x[i - 1])**2 + koeff['d'][i] * (t_round - x[i - 1])**3
