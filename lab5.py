from random import randint
from numpy import *
import numpy
from scipy.stats import t, f
import math


class Critical_values:
    @staticmethod
    def get_cohren_value(size_of_selections, qty_of_selections, significance):
        from _pydecimal import Decimal
        from scipy.stats import f
        size_of_selections += 1
        partResult1 = significance / (size_of_selections - 1)
        params = [partResult1, qty_of_selections, (size_of_selections - 1 - 1) * qty_of_selections]
        fisher = f.isf(*params)
        result = fisher / (fisher + (size_of_selections - 1 - 1))
        return Decimal(result).quantize(Decimal('.0001')).__float__()

    @staticmethod
    def get_student_value(f3, significance):
        from _pydecimal import Decimal
        from scipy.stats import t
        return Decimal(abs(t.ppf(significance / 2, f3))).quantize(Decimal('.0001')).__float__()

    @staticmethod
    def get_fisher_value(f3, f4, significance):
        from _pydecimal import Decimal
        from scipy.stats import f
        return Decimal(abs(f.isf(significance, f4, f3))).quantize(Decimal('.0001')).__float__()


cr = Critical_values()


def dob(*args):
    res = [1 for _ in range(len(args[0]))]
    for i in range(len(args[0])):
        for j in args:
            res[i] *= j[i]
    return res


def getcolumn(arr, n):
    return [i[n] for i in arr]


m = 3
p = 0.95
bzn, yzn, beta = [], [], []
rows = N = 8
x1_min, x1_max = -20, 30
x2_min, x2_max = -20, 40
x3_min, x3_max = -20, -10
x_avarage_max = (x1_max + x2_max + x3_max) / 3
x_avarage_min = (x1_min + x2_min + x3_min) / 3
Ymax = 200 + x_avarage_max
Ymin = 200 + x_avarage_min

# матриця кодованих значень х
matrix_x_cod_for4 = [[+1, -1, -1, -1],
                     [+1, -1, +1, +1],
                     [+1, +1, -1, +1],
                     [+1, +1, +1, -1]]

matrix_x_for4 = [[x1_min, x2_min, x3_min],
                 [x1_min, x2_max, x3_max],
                 [x1_max, x2_min, x3_max],
                 [x1_max, x2_max, x3_min]]
matrix_x_for4 = numpy.array(matrix_x_for4)

# матриця кодованих значень х
matrix_x_cod = [[+1, -1, -1, -1, +1, +1, +1, -1],
                [+1, -1, -1, +1, +1, -1, -1, +1],
                [+1, -1, +1, -1, -1, +1, -1, +1],
                [+1, -1, +1, +1, -1, -1, +1, -1],
                [+1, +1, -1, -1, -1, -1, +1, +1],
                [+1, +1, -1, +1, -1, +1, -1, -1],
                [+1, +1, +1, -1, +1, -1, -1, -1],
                [+1, +1, +1, +1, +1, +1, +1, +1]]

# матриця значень х
matrix_x = [[1, x1_min, x2_min, x3_min, x1_min * x2_min, x1_min * x3_min, x2_min * x3_min, x1_min * x2_min * x3_min],
            [1, x1_min, x2_min, x3_max, x1_min * x2_min, x1_min * x3_max, x2_min * x3_max, x1_min * x2_min * x3_max],
            [1, x1_min, x2_max, x3_min, x1_min * x2_max, x1_min * x3_min, x2_max * x3_min, x1_min * x2_max * x3_min],
            [1, x1_min, x2_max, x3_max, x1_min * x2_max, x1_min * x3_max, x2_max * x3_max, x1_min * x2_max * x3_max],
            [1, x1_max, x2_min, x3_min, x1_max * x2_min, x1_min * x3_min, x2_min * x3_min, x1_min * x2_min * x3_min],
            [1, x1_max, x2_min, x3_max, x1_max * x2_min, x1_max * x3_max, x2_min * x3_max, x1_max * x2_min * x3_max],
            [1, x1_max, x2_max, x3_min, x1_max * x2_max, x1_max * x3_min, x2_max * x3_min, x1_max * x2_max * x3_min],
            [1, x1_max, x2_max, x3_max, x1_max * x2_max, x1_max * x3_max, x2_max * x3_max, x1_max * x2_max * x3_max]]

check = True
while check:

    # матриця рандомних значень у
    random_matrix_y = random.randint(Ymin, Ymax, size=(rows, m))


    # сума середніх значень відгуку функції за рядками
    def sum_rows(random_matrix_y):
        y = numpy.sum(random_matrix_y, axis=1) / m
        return y


    Yavg = sum_rows(random_matrix_y)


    def sum_columns(matrix_x_for4):
        mx = numpy.sum(matrix_x_for4, axis=0) / 4
        return mx


    mx = sum_columns(matrix_x_for4)


    # Нормовані коефіціенти рівняння регресії
    def sum_my(y1, y2, y3, y4):
        my = (y1 + y2 + y3 + y4) / 4
        return my


    my = sum_my(Yavg[0], Yavg[3], Yavg[5], Yavg[6])


    # Нормовані коефіціенти рівняння регресії
    def find_a(a, b, c, d):
        az = (a * Yavg[0] + b * Yavg[3] + c * Yavg[5] + d * Yavg[6]) / 4
        return az


    a1 = find_a(x1_min, x1_min, x1_max, x1_max)
    a2 = find_a(x2_min, x2_max, x2_min, x2_max)
    a3 = find_a(x3_min, x3_max, x3_max, x3_min)


    # Нормовані коефіціенти рівняння регресії
    def find_aa(a, b, c, d):
        aa = (a ** 2 + b ** 2 + c ** 2 + d ** 2) / 4
        return aa


    a11 = find_aa(x1_min, x1_min, x1_max, x1_max)
    a22 = find_aa(x2_min, x2_max, x2_min, x2_max)
    a33 = find_aa(x3_min, x3_max, x3_max, x3_min)

    # Нормовані коефіціенти рівняння регресії
    a12 = a21 = (x1_min * x2_min + x1_min * x2_max + x1_max * x2_min + x1_max * x2_max) / 4
    a13 = a31 = (x1_min * x3_min + x1_min * x3_max + x1_max * x3_max + x1_max * x3_min) / 4
    a23 = a32 = (x2_min * x3_min + x2_max * x3_max + x2_min * x3_max + x2_max * x3_min) / 4

    # Матриця для визначення коефіціентів регресії
    A = [[my, mx[0], mx[1], mx[2]], [a1, a11, a12, a13], [a2, a12, a22, a32], [a3, a13, a23, a33]]
    B = [[1, my, mx[1], mx[2]], [mx[0], a1, a12, a13], [mx[1], a2, a22, a32], [mx[2], a3, a23, a33]]
    C = [[1, mx[0], my, mx[2]], [mx[0], a11, a1, a13], [mx[1], a12, a2, a32], [mx[2], a13, a3, a33]]
    D = [[1, mx[0], mx[1], my], [mx[0], a11, a12, a1], [mx[1], a12, a22, a2], [mx[2], a13, a23, a3]]
    E = [[1, mx[0], mx[1], mx[2]], [mx[0], a11, a12, a13], [mx[1], a12, a22, a32], [mx[2], a13, a23, a33]]
    X = []


    # Коефіціенти регресії
    def coef_regr(a, b):
        b = linalg.det(a) / linalg.det(b)
        return b


    b0 = coef_regr(A, E)
    b1 = coef_regr(B, E)
    b2 = coef_regr(C, E)
    b3 = coef_regr(D, E)
    X.append(round(b0, 2))
    X.append(round(b1, 2))
    X.append(round(b2, 2))
    X.append(round(b3, 2))


    # Нормоване рівняння регресії
    def find_y_norm(a, b, c):
        y_norm = X[0] + X[1] * a + X[2] * b + X[3] * c
        return y_norm


    y_norm1 = find_y_norm(x1_min, x2_min, x3_min)
    y_norm2 = find_y_norm(x1_min, x2_max, x3_max)
    y_norm3 = find_y_norm(x1_max, x2_min, x3_max)
    y_norm4 = find_y_norm(x1_max, x2_max, x3_min)

    # Перевірка однорідності дисперсії за критерієм Кохрена
    # Пошук дисперсій по рядкам
    dispersion_y = [0, 0, 0, 0]

    for i in range(m):
        dispersion_y[0] += ((random_matrix_y[0][i] - Yavg[0]) ** 2) / m
        dispersion_y[1] += ((random_matrix_y[1][i] - Yavg[3]) ** 2) / m
        dispersion_y[2] += ((random_matrix_y[2][i] - Yavg[5]) ** 2) / m
        dispersion_y[3] += ((random_matrix_y[3][i] - Yavg[6]) ** 2) / m

    ajk = dispersion_y[0] + dispersion_y[1] + dispersion_y[2] + dispersion_y[3]

    Gp = 0
    if ajk == 0:
        m += 1
        print("Збільшуємо m на одиницю")
    else:
        Gp = max(dispersion_y) / (ajk)
        f1 = m - 1
        f2 = rows
        q = 1 - p
        Gt = Critical_values.get_cohren_value(f2, f1, q)
        if Gp <= Gt:
            print("Дисперсія однорідна")
            check = False
        else:
            m += 1
            print("Збільшуємо m на одиницю")

f1 = m - 1
f2 = rows
f3 = f1 * f2
Ft = cr.get_student_value(f3, q)
Sb = sum(dispersion_y) / rows
Sbetakvadr = Sb / (rows * m)
Sbeta = sqrt(Sb / (rows * m))


# Визначимо оцінки коефіціентів
def find_beta(a, b, c, d):
    beta = (Yavg[0] * a + Yavg[3] * b + Yavg[5] * c + Yavg[6] * d) / rows
    return beta


beta0 = find_beta(matrix_x_cod[0][0], matrix_x_cod[1][0], matrix_x_cod[2][0], matrix_x_cod[3][0])
beta1 = find_beta(matrix_x_cod[0][1], matrix_x_cod[1][1], matrix_x_cod[2][1], matrix_x_cod[3][1])
beta2 = find_beta(matrix_x_cod[0][2], matrix_x_cod[1][2], matrix_x_cod[2][2], matrix_x_cod[3][2])
beta3 = find_beta(matrix_x_cod[0][3], matrix_x_cod[1][3], matrix_x_cod[2][3], matrix_x_cod[3][3])


# Пошук коефіціента t
def find_t(a, b):
    t = a / b
    return t


t0 = find_t(beta0, Sbeta)
t1 = find_t(beta1, Sbeta)
t2 = find_t(beta2, Sbeta)
t3 = find_t(beta3, Sbeta)
t_list = [fabs(t0), fabs(t1), fabs(t2), fabs(t3)]
b_list = [b0, b1, b2, b3]

tbool = tuple(Ft < i for i in t_list)


# Запишемо рівняння з урахуванням критерію Стьюдента
def find_yj(a, b, c):
    yj = b_list[0] + b_list[1] * a + b_list[2] * b + b_list[3] * c
    return yj


yj1 = find_yj(x1_min, x2_min, x3_min)
yj2 = find_yj(x1_min, x2_max, x3_max)
yj3 = find_yj(x1_max, x2_min, x3_max)
yj4 = find_yj(x1_max, x2_max, x3_min)

# Перевірка умови за критерієм Фішера
d = tbool.count(True)  # кількість значимих коефіціентів
f1 = m - 1
f2 = rows
f4 = rows - d
f3 = f1 * f2
Sad = m * (((yj1 - Yavg[0]) ** 2 + (yj2 - Yavg[3]) ** 2 + (yj3 - Yavg[5]) ** 2 + (yj4 - Yavg[6]) ** 2)) / f4
Fp = Sad / Sbetakvadr
Fp = cr.get_fisher_value(f3, f4, q)
print("m = {}\np = {}\n".format(m, p))
print("Нормоване рівняння регресії y = {:.2f} + {:.2f}*x1 + {:.2f}*x2".format(X[0], X[1], X[2]))
print("{:.1f} + {:.1f} + {:.1f} + {:.1f} = {:.1f}".format(X[0], X[1] * x1_min, X[2] * x2_min, X[3] * x3_min, y_norm1))
print("{:.1f} + {:.1f} + {:.1f} + {:.1f} = {:.1f}".format(X[0], X[1] * x1_min, X[2] * x2_max, X[3] * x3_max, y_norm2))
print("{:.1f} + {:.1f} + {:.1f} + {:.1f} = {:.1f}".format(X[0], X[1] * x1_max, X[2] * x2_min, X[3] * x3_max, y_norm3))
print("{:.1f} + {:.1f} + {:.1f} + {:.1f} = {:.1f}".format(X[0], X[1] * x1_max, X[2] * x2_max, X[3] * x3_min, y_norm4))
print("\n")
print("Перевірка за Кохреном")
print("S²{y2}: ", round(dispersion_y[0], 2))
print("S²{y2}: ", round(dispersion_y[1], 2))
print("S²{y3}: ", round(dispersion_y[2], 2))
print("S²{y4}: ", round(dispersion_y[3], 2))
print("Gp: ", Gp)
print("\n")
print("Перевірка за Стьюдентом")
print("Sb²: {:.2f} \t\tS²(β): {:.2f} \t\tS(β): {:.2f}".format(Sb, Sbetakvadr, Sbeta))
print("β1: {:.2f} \t\t\tβ2: {:.2f} \t\tβ3: {:.2f} \t\tβ4: {:.2f}".format(beta0, beta1, beta2, beta3))
print("t0: {:.2f} \t\t\tt1: {:.2f} \t\t\tt2: {:.2f} \t\tt3: {:.2f}".format(t0, t1, t2, t3))
print("ŷ1: {:.2f} \t\t\tŷ2: {:.2f} \t\t\tŷ3: {:.2f} \t\tŷ4: {:.2f}".format(yj1, yj2, yj3, yj4))
print("\n")
print("Перевірка за Фішером")
print("Sad²: {:.2f} \nFp: {:.2f}".format(Sad, Fp))
print("\n")

if Fp < Ft:
    print("Pівняння регресії адекватно оригіналу при рівні значимості 0.05")
    cont = False
else:
    cont = True
    print("Pівняння регресії неадекватно оригіналу при рівні значимості 0.05")
    l = 1.215

    x_min = [-5, -1, -4]
    x_max = [8, 4, 2]

    x_0 = [(x_min[0] + x_max[0]) / 2,
           (x_min[1] + x_max[1]) / 2,
           (x_min[2] + x_max[2]) / 2]

    x_l = [l * (x_max[0] - x_0[0]) + x_0[0],
           l * (x_max[1] - x_0[1]) + x_0[1],
           l * (x_max[2] - x_0[2]) + x_0[2]]

    x_cp_min = sum(x_min) / 3
    x_cp_max = sum(x_max) / 3

    ymin = round(200 + x_cp_min)
    ymax = round(200 + x_cp_max)

    xnorm = [[-1, -1, -1],
             [-1, 1, 1],
             [1, -1, 1],
             [1, 1, -1],
             [-1, -1, 1],
             [-1, 1, -1],
             [1, -1, -1],
             [1, 1, 1],
             [-l, 0, 0],
             [l, 0, 0],
             [0, -l, 0],
             [0, l, 0],
             [0, 0, -l],
             [0, 0, l],
             [0, 0, 0]]

    x_nat = [[x_min[0], x_min[1], x_min[2]],
             [x_min[0], x_min[1], x_max[2]],
             [x_min[0], x_max[1], x_min[2]],
             [x_min[0], x_max[1], x_max[2]],
             [x_max[0], x_min[1], x_min[2]],
             [x_max[0], x_min[1], x_max[2]],
             [x_max[0], x_max[1], x_min[2]],
             [x_max[0], x_max[1], x_max[2]],
             [-x_l[0], x_0[1], x_0[2]],
             [x_l[0], x_0[1], x_0[2]],
             [x_0[0], -x_l[1], x_0[2]],
             [x_0[0], x_l[1], x_0[2]],
             [x_0[0], x_0[1], -x_l[2]],
             [x_0[0], x_0[1], x_l[2]],
             [x_0[0], x_0[1], x_0[2]]]

    n = 15
    m = 3
    p = 0.95


    def matrix():
        columns = ["x1", "x2", "x3", "x1*x2", "x1*x3", "x2*x3", "x1*x2*x3",
                   "x1²", "x2²", "x3²", "Yi1", "Yi2", "Yi3", "Ys", "Ye"]

        for i in range(len(columns)):
            print(" {:^9} ".format(columns[i]), end="")
        print()

        for i in range(n):
            for j in range(1, 11):
                print(" {:^9} ".format(round(cmb(x_nat[i])[j], 2)), end="")

            for j in y[i][:-1]:
                print(" {:^9}  ".format(j), end="")
            print(" {:6.2f}   {:6.2f}  "
                  .format(y[i][-1],
                          sum([cmb(x_nat[i])[j] * b[j] * stud[j] for j in range(10)])), end="")
            print()

        print("Функція відгуку:\nY = ", end="")
        if stud[0] != 0:
            print("{:.3f}".format(b[0]), end="")
        for i in range(1, 10):
            if stud[i] != 0:
                print(" + {:.3f}*{}".format(b[i], columns[i]), end="")
        print()


    def geny(n, m, y_max, y_min):
        mat_y = [[randint(y_min, y_max) for j in range(m)] for i in range(n)]
        for elem in mat_y:
            elem.append(sum(elem) / len(elem))
        return mat_y


    def kohren(mat_y, m, n):
        s = []
        for i in range(n):
            ks = 0
            for j in range(m):
                ks += (mat_y[i][-1] - mat_y[i][j]) ** 2
            s.append(ks / m)
        gp = max(s) / sum(s)
        fisher = table_fisher(0.95, n, m, 1)
        gt = fisher / (fisher + (m - 1) - 2)
        return gp < gt


    def cmb(arr):
        return [1, *arr,
                arr[0] * arr[1],
                arr[0] * arr[2],
                arr[1] * arr[2],
                arr[0] * arr[1] * arr[2],
                arr[0] * arr[0],
                arr[1] * arr[1],
                arr[2] * arr[2]]


    def get_b(lmaty):
        a00 = [[],
               [x_nat_mod[0]], [x_nat_mod[1]], [x_nat_mod[2]],
               [x_nat_mod[0], x_nat_mod[1]],
               [x_nat_mod[0], x_nat_mod[2]],
               [x_nat_mod[1], x_nat_mod[2]],
               [x_nat_mod[0], x_nat_mod[1], x_nat_mod[2]],
               [x_nat_mod[0], x_nat_mod[0]],
               [x_nat_mod[1], x_nat_mod[1]],
               [x_nat_mod[2], x_nat_mod[2]]]

        def calcxi(n, listx):
            sumxi = 0
            for i in range(n):
                lsumxi = 1
                for j in range(len(listx)):
                    lsumxi *= listx[j][i]
                sumxi += lsumxi
            return sumxi

        a0 = [15]
        for i in range(10):
            a0.append(calcxi(n, a00[i + 1]))

        a1 = [calcxi(n, a00[i] + a00[1]) for i in range(len(a00))]
        a2 = [calcxi(n, a00[i] + a00[2]) for i in range(len(a00))]
        a3 = [calcxi(n, a00[i] + a00[3]) for i in range(len(a00))]
        a4 = [calcxi(n, a00[i] + a00[4]) for i in range(len(a00))]
        a5 = [calcxi(n, a00[i] + a00[5]) for i in range(len(a00))]
        a6 = [calcxi(n, a00[i] + a00[6]) for i in range(len(a00))]
        a7 = [calcxi(n, a00[i] + a00[7]) for i in range(len(a00))]
        a8 = [calcxi(n, a00[i] + a00[8]) for i in range(len(a00))]
        a9 = [calcxi(n, a00[i] + a00[9]) for i in range(len(a00))]
        a10 = [calcxi(n, a00[i] + a00[10]) for i in range(len(a00))]

        a = numpy.array([[a0[0], a0[1], a0[2], a0[3], a0[4], a0[5],
                          a0[6], a0[7], a0[8], a0[9], a0[10]],
                         [a1[0], a1[1], a1[2], a1[3], a1[4], a1[5],
                          a1[6], a1[7], a1[8], a1[9], a1[10]],
                         [a2[0], a2[1], a2[2], a2[3], a2[4], a2[5],
                          a2[6], a2[7], a2[8], a2[9], a2[10]],
                         [a3[0], a3[1], a3[2], a3[3], a3[4], a3[5],
                          a3[6], a3[7], a3[8], a3[9], a3[10]],
                         [a4[0], a4[1], a4[2], a4[3], a4[4], a4[5],
                          a4[6], a4[7], a4[8], a4[9], a4[10]],
                         [a5[0], a5[1], a5[2], a5[3], a5[4], a5[5],
                          a5[6], a5[7], a5[8], a5[9], a5[10]],
                         [a6[0], a6[1], a6[2], a6[3], a6[4], a6[5],
                          a6[6], a6[7], a6[8], a6[9], a6[10]],
                         [a7[0], a7[1], a7[2], a7[3], a7[4], a7[5],
                          a7[6], a7[7], a7[8], a7[9], a7[10]],
                         [a8[0], a8[1], a8[2], a8[3], a8[4], a8[5],
                          a8[6], a8[7], a8[8], a8[9], a8[10]],
                         [a9[0], a9[1], a9[2], a9[3], a9[4], a9[5],
                          a9[6], a9[7], a9[8], a9[9], a9[10]],
                         [a10[0], a10[1], a10[2], a10[3], a10[4], a10[5],
                          a10[6], a10[7], a10[8], a10[9], a10[10]]])
        c0 = [calcxi(n, [lmaty])]
        for i in range(len(a00) - 1):
            c0.append(calcxi(n, a00[i + 1] + [lmaty]))
        c = numpy.array([c0[0], c0[1], c0[2], c0[3], c0[4], c0[5],
                         c0[6], c0[7], c0[8], c0[9], c0[10]])
        b = numpy.linalg.solve(a, c)

        return b


    def table_student(prob, n, m):
        x_vec = [i * 0.0001 for i in range(int(5 / 0.0001))]
        par = 0.5 + prob / 0.1 * 0.05
        f3 = (m - 1) * n
        for i in x_vec:
            if abs(t.cdf(i, f3) - par) < 0.000005:
                return i


    def table_fisher(prob, n, m, d):
        x_vec = [i * 0.001 for i in range(int(10 / 0.001))]
        f3 = (m - 1) * n
        for i in x_vec:
            if abs(f.cdf(i, n - d, f3) - prob) < 0.0001:
                return i


    def student(n, m, mat_y):
        disp = []
        for i in mat_y:
            s = 0
            for k in range(m):
                s += (i[-1] - i[k]) ** 2
            disp.append(s / m)

        sbt = (sum(disp) / n / n / m) ** (0.5)

        bs = []
        for i in range(11):
            ar = []
            for j in range(len(mat_y)):
                ar.append(mat_y[j][-1] * cmb(xnorm[j])[i] / n)
            bs.append(sum(ar))

        t = [(bs[i] / sbt) for i in range(11)]
        tt = table_student(0.95, n, m)
        st = [i > tt for i in t]
        return st


    def kohren(mat_y, m, n):
        s = []
        for i in range(n):
            ks = 0
            for j in range(m):
                ks += (mat_y[i][-1] - mat_y[i][j]) ** 2
            s.append(ks / m)
        gp = max(s) / sum(s)
        fisher = table_fisher(0.95, n, m, 1)
        gt = fisher / (fisher + (m - 1) - 2)
        return gp < gt


    def fisher(b_0, x_mod, n, m, d, mat_y):
        if d == n:
            return True
        disp = []
        for i in mat_y:
            s = 0
            for k in range(m):
                s += (i[-1] - i[k]) ** 2
            disp.append(s / m)

        sad = sum([(sum([cmb(x_nat[i])[j] * b_0[j] for j in range(11)]) - mat_y[i][-1]) ** 2 for i in range(n)])
        sad = sad * m / (n - d)
        fp = sad / sum(disp) / n
        ft = table_fisher(0.95, n, m, d)
        return fp < ft


    while True:
        while True:
            print("m = {0}\nN = {1}\np = {2}\n".format(m, n, p))
            x_nat_mod = [[x_nat[i][j] for i in range(15)] for j in range(3)]
            y = geny(n, m, ymax, ymin)
            matymod = [y[i][-1] for i in range(len(y))]

            kohren_flag = kohren(y, 3, 15)
            print("Перевірка за Кохреном:\nДисперсія {}однорідна\n"
                  .format("" if kohren_flag else "не "))
            if kohren_flag:
                break
            else:
                m += 1

        b = get_b(matymod)

        stud = student(n, m, y)
        d = sum(stud)

        fisher_ = fisher(b, x_nat_mod, n, m, d, y)
        print("Перевірка за Фішером:\nРівняння {}адекватне\n"
              .format("" if fisher_ else "не "))
        matrix()
        if fisher_:
            break
