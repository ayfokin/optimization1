import math


NUMBER_OF_EPSILON = 10000


def sign(x):
    if x > 0:
        return 1
    if x < 1:
        return -1
    return 0


def f1(x):
    return -5 * x ** 5 + 4 * x ** 4 - 12 * x ** 3 + 11 * x ** 2 - 2 * x + 1


def f2(x):
    return math.log10(x - 2) ** 2 + math.log10(10 - x) ** 2 - x ** 0.2


def f3(x):
    return -3 * x * math.sin(0.75 * x) + math.exp(-2 * x)


def f4(x):
    return math.exp(3 * x) + 5 * math.exp(- 2 * x)


def f5(x):
    return 0.2 * x * math.log10(x) + (x - 2.3) ** 2


def dichotomy(function, epsilon):
    iteration = 0
    delta = epsilon * 0.49
    left_border = function[1]
    right_border = function[2]

    x1 = 0
    while right_border - left_border < epsilon:
        iteration += 1
        middle = (left_border + right_border) / 2
        x1 = middle - delta
        x2 = middle + delta
        fx1 = function[0](x1)
        fx2 = function[0](x2)

        if fx1 > fx2:
            left_border = x1

        if fx1 < fx2:
            right_border = x2

        if fx1 == fx2:
            left_border = x1
            right_border = x2

    return x1


def golden_ratio(function, epsilon):
    iteration = 0
    left_border = function[1]
    right_border = function[2]

    x1 = 0
    while right_border - left_border < epsilon:
        iteration += 1
        x1 = left_border + 0.381966011 * (right_border - left_border)
        x2 = right_border - 0.381966011 * (right_border - left_border)
        fx1 = function[0](x1)
        fx2 = function[0](x2)

        if fx1 > fx2:
            left_border = x1

        if fx1 < fx2:
            right_border = x2

        if fx1 == fx2:
            left_border = x1
            right_border = x2

    return x1


def fibonacci(function, epsilon):
    n = get_n(function[1], function[2], epsilon)
    a = function[1]
    b = function[2]
    x1 = a + get_fibonacci(n) / get_fibonacci(n + 2) * (b - a)
    x2 = a + get_fibonacci(n + 1) / get_fibonacci(n + 2) * (b - a)

    for k in range(n):
        k += 1
        fx1 = function[0](x1)
        fx2 = function[0](x2)

        if fx1 < fx2:
            b = x2
            x2 = x1
            x1 = a + get_fibonacci(n - k + 1) / get_fibonacci(n - k + 3) * (b - a)

        if fx1 > fx2:
            a = x1
            x1 = x2
            x2 = a + get_fibonacci(n - k + 2) / get_fibonacci(n - k + 3) * (b - a)

        if fx1 == fx2:
            return x1
    return x1


def combined_brent(function, epsilon):
    a = function[1]
    c = function[2]

    x = w = v = (a + c) / 2

    fx = fw = fv = function[0](x)

    d = e = c - a

    while (c - a) < epsilon:
        g = e
        e = d

        u = 0
        if x != w and w != v and x != v:
            f2 = function[0](x)
            f1 = function[0](v)
            f3 = function[0](w)
            u = x - 0.5 * ((f2 - f3) * (x - w) ** 2 - (f2 - f1) * (x - v) ** 2) /\
                ((x - w) * (f2 - f3) - (x - v) * (f2 - f1))

        if (a + epsilon < u < c - epsilon) and (abs(u - x) < g / 2):
            d = abs(u - x)

        else:
            if x < (c - a) / 2:
                u = x + 0.381966011 * (c - x)
                d = c - x
            else:
                u = x - 0.381966011 * (x - a)
                d = x - a

        if abs(u - x) < epsilon:
            u = x + sign(u - x) * epsilon

        fu = function[0](u)

        if fu < fx:

            if u >= x:
                a = x
            else:
                c = x

            v = w
            w = x
            x = u
            fv = fw
            fw = fx
            fx = fu

        else:
            if u >= x:
                c = u
            else:
                a = u

            if fu <= fw or w == x:
                v = w
                w = u
                fv = fw
                fw = fu
            elif fu <= fv or v == x or v == w:
                v = u
                fv = fu
    return x
                    


def get_fibonacci(n):
    return round(1 / math.sqrt(5) * (((1 + math.sqrt(5)) / 2) ** n - ((1 - math.sqrt(5)) / 2) ** n))


def get_n(left_border, right_border, epsilon):
    argument = math.sqrt(5) * (((right_border - left_border) / epsilon) - 0.5)
    return math.ceil(math.log(argument, ((1 + math.sqrt(5)) / 2)))


functions = [[f1, -0.5, 0.5], [f2, 6, 9.9], [f3, 0, 2 * math.pi], [f4, 0, 1], [f5, 0.5, 2.5]]
print(get_n(-0.5, 0.5, 0.00001))

#for function in functions:
 #   iteration = 0
  #  for epsilon in range(NUMBER_OF_EPSILON):
   #     epsilon = epsilon * 1 / NUMBER_OF_EPSILON + 1 / NUMBER_OF_EPSILON


def parabolic(function, epsilon):
    x1 = function[1]
    x3 = function[2]
    f1 = function[0](x1)
    f3 = function[0](x3)

    while x3 - x1 > epsilon:
        x2 = (x1 + x3) / 2
        f2 = function[0](x2)
        u = x2 - 0.5*((f2 - f3)*(x2 - x1)**2 - (f2 - f1)*(x2 - x3)**2)/((x2 - x1)*(f2 - f3) - (x2 - x3)*(f2 - f1))
        fu = function[0](u)
        if x2 < u:
            x2, u, f2, fu = u, x2, fu, f2
        if f2 < fu:
            x1 = u
            f1 = fu
        if f2 > fu:
            x3 = x2
            f3 = f2
        if f2 == fu:
            x1 = u
            x3 = x2
            f1 = fu
            f3 = f2

    return (x1 + x3) / 2


