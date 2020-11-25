import math


NUMBER_OF_EPSILON = 10000


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
    while True:
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
        if right_border - left_border < epsilon:
            break


def golden_ratio(function, epsilon):
    iteration = 0
    left_border = function[1]
    right_border = function[2]
    while True:
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
        if right_border - left_border < epsilon:
            break


def fibonacci(function, epsilon):
    n = get_n(function[1], function[2], epsilon)
    x1 = function[1]
    for i in range(n):

    return 1


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