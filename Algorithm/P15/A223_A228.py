n = 2
while 2 ** n < n ** 3:
    n += 1
    x = 2 ** n - n ** 3

some_list = [1,2,3,4,5]


def func_constant(values):
    print(values[0])


func_constant(some_list)


def func_linear(some_list):
    for ch in some_list:
        print(ch)


func_linear(some_list)


def func_nested_for_loop(some_list):
    for i in some_list:
        for j in some_list:
            print(i, j)


func_nested_for_loop(some_list)
