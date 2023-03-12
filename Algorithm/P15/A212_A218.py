# ctrl + alt + L = 调试代码位置
# shift + F9 -> F8 Debug模式
import timeit
def some_func_to_sum(n):
    final_result = 0
    for x in range(n + 1):
        final_result += x
    return final_result


def another_func_to_sum(n):
    final_result = 0.5 * n * (1 + n)
    return final_result


another_func_to_sum(5)
some_func_to_sum(5)
