some_list = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]


def once(some_list):
    for elem in some_list:
        print(elem)


once(some_list)


def many_times(some_list):
    for elem in some_list:
        print(elem)
    for elem in some_list:
        print(elem)
    for elem in some_list:
        print(elem)


many_times(some_list)

import numpy as np
def find_the_match(search_list, match_elem):
    for elem in search_list:
        if elem == match_elem:
            return True
    return False


some_list_O2 = np.arange(1, 21)  # [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20]
find_the_match(some_list_O2, 20)


def time_and_space_complexity(n):
    for x in range(n):
        print('something that show only once')


time_and_space_complexity(8)


