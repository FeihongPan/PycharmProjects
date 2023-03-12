import sys

n = 10
test_data = []
for i in range(n):
    x = len(test_data)
    y = sys.getsizeof(test_data)

    print('the length of the data list is :{}: The size in memory(byte) is :{}'.format(x, y))
    test_data.append(n)

from IPython.display import Image
Image(filename='A239.png')
