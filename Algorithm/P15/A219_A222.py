def BigO(n: int) -> int:
    return 88 * n ** 3 + 99 * n ** 2 + 100


stack = []
for i in range(10):
    stack.append(BigO(i + 1))

for n in range(10):
    x = 88 * n ** 3 - 99 * n ** 2
    while (88 * n ** 3) <= (99 * n ** 2):
        n += 1
        y = 88 * n ** 3 - 99 * n ** 2

from math import log
import numpy as np
import matplotlib.pyplot as plt

plt.style.use('bmh')
n = np.linspace(1, 10)
labels = ['Constant', 'Logrithmic', 'Linear', 'Loglienar', 'Quadratic', 'Cubic', 'exponential']

big_o = [np.ones(n.shape), np.log(n), n, n * np.log(n), n ** 2, n ** 3, 2 ** n]

plt.figure(figsize=(10, 8))
plt.ylim(0, 50)

for i in range(len(big_o)):
    plt.plot(n, big_o[i], label=labels[i])

plt.legend(loc=0)
plt.ylabel('(Relative) Runtime fpr Comparison')
plt.xlabel('n')
plt.show()
