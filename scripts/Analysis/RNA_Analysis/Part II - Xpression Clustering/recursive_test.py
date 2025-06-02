import numpy as np
import matplotlib.pyplot as plt

np.random.seed(33)

# Between -1 and 1 (normal)
# def recurs_test(loc, size):
#     norm_correlation_values = list(np.random.normal(loc = loc, scale = 0.1, size = size))
#     norm_naked_sum = 0
#     for r in norm_correlation_values:
#         w = 2 if r > 0 else 1
#         norm_naked_sum += w * (r**2)
#     result = norm_naked_sum / len(norm_correlation_values)
#     return result


# print(recurs_test(0.2, 1000))
# print(recurs_test(-0.2, 1000 ))


# fig1 = plt.figure()
# plt.xlabel('Size')
# plt.ylabel('Coef Score')
# xvals  = []
# yvals = []
# for i in range(10,1000,2): 
#     print(recurs_test(i))
#     yvals.append(recurs_test(i))
#     xvals.append(i)
# plt.plot(xvals, yvals)
# plt.show()





norm_correlation_values = list(np.random.normal(loc = 0, scale = 0.1, size = 10))
# norm_naked_sum = 0
# for r in norm_correlation_values:
# w = 2 if r > 0 else 1
# norm_naked_sum += w * (r**2)
# result = norm_naked_sum / len(norm_correlation_values)
# return result



print(norm_correlation_values)