
list_minima = ((y[1:-1] < y[:-2]) & (y[1:-1] <= y[2:])).nonzero()[0]
list_maxima = ((y[1:-1] > y[:-2]) & (y[1:-1] >= y[2:])).nonzero()[0]

values_min = y[list_minima[1:] + 1]  # Exclude first minima, as it is the only one that is above zero
values_max = y[list_maxima + 1]

hist_min = np.histogram(values_min, bins = 500)
hist_max = np.histogram(values_max, bins = 500)

plt.figure()
plt.plot(hist_min[1][1:], hist_min[0])
plt.plot(hist_max[1][1:], hist_max[0])