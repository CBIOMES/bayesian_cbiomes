import scipy.stats as stats
import matplotlib.pyplot as plt

x = np.arange(-4,4.1,0.1)
pdfx = stats.norm.pdf(x)

fig, ax = plt.subplots()
ax.plot(x, pdfx, marker='o', linestyle='none')
ax.set(xlabel='x', ylabel='pdf')

####################################################################

sampx = stats.norm.rvs(size=500, loc=0, scale=1)

# an alternative is
# sampx = numpy.random.normal(size=500, loc=0, scale=1)

fig, ax = plt.subplots()
ax.plot(sampx, marker='o', linestyle='none')

#####################################################################

fig, ax = plt.subplots()
ax.hist(sampx)

#####################################################################

ksd = stats.gaussian_kde(sampx)

fig, ax = plt.subplots()
ax.plot(x, ksd(x))
ax.set(xlabel='x', ylabel='density')

