import matplotlib.pyplot as plt
import sys
import csv
import numpy as np
from scipy.interpolate import interp1d

N = 1000

theory = lambda:None
theory.x = []
theory.y = []
with open("../data/ICING BE SI Data.csv", 'r') as csvf:
    reader = csv.reader(csvf, delimiter=',')
    for x, y in reader:
        theory.x.append(float(x))
        theory.y.append(float(y))
        
# Convert CDF to PDF
#theory.x = np.append(np.array([0]), np.diff(theory.x))
use_interp = False
if not use_interp:
    theory.y = np.append(np.array([0]), np.diff(theory.y) / np.diff(theory.x))
else:
    theory.y = np.diff(theory.y) / np.diff(theory.x)
# Smooth it out? https://stackoverflow.com/questions/5283649/plot-smooth-line-with-pyplot
if use_interp:
    interpx = np.linspace(0, 0.9, num=N, endpoint=True)
    f_cubic = interp1d(theory.x[:-1], theory.y, kind='linear', fill_value='extrapolate')

data = [float(line) for line in sys.stdin.readlines()]
# Plot the PDF, not the histogram
plt.hist(data, N, density=True)
plt.step(theory.x, theory.y,label="Theoretical", where="pre")
plt.title("Sampled PDF")
plt.show()
