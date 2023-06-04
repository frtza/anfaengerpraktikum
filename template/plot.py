# run shell commands
import subprocess
# essential libraries
import matplotlib.pyplot as plt
import numpy as np
# additional libraries
from scipy.optimize import curve_fit
import uncertainties.unumpy as unp # explicit correlation analysis not included

# read columns of data from txt file
x_data, y_data = np.genfromtxt('data.txt', unpack=True) 

# define expected curve to fit data to with parameters a and b
def f(t, a, b):
	return a * np.exp(b * t)

# parameter array and covariance matrix for linear regression or polynomials
par, cov = np.polyfit(x_data, y_data, deg=1, cov=True) 
# parameter array and covariance matrix for least squares curve fit
par, cov = curve_fit(f, x_data, y_data)

# calculate average of measurements
avg = np.average(y_data/x_data)

# steps for x axis
x = np.linspace(0, 10, 10000) 
# calculate values for y axis
y = x ** 2 

# overwrite matplotlibrc figure dimensions with width and height in inches
plt.figure(figsize=(2, 1))
# multiple plots next to each other with vertical and horizontal number
plt.subplot(1, 1, 1)
plt.plot(x, y, label='')
# labels in latex code
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
# axis limits
plt.xlim(0, 10)
plt.ylim(0, 10)
# display legend
plt.legend(loc='best')
# save figure as vector graphic
plt.savefig('build/plot.pdf')

# format latex table
table_header = r''' \begin{tabular}
		{S[table-format=1.2]
		 S[table-format=1.2]}
		\toprule
		{$x \mathbin{/} \unit{\second}$} & 
		{$y \mathbin{/} \unit{\meter}$} \\
		\midrule
'''
table_footer = r''' 	\bottomrule
	\end{tabular}
'''
row_template = r'		{0:1.2f} & {1:1.2f} \\'
# export table to build directory
with open('build/table.tex', 'w') as f:
	f.write(table_header)
	for row in zip(x_data, y_data):
		f.write(row_template.format(*row))
		f.write('\n')
	f.write(table_footer)

# format and export calculated values to build directory
with open('build/avg.tex', 'w') as f: 
	f.write(r'\num{')
	f.write(f'{avg:.2f}')
	f.write(r'}')

# automatically calculate error propagation, format and save as latex
with open('build/gauss_err.tex', 'w') as f:
	process = subprocess.Popen(['../template/extension/python/gauss_err.py', 'x**2 + y**2'], stdout=f)
