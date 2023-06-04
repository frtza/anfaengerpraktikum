# essential libraries
import matplotlib.pyplot as plt
import numpy as np
# additional libraries
from scipy.optimize import curve_fit


# read columns of data from txt file
f, U_B = np.genfromtxt('data.txt', unpack=True) 

U_S = 1 # volt
R = 664 # ohm
C = 450 * 10 ** -9 # farad
f_0 = (2 * np.pi * R * C) ** -1 # hertz

# calculate theory curve
def theory(t):
	return np.sqrt( (1 / 9) * (((t ** 2 - 1) ** 2) / ((1 - t ** 2) ** 2 + 9 * t ** 2)))

# plot data and theory curve
plt.figure(figsize=(5.78, 5.78)) # default (5.78, 3.57)
plt.subplot(2, 1, 1)
x = np.linspace(0.01, 60, 5999)
plt.plot(x, theory(x), '-', c='olivedrab', label='Theoriekurve')
plt.plot(f/f_0, U_B/U_S, 'kx', markersize=3, label='Messdaten')
plt.xscale('linear')
plt.xlim(plt.xlim()[0], 31.5)
plt.xlabel(r'$\symup{\Omega} = \pfrac{\nu}{\raisebox{1ex}{\(\nu_0\)}}$')
plt.ylabel(r'$\symup{\pfrac{U_{\! Br}}{U_{\! S}}}$')
leg = plt.legend(loc='best')
leg.get_frame().set_edgecolor('k')
leg.get_frame().set_facecolor('none')
leg.get_frame().set_linewidth(0.25)
plt.subplot(2, 1, 2)
x = np.linspace(0, 80, 8000) 
plt.plot(x, theory(x), '-', c='olivedrab', label='Theoriekurve')
plt.plot(f/f_0, U_B/U_S, 'kx', markersize=3, label='Messdaten')
plt.xscale('log')
plt.xlabel(r'$\symup{\Omega} = \pfrac{\nu}{\raisebox{1ex}{\(\nu_0\)}}$')
plt.ylabel(r'$\symup{\pfrac{U_{\! Br}}{U_{\! S}}}$')
leg = plt.legend(loc='best')
leg.get_frame().set_edgecolor('k')
leg.get_frame().set_facecolor('none')
leg.get_frame().set_linewidth(0.25)
plt.savefig('build/plot.pdf')
plt.close()

# calculate curve fit
def fit(t, k):
	return (1 / 3) * (k + 1) - theory(t) - k * theory(2 * t)
# evaluate fit parameters
par, cov = curve_fit(fit, f/f_0, U_B/U_S)
k = par[0]

# fit curve
plt.figure(figsize=(5.78, 2.89))
plt.xscale('log')
x = np.linspace(0, 65, 10000)
plt.plot(x, theory(x) + k * theory(2 * x), '--', dashes=(3.1, 1.1), c='sandybrown', label=r'$F_0$')
plt.plot(x, fit(x, k), '-', c='sandybrown', label=r'$F_1$')
plt.plot(f/f_0, U_B/U_S, 'kx', markersize=3, label='Messdaten')
plt.xlabel(r'$\symup{\Omega} = \pfrac{\nu}{\raisebox{1ex}{\(\nu_0\)}}$')
plt.ylabel(r'$\symup{\pfrac{U_{\! Br}}{U_{\! S}}}$')
leg = plt.legend(loc='best', handlelength=2.075)
leg.get_frame().set_edgecolor('k')
leg.get_frame().set_facecolor('none')
leg.get_frame().set_linewidth(0.25)
plt.savefig('build/fit.pdf')

# save parameter
with open('build/k.tex', 'w') as file:
	file.write(r'\num{')
	file.write(f'{k:.2f}')
	file.write(r'}')

# format latex table
table_header = r''' \begin{tabular}
		{S[table-format=5.0]
		 S[table-format=1.3]
		 S[table-format=1.2]
		 S[table-format=1.3]
		 S[table-format=1.3]}
		\toprule
		{$\nu \mathbin{/} \unit{\hertz}$} &
		{$\symup{U_{\! Br}} \mathbin{/} \unit{\volt}$} &
		{$\sfrac{\nu}{\nu_0}$} &
		{$\symup{\sfrac{U_{\! Br}}{\, U_{\! S}}}$} &
		{Theorie} \\
		\midrule
'''
table_footer = r''' 	\bottomrule
	\end{tabular}
'''
row_template = r'		{0:5.0f} & {1:1.3f} & {2:1.2f} & {3:1.3f} & {4:1.3f} \\'
# export table to build directory
with open('build/table_e.tex', 'w') as file:
	file.write(table_header)
	for row in zip(f, U_B, f/f_0, U_B/U_S, theory(f/f_0)):
		file.write(row_template.format(*row))
		file.write('\n')
	file.write(table_footer)
