import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import AutoMinorLocator
from scipy.optimize import curve_fit

def format_func(value, tick_number):
    N = int(np.round(2 * value / np.pi))
    if N == 0:
        return "0"
    elif N == 1:
        return r"$\pi \mkern-0.5\thinmuskip /2$"
    elif N == 2:
        return r"$\vphantom{0}\pi$"
    elif N % 2 > 0:
        return r"${0}\pi \mkern-0.5\thinmuskip /2$".format(N)
    else:
        return r"${0}\pi$".format(N // 2)

phi, U_a, U_b = np.genfromtxt('data_1.txt', unpack=True) 
r, U = np.genfromtxt('data_2.txt', unpack=True) 

def fit_1(t, a, b, c, d):
	return a * np.cos(b * t + c) + d

par_a, cov_a = curve_fit(fit_1, phi, U_a)
par_b, cov_b = curve_fit(fit_1, phi, U_b)

x = np.linspace(- np.pi / 10, 16 * np.pi / 10, 10000) 

ax = plt.subplot(1,1,1)
plt.plot(x, fit_1(x, *par_a), c='olivedrab', label='Fit-Kurve')
plt.plot(np.pi * phi / 180, U_a, 'kx', markersize=3.21, label='Messdaten')
plt.xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2],
		   [r'$\vphantom{\pfrac{0}{0}}0$',
			r'$\pfrac{1}{\raisebox{1ex}{\(2\)}}\pi$',
			r'$\vphantom{\pfrac{0}{0}}\pi$',
			r'$\pfrac{3}{\raisebox{1ex}{\(2\)}}\pi$'])
plt.xlabel(r'$\phi \mathbin{/} \unit{\radian}$')
plt.ylabel(r'$U \mathbin{/} \unit{\milli\volt} $')
leg = plt.legend(loc='best', edgecolor='k', facecolor='none')
leg.get_frame().set_linewidth(0.25)
ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
plt.savefig('build/plot_1a.pdf')
plt.close()

ax = plt.subplot(1,1,1)
plt.plot(x, fit_1(x, *par_b), c='olivedrab', label='Fit-Kurve')
plt.plot(np.pi * phi / 180, U_b, 'kx', markersize=3.21, label='Messdaten')
plt.xticks([0, np.pi / 2, np.pi, 3 * np.pi / 2],
		   [r'$\vphantom{\pfrac{0}{0}}0$',
			r'$\pfrac{1}{\raisebox{1ex}{\(2\)}}\pi$',
			r'$\vphantom{\pfrac{0}{0}}\pi$',
			r'$\pfrac{3}{\raisebox{1ex}{\(2\)}}\pi$'])
plt.xlabel(r'$\phi \mathbin{/} \unit{\radian}$')
plt.ylabel(r'$U \mathbin{/} \unit{\milli\volt} $')
leg = plt.legend(loc='best', edgecolor='k', facecolor='none')
leg.get_frame().set_linewidth(0.25)
ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
plt.savefig('build/plot_1b.pdf')
plt.close()

def fit_2(t, a, b):
	return a / t ** 2 + b

par_1, cov_1 = curve_fit(fit_2, r[0:len(r) - 4], U[0:len(r) - 4] / 1000)
par_2, cov_2 = curve_fit(fit_2, r, U / 1000)

ax = plt.subplot(1,1,1)
x = np.linspace(9.9, 150, 10000)
plt.plot(x, fit_2(x, *par_2), c='olivedrab', alpha=0.5, label='Regulärer Fit')
x = np.linspace(12.18, 150, 10000)
plt.plot(x, fit_2(x, *par_1), c='olivedrab', label='Angepasster Fit')
plt.plot(r, U / 1000, 'kx', markersize=3.21, label='Messdaten')
plt.xlabel(r'$r \mathbin{/} \unit{\centi\meter}$')
plt.ylabel(r'$U \mathbin{/} \unit{\volt} $')
leg = plt.legend(loc='best', edgecolor='k', facecolor='none')
leg.get_frame().set_linewidth(0.25)
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
plt.savefig('build/plot_2.pdf')

table_header_1 = r''' \begin{tabular}
		{S[table-format=3.0]
		 S[table-format=1.2]
		 S[table-format=3.1]}
		\toprule
		{$\phi \mathbin{/} \unit{\degree}$} & 
		{$\phi \mathbin{/} \unit{\radian}$} & 
		{$U \mathbin{/} \unit{\milli\volt}$} \\
		\midrule
'''
table_header_2 = r''' \begin{tabular}
		{S[table-format=2.0]
		 @{\hspace{2ex}}
		 S[table-format=4.0]
		 @{\hspace{4ex}}
		 S[table-format=3.0]
		 @{\hspace{2ex}}
		 S[table-format=2.1]}
		\toprule
		{$r \mathbin{/} \unit{\centi\meter}$} & 
		{$U \mathbin{/} \unit{\milli\volt}$}  &
		{$r \mathbin{/} \unit{\centi\meter}$} & 
		{$U \mathbin{/} \unit{\milli\volt}$}  \\
		\midrule
'''
table_footer = r''' 	\bottomrule
	\end{tabular}
'''
row_template_1 = r'		{0:3.0f} & {1:1.2f} & {2:3.1f} \\'
row_template_2 = r'		{0:2.0f} & {1:4.0f} & {2:3.0f} & {3:2.1f} \\'

rad = [r'$0$',
	   r'$\pi \mkern-0.5\thinmuskip /6$',
	   r'$\pi \mkern-0.5\thinmuskip /3$',
	   r'$\pi \mkern-0.5\thinmuskip /2$',
	   r'$2\pi \mkern-0.5\thinmuskip /3$',
	   r'$5\pi \mkern-0.5\thinmuskip /6$',
	   r'$\pi$',
	   r'$7\pi \mkern-0.5\thinmuskip /6$',
	   r'$4\pi \mkern-0.5\thinmuskip /3$',
	   r'$3\pi \mkern-0.5\thinmuskip /2$']

with open('build/table_1a.tex', 'w') as file:
	file.write(table_header_1)
	for row in zip(phi, np.pi * phi / 180, U_a):
		file.write(row_template_1.format(*row))
		file.write('\n')
	file.write(table_footer)

with open('build/table_1b.tex', 'w') as file:
	file.write(table_header_1)
	for row in zip(phi, np.pi * phi / 180, U_b):
		file.write(row_template_1.format(*row))
		file.write('\n')
	file.write(table_footer)

with open('build/table_2.tex', 'w') as file:
	file.write(table_header_2)
	for row in zip(r[11:23][::-1], U[11:23][::-1], r[0:11][::-1], U[0:11][::-1]):
		file.write(row_template_2.format(*row))
		file.write('\n')
	file.write(f' 	    {r[11]:2.0f} & {U[11]:4.0f} & & ')
	file.write(r'\\')
	file.write('\n')
	file.write(table_footer)

with open('build/a1.tex', 'w') as file:
	file.write(r'\qty{')
	file.write(f'{par_a[0]:.2f}')
	file.write(r'(')
	file.write(f'{np.sqrt(np.diag(cov_a))[0]:.2f}')
	file.write(r')')
	file.write(r'}{\milli\volt}')

with open('build/b1.tex', 'w') as file:
	file.write(r'\num{')
	file.write(f'{par_a[1]:.2f}')
	file.write(r'(')
	file.write(f'{np.sqrt(np.diag(cov_a))[1]:.2f}')
	file.write(r')')
	file.write(r'}')

with open('build/c1.tex', 'w') as file:
	file.write(r'\num{')
	file.write(f'{par_a[2]:.2f}')
	file.write(r'(')
	file.write(f'{np.sqrt(np.diag(cov_a))[2]:.2f}')
	file.write(r')')
	file.write(r'}')

with open('build/d1.tex', 'w') as file:
	file.write(r'\qty{')
	file.write(f'{par_a[3]:.2f}')
	file.write(r'(')
	file.write(f'{np.sqrt(np.diag(cov_a))[3]:.2f}')
	file.write(r')')
	file.write(r'}{\milli\volt}')

with open('build/a2.tex', 'w') as file:
	file.write(r'\qty{')
	file.write(f'{par_b[0]:.2f}')
	file.write(r'(')
	file.write(f'{np.sqrt(np.diag(cov_b))[0]:.2f}')
	file.write(r')')
	file.write(r'}{\milli\volt}')

with open('build/b2.tex', 'w') as file:
	file.write(r'\num{')
	file.write(f'{par_b[1]:.2f}')
	file.write(r'(')
	file.write(f'{np.sqrt(np.diag(cov_b))[1]:.2f}')
	file.write(r')')
	file.write(r'}')

with open('build/c2.tex', 'w') as file:
	file.write(r'\num{')
	file.write(f'{par_b[2]:.2f}')
	file.write(r'(')
	file.write(f'{np.sqrt(np.diag(cov_b))[2]:.2f}')
	file.write(r')')
	file.write(r'}')

with open('build/d2.tex', 'w') as file:
	file.write(r'\qty{')
	file.write(f'{par_b[3]:.2f}')
	file.write(r'(')
	file.write(f'{np.sqrt(np.diag(cov_b))[3]:.2f}')
	file.write(r')')
	file.write(r'}{\milli\volt}')

with open('build/r1.tex', 'w') as file:
	file.write(r'\qty{')
	file.write(f'{par_1[0]:.2f}')
	file.write(r'(')
	file.write(f'{np.sqrt(np.diag(cov_1))[0]:.2f}')
	file.write(r')')
	file.write(r'}{\milli\volt\centi\meter\squared}')

with open('build/s1.tex', 'w') as file:
	file.write(r'\qty{')
	file.write(f'{par_1[1]:.2f}')
	file.write(r'(')
	file.write(f'{np.sqrt(np.diag(cov_1))[1]:.2f}')
	file.write(r')')
	file.write(r'}{\milli\volt}')

with open('build/r2.tex', 'w') as file:
	file.write(r'\qty{')
	file.write(f'{par_2[0]:.2f}')
	file.write(r'(')
	file.write(f'{np.sqrt(np.diag(cov_2))[0]:.2f}')
	file.write(r')')
	file.write(r'}{\milli\volt\centi\meter\squared}')

with open('build/s2.tex', 'w') as file:
	file.write(r'\qty{')
	file.write(f'{par_2[1]:.2f}')
	file.write(r'(')
	file.write(f'{np.sqrt(np.diag(cov_2))[1]:.2f}')
	file.write(r')')
	file.write(r'}{\milli\volt}')
