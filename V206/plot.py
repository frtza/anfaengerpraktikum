# essential libraries
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
# additional libraries
import scipy.constants as const
from matplotlib.ticker import AutoMinorLocator
from scipy.optimize import curve_fit

# update latex settings
pgf_latex = {'pgf.preamble': '\n'.join([r'\usepackage[per-mode=reciprocal]{siunitx}'])}
mpl.rcParams.update(pgf_latex)

# read columns of data from txt file
T_2, p_a, p_b, T_1, P, t = np.genfromtxt('data.txt', unpack=True) 

# define polynomial temperature curve
def T(t, a, b, c):
	return a*t**2+b*t+c
par_1, cov_1 = np.polyfit(t*60, const.convert_temperature(T_1, 'C', 'K'), deg=2, cov=True)
err_1 = np.sqrt(np.diag(cov_1))
par_2, cov_2 = np.polyfit(t*60, const.convert_temperature(T_2, 'C', 'K'), deg=2, cov=True)
err_2 = np.sqrt(np.diag(cov_2))

# steps for x axis
x = np.linspace(-0.475*60, 24.275*60, 10000) 

# define axis object
ax = plt.subplot(1,1,1)
# set axis tick markers
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
# plot fitted curve
f_1, = plt.plot(x, T(x, *par_1), '-', color='olivedrab', alpha=0.4, label='Fit an $T_1$')
f_2, = plt.plot(x, T(x, *par_2), '-', color='steelblue', alpha=0.4, label='Fit an $T_2$')
# plot measurement data
p_1, = plt.plot(t*60, const.convert_temperature(T_1, 'C', 'K'), 'x', markersize=3.21,
				color='olivedrab', label='Daten für $T_1$')
p_2, = plt.plot(t*60, const.convert_temperature(T_2, 'C', 'K'), 'x', markersize=3.21,
				color='steelblue', label='Daten für $T_1$')
# placeholder plot for group labels
d, = plt.plot(200, 300, marker='None', linestyle='None', label='dummy 1')
# labels in latex code
plt.xlabel(r'$t \mathbin{/} \unit{\second}$')
plt.ylabel(r'$T \mathbin{/} \unit{\kelvin}$')
# display legend
leg = plt.legend([d, p_1, f_1, d, p_2, f_2],
				 ['$T_1$', 'Messdaten', 'Regression', '$T_2$', 'Messdaten', 'Regression'],
				 loc='best', ncol=2, edgecolor='k', facecolor='none')
leg.get_frame().set_linewidth(0.25)
# save figure as vector graphic
plt.savefig('build/plot_temperature.pdf')
# clear plot
plt.close()

# define exponential steam pressure curve
def p(temp, p_0, q_0):
	return p_0*np.exp(-q_0*temp)
par_3, cov_3 = curve_fit(p, 1/const.convert_temperature(T_1[:-1], 'C', 'K'), (p_b[:-1]+1)*100)
err_3 = np.sqrt(np.diag(cov_3))

# steps for x axis
x = np.linspace(3.125/1000, 3.4125/1000, 10000) 

# define axis object
ax = plt.subplot(1,1,1)
# set axis tick markers
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
# plot fitted curve
f_3, = plt.plot(x, p(x, *par_3), '-', color='salmon', alpha=0.4, label='Regression')
# plot measurement data
p_3, = plt.plot(1/const.convert_temperature(T_1[:-1], 'C', 'K'), (p_b[:-1]+1)*100, 'x',
				markersize=3.21, color='salmon', label='Messdaten')
# labels in latex code
plt.xlabel(r'$T^{\, -1} \! \mathbin{/} \unit{\per\kelvin}$')
plt.ylabel(r'$p \mathbin{/} \unit{\kilo\pascal}$')
# display legend
leg = plt.legend([p_3, f_3], ['Messdaten', 'Regression'],
				 loc='best', edgecolor='k', facecolor='none')
leg.get_frame().set_linewidth(0.25)
# save figure as vector graphic
plt.savefig('build/plot_pressure.pdf')
# clear plot
plt.close()

# format latex table
table_header = r''' 	\begin{tabular}
		{S[table-format=2.1]
		 S[table-format=2.1]
		 S[table-format=3.1]
		 S[table-format=2.1]
		 S[table-format=3.1]
		 S[table-format=1.1]
		 S[table-format=3.0]
		 S[table-format=2.1]
		 S[table-format=4.0]
		 S[table-format=3.0]}
		\toprule
		{$t \mathbin{/} \unit{\minute}$} & 
		{$T_1 \mathbin{/} \unit{\celsius}$} & 
		{$T_1 \mathbin{/} \unit{\kelvin}$} & 
		{$T_2 \mathbin{/} \unit{\celsius}$} & 
		{$T_2 \mathbin{/} \unit{\kelvin}$} & 
		{$\tilde{p}_a \mathbin{/} \unit{\bar}$} & 
		{$p_a \mathbin{/} \unit{\kilo\pascal}$} & 
		{$\tilde{p}_b \mathbin{/} \unit{\bar}$} & 
		{$p_b \mathbin{/} \unit{\kilo\pascal}$} & 
		{$P \mathbin{/} \unit{\watt}$} \\
		\midrule
'''
table_footer = r''' 		\bottomrule
	\end{tabular}
'''
row_template = r'		{0:2.1f} & {1:2.1f} & {2:3.1f} & {3:2.1f} & {4:3.1f} & {5:1.1f} & {6:3.0f} & {7:2.1f} & {8:4.0f} & {9:3.0f} \\'

# export table to build directory
with open('build/table.tex', 'w') as f:
	f.write(table_header)
	for row in zip(t, T_1, const.convert_temperature(T_1, 'C', 'K'), T_2, const.convert_temperature(T_2, 'C', 'K'), p_a, (p_a+1)*100, p_b, (p_b+1)*100, P):
		f.write(row_template.format(*row))
		f.write('\n')
	f.write(table_footer)

# export calculated temperature values
with open('build/A_1.tex', 'w') as f:
	f.write(r'\qty[per-mode=reciprocal]{')
	f.write(f'{par_1[0]*1000000:.2f}({err_1[0]*1000000:.2f})')
	f.write(r'}{\micro\kelvin\per\second\squared}')
with open('build/B_1.tex', 'w') as f:
	f.write(r'\qty[per-mode=reciprocal]{')
	f.write(f'{par_1[1]*1000:.2f}({err_1[1]*1000:.2f})')
	f.write(r'}{\milli\kelvin\per\second}')
with open('build/C_1.tex', 'w') as f:
	f.write(r'\qty[per-mode=reciprocal]{')
	f.write(f'{par_1[2]:.2f}({err_1[2]:.2f})')
	f.write(r'}{\kelvin}')
with open('build/A_2.tex', 'w') as f:
	f.write(r'\qty[per-mode=reciprocal]{')
	f.write(f'{par_2[0]*1000000:.2f}({err_2[0]*1000000:.2f})')
	f.write(r'}{\micro\kelvin\per\second\squared}')
with open('build/B_2.tex', 'w') as f:
	f.write(r'\qty[per-mode=reciprocal]{')
	f.write(f'{par_2[1]*1000:.2f}({err_2[1]*1000:.2f})')
	f.write(r'}{\milli\kelvin\per\second}')
with open('build/C_2.tex', 'w') as f:
	f.write(r'\qty[per-mode=reciprocal]{')
	f.write(f'{par_2[2]:.2f}({err_2[2]:.2f})')
	f.write(r'}{\kelvin}')

# export calculated pressure values
with open('build/R.tex', 'w') as f:
	f.write(r'\qty[per-mode=reciprocal]{')
	f.write(f'{const.R}')
	f.write(r'}{\joule\per\mole\per\kelvin}')
with open('build/L_R.tex', 'w') as f:
	f.write(r'\qty[per-mode=reciprocal]{')
	f.write(f'{par_3[1]:.2f}({err_3[1]:.2f})')
	f.write(r'}{\kelvin}')
with open('build/L.tex', 'w') as f:
	f.write(r'\qty[per-mode=reciprocal]{')
	f.write(f'{par_3[1]*const.R:.2f}({err_3[1]*const.R:.2f})')
	f.write(r'}{\joule\per\mole}')
with open('build/p.tex', 'w') as f:
	f.write(r'\qty[per-mode=reciprocal]{')
	f.write(f'{par_3[0]/1000000:.2f}({err_3[0]/1000000:.2f})')
	f.write(r'}{\giga\pascal}')
