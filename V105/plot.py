import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np
import scipy.constants as const
import uncertainties.unumpy as unp
from uncertainties import ufloat
# Exportiere Konstante
mu_0 = const.physical_constants["vacuum mag. permeability"][0] # Feldkonstante [N/A^2]
with open('build/mu_0.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{mu_0 * 10**6:.3f}e-6')
	f.write(r'}{\newton\per\ampere\squared}')
g = const.physical_constants["standard acceleration of gravity"][0]
with open('build/g.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{g:.3f}')
	f.write(r'}{\meter\per\second\squared}')
# Untersuchung des Magnetfeldverlaufs
N = 195 # Windungszahl
I = 1 # Stromstärke [A]
R = 0.109 # Radius [m]
d = 0.138 # Abstand [m]
r = 0.025 # Kugelradius [m]
M = 0.150 # Kugelmasse [kg]
m = 0.0014 # Testmasse [kg]
# Aufstellen der Feldfunktionen
def B_0(x):
	return mu_0/2 * N * I * R**2 / (R**2 + x**2)**(3/2)
def dB_0(x):
	return - 3/2 * mu_0 * N * I * R**2 * x / (R**2 + x**2)**(5/2)
def B(x):
	return B_0(x + d/2) + B_0(x - d/2)
def dB(x):
	return dB_0(x + d/2) + dB_0(x - d/2)
# Exportierte Werte
B_null = B(0)
with open('build/B_null.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{B_null * 1000:.3f}e-3')
	f.write(r'}{\tesla}')
J = 2/5 * M * r**2
with open('build/J.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{J * 100000:.2f}e-5')
	f.write(r'}{\kilo\gram\meter\squared}')
# Graphische Darstellung
x = np.linspace(-0.25,0.25,10000)
fig, ax1 = plt.subplots()
ax1.xaxis.set_minor_locator(AutoMinorLocator(5))
ax1.yaxis.set_minor_locator(AutoMinorLocator(2))
ax1.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax1.set_xlabel(r'$x \mathbin{/} \unit{\meter}$')
ax1.set_ylabel(r'$B \mathbin{/} \unit{\tesla}$', c='olivedrab')
p1, = ax1.plot(x, B(x), c='olivedrab')
ax2 = ax1.twinx()
ax2.yaxis.set_minor_locator(AutoMinorLocator(4))
ax2.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
ax2.set_ylabel(r'$\pfrac{\symup dB}{\symup dx} \mathbin{/} \unit{\tesla\per\meter}$', c='olivedrab', alpha=0.6, labelpad=-1)
p2, = ax2.plot(x, dB(x), c='olivedrab', alpha=0.5)
leg = ax1.legend([p1, p2], ['Magnetfeld', 'Feldgradient'], edgecolor='k', facecolor='none')
leg.get_frame().set_linewidth(0.25)
ax1.set_zorder(ax2.get_zorder()+1)
ax1.set_frame_on(False)
fig.set_size_inches(5.78, 3.07)
fig.savefig('build/plot_feld.pdf')
plt.close()
# Messdaten einlesen
r_grav, I_grav = np.genfromtxt('data/gravitation.txt', unpack=True)
r_grav = r_grav / 1000
B_grav = B_null * I_grav
I_schw, T_schw = np.genfromtxt('data/schwingung.txt', unpack=True)
B_schw = B_null * I_schw
T_schw = T_schw / 10
I_präz, T_präz = np.genfromtxt('data/präzession.txt', unpack=True)
B_präz = B_null * I_präz
# Graphik zur Gravitation
par, cov = np.polyfit(B_grav, r_grav, deg=1, cov=True)
err = np.sqrt(np.diag(cov))
a_grav = ufloat(par[0],err[0])
with open('build/a_grav.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{a_grav.n:.2f}({a_grav.s:.2f})')
	f.write(r'}{\meter\per\tesla}')
b_grav = ufloat(par[1],err[1])
with open('build/b_grav.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{b_grav.n * 1000:.2f}({b_grav.s * 1000:.2f})e-3')
	f.write(r'}{\meter}')
x = np.linspace(0.00155,0.00345,2)
ax = plt.subplot(111)
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.plot(x, par[0] * x + par[1], '-', c='olivedrab', label='Ausgleichsgerade')
plt.plot(B_grav, r_grav, 'kx', ms=3.21, label='Messergebnisse')
plt.xlabel(r'$B \mathbin{/} \unit{\tesla}$')
plt.ylabel(r'$r \mathbin{/} \unit{\meter}$')
leg = plt.legend(edgecolor='k', facecolor='none')
leg.get_frame().set_linewidth(0.25)
plt.savefig('build/plot_grav.pdf')
plt.close()
# Graphik zur Schwingung
par, cov = np.polyfit(1/B_schw, T_schw**2, deg=1, cov=True)
err = np.sqrt(np.diag(cov))
a_schw = ufloat(par[0],err[0])
with open('build/a_schw.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{a_schw.n * 1000:.2f}({a_schw.s * 1000:.2f})e-3')
	f.write(r'}{\tesla\second\squared}')
b_schw = ufloat(par[1],err[1])
with open('build/b_schw.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{b_schw.n:.2f}({b_schw.s:.2f})')
	f.write(r'}{\second\squared}')
x = np.linspace(250,750,2)
ax = plt.subplot(111)
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
plt.plot(x, par[0] * x + par[1], '-', c='olivedrab', label='Ausgleichsgerade')
plt.plot(1/B_schw, T_schw**2, 'kx', ms=3.21, label='Messergebnisse')
plt.xlabel(r'$B^{-1} \hspace{-0.5ex} \mathbin{/} \unit{\per\tesla}$')
plt.ylabel(r'$T^2 \hspace{-0.25ex} \mathbin{/} \unit{\second\squared}$')
leg = plt.legend(edgecolor='k', facecolor='none')
leg.get_frame().set_linewidth(0.25)
plt.savefig('build/plot_schw.pdf')
plt.close()
# Graphik zur Präzession
par_1, cov_1 = np.polyfit(B_präz, 1/T_präz, deg=1, cov=True)
err_1 = np.sqrt(np.diag(cov_1))
a_präz_1 = ufloat(par_1[0],err_1[0])
with open('build/a_präz_1.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{a_präz_1.n:.2f}({a_präz_1.s:.2f})')
	f.write(r'}{\per\second\per\tesla}')
b_präz_1 = ufloat(par_1[1],err_1[1])
with open('build/b_präz_1.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{b_präz_1.n * 1000:.2f}({b_präz_1.s * 1000:.2f})e-3')
	f.write(r'}{\per\second}')
par_2 = np.polyfit(np.take(B_präz,[2,10]), 1/np.take(T_präz,[2,10]), deg=1, cov=False)
a_präz_2 = par_2[0]
with open('build/a_präz_2.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{a_präz_2:.2f}')
	f.write(r'}{\per\second\per\tesla}')
b_präz_2 = par_2[1]
with open('build/b_präz_2.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{b_präz_2 * 1000:.2f}e-3')
	f.write(r'}{\per\second}')
par_3, cov_3 = np.polyfit(np.take(B_präz,[5,6,7]), 1/np.take(T_präz,[5,6,7]), deg=1, cov=True)
err_3 = np.sqrt(np.diag(cov_1))
a_präz_3 = ufloat(par_3[0],err_3[0])
with open('build/a_präz_3.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{a_präz_3.n:.2f}({a_präz_3.s:.2f})')
	f.write(r'}{\per\second\per\tesla}')
b_präz_3 = ufloat(par_3[1],err_3[1])
with open('build/b_präz_3.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{b_präz_3.n * 1000:.2f}({b_präz_3.s * 1000:.2f})e-3')
	f.write(r'}{\per\second}')
x = np.linspace(0.0006,0.005,2)
z = np.linspace(0.0028,0.0037,2)
ax = plt.subplot(111)
ax.xaxis.set_minor_locator(AutoMinorLocator(10))
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.plot(z, par_3[0] * z + par_3[1], '-', c='olivedrab', alpha=0.5)
plt.plot(x, par_1[0] * x + par_1[1], '-', c='olivedrab', alpha=0.5)
plt.plot(x, par_2[0] * x + par_2[1], '-', c='olivedrab', label='Ausgleichsgeraden')
plt.plot(B_präz[0:2], 1/T_präz[0:2], 'x', c='#808080', ms=3.21)
plt.plot(B_präz[3:10], 1/T_präz[3:10], 'x', c='#808080', ms=3.21)
plt.plot(B_präz[11:], 1/T_präz[11:], 'x', c='#808080', ms=3.21)
plt.plot(B_präz[2], 1/T_präz[2], 'kx', ms=3.21, label='Messergebnisse')
plt.plot(B_präz[10], 1/T_präz[10], 'kx', ms=3.21)
plt.xlabel(r'$B \mathbin{/} \unit{\tesla}$')
plt.ylabel(r'$T^{-1} \hspace{-0.5ex} \mathbin{/} \unit{\per\second}$')
leg = plt.legend(edgecolor='k', facecolor='none')
leg.get_frame().set_linewidth(0.25)
plt.savefig('build/plot_präz.pdf')
plt.close()
# Tabelle zur Gravitation
table_header = r'''	\begin{tabular}
		{S[table-format=2.0]
		 S[table-format=1.2]
		 S[table-format=1.2]}
		\toprule	
		{$r \mathbin{/} \unit{\milli\meter}$} &
		{$I \mathbin{/} \unit{\ampere}$} &
		{$B \mathbin{/} \unit{\milli\tesla}$} \\
		\midrule
'''
table_footer = r'''		\bottomrule
	\end{tabular}
'''
row_template = r'		{0:3.0f} & {1:1.2f} & {2:1.2f} \\'
with open('build/tab_grav.tex', 'w') as f:
	f.write(table_header)
	for row in zip(r_grav * 1000, I_grav, B_grav * 1000):
		f.write(row_template.format(*row))
		f.write('\n')
	f.write(table_footer)
# Tabelle zur Schwingung
table_header = r'''	\begin{tabular}
		{S[table-format=1.1]
		 S[table-format=1.2]
		 S[table-format=3.0]
		 S[table-format=1.2]
		 S[table-format=1.2]}
		\toprule	
		{$I \mathbin{/} \unit{\ampere}$} &
		{$B \mathbin{/} \unit{\milli\tesla}$} &
		{$B^{-1} \hspace{-0.5ex} \mathbin{/} \unit{\per\tesla}$} &
		{$T \hspace{-0.05ex} \mathbin{/} \unit{\second}$} &
		{$T^2 \hspace{-0.25ex} \mathbin{/} \unit{\second\squared}$} \\
		\midrule
'''
table_footer = r'''		\bottomrule
	\end{tabular}
'''
row_template = r'		{0:1.1f} & {1:1.2f} & {2:3.0f} & {3:1.2f} & {4:1.2f} \\[-0.15ex]'
with open('build/tab_schw.tex', 'w') as f:
	f.write(table_header)
	for row in sorted( zip(I_schw, B_schw * 1000, 1/B_schw, T_schw, T_schw**2)):
		f.write(row_template.format(*row))
		f.write('\n')
	f.write(table_footer)
# Tabelle zur Präzession
table_header = r'''	\begin{tabular}
		{S[table-format=1.2]
		 S[table-format=1.2]
		 S[table-format=2.2]
		 S[table-format=3.2]}
		\toprule	
		{$I \mathbin{/} \unit{\ampere}$} &
		{$B \mathbin{/} \unit{\milli\tesla}$} &
		{$T \hspace{-0.05ex} \mathbin{/} \unit{\second}$} &
		{$T^{-1} \hspace{-0.5ex} \mathbin{/} \unit{\per\second}$} \\
		\midrule
'''
table_footer = r'''		\bottomrule
	\end{tabular}
'''
row_template = r'		{0:1.2f} & {1:1.2f} & {2:2.2f} & {3:3.2f} \\[-0.15ex]'
with open('build/tab_präz.tex', 'w') as f:
	f.write(table_header)
	for row in sorted( zip(I_präz, B_präz * 1000, T_präz, 1/T_präz)):
		f.write(row_template.format(*row))
		f.write('\n')
	f.write(table_footer)
# Magnetische Momente berechnen
mu_grav = a_grav * m * g
with open('build/mu_grav.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{mu_grav.n:.3f}({mu_grav.s:.3f})')
	f.write(r'}{\ampere\meter\squared}')
mu_schw = 4 * np.pi**2 * J / a_schw
with open('build/mu_schw.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{mu_schw.n:.3f}({mu_schw.s:.3f})')
	f.write(r'}{\ampere\meter\squared}')
L_0 = J * 2 * np.pi * 6 # (f = 6 Hz)
with open('build/L_0.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{L_0 * 1000:.2f}e-3')
	f.write(r'}{\newton\meter\second}')
mu_präz_1 = a_präz_1 * 4 * np.pi**2 * J * 6
with open('build/mu_präz_1.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{mu_präz_1.n:.3f}({mu_präz_1.s:.3f})')
	f.write(r'}{\ampere\meter\squared}')
mu_präz_2 = a_präz_2 * 4 * np.pi**2 * J * 6
with open('build/mu_präz_2.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{mu_präz_2:.3f}')
	f.write(r'}{\ampere\meter\squared}')
mu_präz_3 = a_präz_3 * 4 * np.pi**2 * J * 6
with open('build/mu_präz_3.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{mu_präz_3.n:.3f}({mu_präz_3.s:.3f})')
	f.write(r'}{\ampere\meter\squared}')
f_präz_3 = mu_präz_2 / (a_präz_3 * 4 * np.pi**2 * J)
with open('build/f_präz_3.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{f_präz_3.n:.1f}({f_präz_3.s:.1f})')
	f.write(r'}{\hertz}')

