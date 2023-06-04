# bibliotheken laden
import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
from uncertainties.unumpy import nominal_values as noms, std_devs as stds
from uncertainties import ufloat
from scipy.constants import convert_temperature
from scipy.optimize import curve_fit
from matplotlib.ticker import AutoMinorLocator
from matplotlib import container

# dimensionen in mm und g einlesen
d_k, m_k, d_g, m_g = np.genfromtxt('data/dimension.txt', unpack=True)

# tabelle mit messwerten
table_header = r'''	\begin{tabular}
		{S[table-format=2.2]
		 S[table-format=1.2]
		 c
		 S[table-format=2.2]
		 S[table-format=1.2]}
		\toprule
		\multicolumn{2}{c}{\textsc{Kleine Kugel}} & & \multicolumn{2}{c}{\textsc{Große Kugel}} \\
		\cmidrule(lr){1-2}\cmidrule(lr){4-5}
		{$d \mathbin{/} \unit{\milli\meter}$} &
		{$m \mathbin{/} \unit{\gram}$} &
		&
		{$d \mathbin{/} \unit{\milli\meter}$} &
		{$m \mathbin{/} \unit{\gram}$} \\
		\midrule
'''
table_footer = r'''		\bottomrule
	\end{tabular}
'''
row_template = r'		{0:2.2f} & {1:1.2f} & & {2:2.2f} & {3:1.2f} \\'
with open('build/meas_dim_tab.tex', 'w') as f:
	f.write(table_header)
	for row in zip(d_k, m_k, d_g, m_g):
		f.write(row_template.format(*row))
		f.write('\n')
	f.write(table_footer)

# fehler der dimensionen bestimmen
d_k = ufloat(np.mean(d_k),np.std(d_k))
m_k = ufloat(np.mean(m_k),np.std(m_k))
d_g = ufloat(np.mean(d_g),np.std(d_g))
m_g = ufloat(np.mean(m_g),np.std(m_g))

# volumen in mm^3 berechnen
V_k = 4/3 * np.pi * (d_k/2)**3
V_g = 4/3 * np.pi * (d_g/2)**3

# dichte in g/mm^3 berechnen
rho_k = m_k / V_k
rho_g = m_g / V_g

# tabelle mit berechneten dimensionen
table_header = r'''	\sisetup{retain-zero-uncertainty}
	\begin{tabular}{r c c c c}
		\toprule
		{$ $} &
		{$d \mathbin{/} \unit{\milli\meter}$} &
		{$m \mathbin{/} \unit{\gram}$} &
		{$V \mathbin{/} \unit{\centi\meter\cubed}$} &
		{$\rho \mathbin{/} \unit{\gram\per\centi\meter\cubed}$} \\
		\cmidrule[\lightrulewidth]{2-5}
'''
with open('build/calc_dim_tab.tex', 'w') as f:
	f.write(table_header)
	f.write(r'		\textsc{Kleine Kugel} & ')
	f.write(r'\num{')
	f.write(f'{d_k.n:.2f}({d_k.s:.2f})')
	f.write(r'} & \num{')
	f.write(f'{m_k.n:.2f}({m_k.s:.2f})')
	f.write(r'} & \num{')
	f.write(f'{V_k.n/1000:.2f}({V_k.s/1000:.2f})')
	f.write(r'} & \num{')
	f.write(f'{1000*rho_k.n:.2f}({1000*rho_k.s:.2f})')
	f.write(r'} \\')
	f.write('\n')
	f.write(r'		\textsc{Große Kugel} & ')
	f.write(r'\num{')
	f.write(f'{d_g.n:.2f}({d_g.s:.2f})')
	f.write(r'} & \num{')
	f.write(f'{m_g.n:.2f}({m_g.s:.2f})')
	f.write(r'} & \num{')
	f.write(f'{V_g.n/1000:.2f}({V_g.s/1000:.2f})')
	f.write(r'} & \num{')
	f.write(f'{1000*rho_g.n:.2f}({1000*rho_g.s:.2f})')
	f.write(r'} \\')
	f.write('\n')
	f.write(table_footer)

# kleine kugel zeit einlesen
t_k, T_k = np.genfromtxt('data/klein_konstant.txt', unpack=True)

# große kugel zeit einlesen
t_g, T_g = np.genfromtxt('data/groß_konstant.txt', unpack=True)

# tabelle kugel fallzeit konstante temperatur
table_header = r'''	\begin{tabular}
		{S[table-format=2.2]
		 S[table-format=2.2]
		 c
		 S[table-format=2.2]
		 S[table-format=3.2]}
		\toprule
		\multicolumn{5}{c}{$t \mathbin{/} \unit{\second}$} \\
		\midrule
		\multicolumn{2}{c}{\textsc{Kleine Kugel}} & & \multicolumn{2}{c}{\textsc{Große Kugel}} \\
		\cmidrule(lr){1-2}\cmidrule(lr){4-5}
'''
row_template = r'		{0:2.2f} & {1:2.2f} & & {2:2.2f} & {3:3.2f} \\'
with open('build/const_tab.tex', 'w') as f:
	f.write(table_header)
	for row in zip(t_k[:5], t_k[5:], np.take(t_g,[0,1,4,8,9]), np.take(t_g,[5,6,7,2,3])):
		f.write(row_template.format(*row))
		f.write('\n')
	f.write(table_footer)

# fehler kleine kugel bestimmen
t_kl = ufloat(np.mean(t_k),np.std(t_k))
T_k = ufloat(np.mean(T_k),np.std(T_k))

# kleine kugel werte schreiben
with open('build/t_kl.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{t_kl.n:.2f}({t_kl.s:.2f})')
	f.write(r'}{\second}')
with open('build/T_kl.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{T_k.n:.0f}({T_k.s:.0f})')
	f.write(r'}{\celsius}')

# fehler große kugel bestimmen
t_gr = ufloat(np.mean(t_g),np.std(t_g))
T_g = ufloat(np.mean(T_g),np.std(T_g))

# große kugel werte schreiben
with open('build/t_gr.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{t_gr.n:.2f}({t_gr.s:.2f})')
	f.write(r'}{\second}')
with open('build/T_gr.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{T_g.n:.0f}({T_g.s:.0f})')
	f.write(r'}{\celsius}')

# temperaturbedingter verlauf der fallzeit für die große kugel erstellen
T_C, t_1, t_2 = np.genfromtxt('data/groß_zeit.txt', unpack=True)
t = unp.uarray(np.mean([t_1, t_2], axis = 0), np.std([t_1, t_2], axis = 0))

# tabelle große kugel fallzeit variable temperatur
table_header = r'''	\begin{tabular}
		{S[table-format=2.0]
		 S[table-format=2.2]
		 S[table-format=2.2]
		 @{\hspace{6ex}}
		 S[table-format=2.0]
		 S[table-format=2.2]
		 S[table-format=2.2]
		 @{\hspace{6ex}}
		 S[table-format=2.0]
		 S[table-format=2.2]
		 S[table-format=2.2]}
		\toprule
		{$T \mathbin{/} \unit{\celsius}$} &
		\multicolumn{2}{c}{$t \mathbin{/} \unit{\second}$} \hspace{3ex} &
		{$T \mathbin{/} \unit{\celsius}$} &
		\multicolumn{2}{c}{$t \mathbin{/} \unit{\second}$} \hspace{3ex} &
		{$T \mathbin{/} \unit{\celsius}$} &
		\multicolumn{2}{c}{$t \mathbin{/} \unit{\second}$} \\
		\midrule
'''
row_template = r'		{0:2.0f} & {1:2.2f} & {2:2.2f} & {3:2.0f} & {4:2.2f} & {5:2.2f} & {6:2.0f} & {7:2.2f} & {8:2.2f} \\'
with open('build/var_tab.tex', 'w') as f:
	f.write(table_header)
	for row in zip(T_C[:3], t_1[:3], t_2[:3], T_C[3:6], t_1[3:6], t_2[3:6], T_C[6:], t_1[6:], t_2[6:]):
		f.write(row_template.format(*row))
		f.write('\n')
	f.write(table_footer)

# konstanten schreiben
m_kl = 4.4531
with open('build/m_kl.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{m_kl:1.4f}')
	f.write(r'}{\gram}')
K_kl = 0.0764
with open('build/K_kl.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{K_kl:1.4f}')
	f.write(r'}{\milli\pascal\centi\meter\cubed\per\gram}')

# viskosität berechnen
rho_fl = np.array([0.998595,
				   0.997991,
				   0.997295,
				   0.996511,
				   0.995645,
				   0.994700,
				   0.993681,
				   0.992591,
				   0.990208,
				   0.988030])
with open('build/rho_fl.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{rho_fl[0]:1.6f}')
	f.write(r'}{\gram\per\centi\meter\cubed}')
rho_kl = 1000 * rho_k
eta = K_kl * (rho_kl - rho_fl[0]) * t_kl
with open('build/eta.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{eta.n:.3f}({eta.s:.3f})')
	f.write(r'}{\milli\pascal\second}')

# konstante große kugel berechnen
rho_gr = 1000 * rho_g
K_gr = eta / ((rho_gr - rho_fl[0]) * t_gr)
with open('build/K_gr.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{K_gr.n:.3f}({K_gr.s:.3f})')
	f.write(r'}{\milli\pascal\centi\meter\cubed\per\gram}')

# arrays ergänzen
T_C = np.insert(T_C, 0, T_g.n)
t = np.insert(t, 0, t_gr)

# array konvertieren
T_K = convert_temperature(T_C, 'C', 'K')

# arrays berechnen
eta = K_gr * (rho_gr - rho_fl) * t
v = 100 / t
Re = (rho_fl * 1000) * (v / 1000) * (d_g / 1000) / (eta / 1000)

# reynoldszahl kleine kugel
v_kl = 100 / t_kl
with open('build/v_kl.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{v_kl.n:.2f}({v_kl.s:.2f})')
	f.write(r'}{\milli\meter\per\second}')
Re_kl = (rho_fl[0] * 1000) * (v_kl / 1000) * (d_k / 1000) / (eta[0] / 1000)
with open('build/Re_kl.tex', 'w') as f:
	f.write(r'\num{')
	f.write(f'{Re_kl.n:.2f}({Re_kl.s:.2f})')
	f.write(r'}')

# tabelle große kugel komplett
table_header = r'''	\begin{tabular}
		{S[table-format=2.0]
		 S[table-format=3.2]
		 @{${}\pm{}$}
		 S[table-format=1.2]
		 S[table-format=1.6]
		 S[table-format=1.3]
		 @{${}\pm{}$}
		 S[table-format=1.3]
		 S[table-format=1.2]
		 @{${}\pm{}$}
		 S[table-format=1.2]
		 S[table-format=2.2]
		 @{${}\pm{}$}
		 S[table-format=1.2]}
		\toprule
		{$T \mathbin{/} \unit{\celsius}$} &
		\multicolumn{2}{c}{$t \mathbin{/} \unit{\second}$} &
		{$\rho_\text{f\hspace{0.15ex}l} \mathbin{/} \unit{\gram\per\centi\meter\cubed}$} &
		\multicolumn{2}{c}{$\eta \mathbin{/} \unit{\milli\pascal\second}$} &
		\multicolumn{2}{c}{$v \mathbin{/} \unit{\milli\meter\per\second}$} &
		\multicolumn{2}{c}{$\symit{Re}$} \\
		\midrule
'''
row_template = r'		{0:2.0f} & {1:3.2f} & {2:1.2f} & {3:1.6f} & {4:1.3f} & {5:1.3f} & {6:1.2f} & {7:1.2f} & {8:2.2f} & {9:1.2f} \\'
with open('build/tab.tex', 'w') as f:
	f.write(table_header)
	for row in zip(T_C, noms(t), stds(t), rho_fl, noms(eta), stds(eta), noms(v), stds(v), noms(Re), stds(Re)):
		f.write(row_template.format(*row))
		f.write('\n')
	f.write(table_footer)

# literaturwerte einlesen
T_lit, eta_lit = np.genfromtxt('data/literatur.txt', unpack=True)

# fit funktion
def fit(T, a, b, c, d):
	return a + b/T + c*T + d*T**2

# curve fits
par, cov = curve_fit(fit, T_C, noms(eta))
err = np.sqrt(np.diag(cov))
par_lit, cov_lit = curve_fit(fit, T_lit, eta_lit)
err_lit = np.sqrt(np.diag(cov_lit))

# export help properties
with open('build/a.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{par[0]:.2f}({err[0]:.2f})')
	f.write(r'}{\milli\pascal\second}')
with open('build/b.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{par[1]:.2f}({err[1]:.2f})')
	f.write(r'}{\milli\pascal\second\celsius}')
with open('build/c.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{par[2]:.2f}({err[2]:.2f})')
	f.write(r'}{\milli\pascal\second\per\celsius}')
with open('build/d.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{1000*par[3]:.2f}({1000*err[3]:.2f})')
	f.write(r'}{\micro\pascal\second\per\celsius\squared}')
with open('build/a_lit.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{par_lit[0]:.2f}({err_lit[0]:.2f})')
	f.write(r'}{\milli\pascal\second}')
with open('build/b_lit.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{par_lit[1]:.2f}({err_lit[1]:.2f})')
	f.write(r'}{\milli\pascal\second\celsius}')
with open('build/c_lit.tex', 'w') as f:
	f.write(r'\qty[retain-zero-uncertainty]{')
	f.write(f'{par_lit[2]:.2f}({err_lit[2]:.2f})')
	f.write(r'}{\milli\pascal\second\per\celsius}')
with open('build/d_lit.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{1000*par_lit[3]:.2f}({1000*err_lit[3]:.2f})')
	f.write(r'}{\micro\pascal\second\per\celsius\squared}')

# graphischer viskositätsverlauf
T_x = np.linspace(17,51,10000)
ax = plt.subplot(111)
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(4))
plt.plot(T_x, fit(T_x, *par), '-', c='olivedrab', alpha=0.4)
plt.plot(T_x, fit(T_x, *par_lit), '-', c='steelblue', alpha=0.4)
plt.errorbar(T_C, noms(eta), yerr=stds(eta), fmt='.', ms=4.5, c='olivedrab', label='Messdaten')
plt.plot(T_lit[1:5], eta_lit[1:5], 'x', ms=3.21, c='steelblue', label='Literaturwerte')
plt.xlabel(r'$T \mathbin{/} \unit{\celsius}$')
plt.ylabel(r'$\eta \mathbin{/} \unit{\milli\pascal\second}$')
handles, labels = ax.get_legend_handles_labels()
handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles]
leg = plt.legend(handles, labels, edgecolor='k', facecolor='none')
leg.get_frame().set_linewidth(0.25)
plt.savefig('build/plot_eta.pdf')
plt.close()

# ausgleichgerade bestimmen
par, cov = np.polyfit(1/T_K, np.log(noms(eta)), deg=1, cov=True)
err = np.sqrt(np.diag(cov))
A = unp.exp(ufloat(par[1], err[1]))
B = par[0]

# parameter exportieren
with open('build/A.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{1000*noms(A):.2f}({1000*stds(A):.2f})')
	f.write(r'}{\micro\pascal\second}')
with open('build/B.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{B:.2f}({err[0]:.2f})')
	f.write(r'}{\kelvin}')

# dimensionslose viskosität logarithmisieren
eta = unp.log(eta/noms(A))

# linearisierte darstellung
T_x = np.linspace(0.00308,0.00345,10000)
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
ax = plt.subplot(111)
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
plt.plot(T_x, B*T_x, '-', c='olivedrab', alpha=0.4, label='Regression')
plt.errorbar(1/T_K, noms(eta), yerr=stds(eta), fmt='.', ms=4.5, c='olivedrab', label='Messdaten')
plt.xlabel(r'$T^{\hspace{0.15ex}-1} \mathbin{/} \unit{\per\kelvin}$')
plt.ylabel(r'$\ln \left( \eta / \hspace{-0.3ex} A \right)$')
handles, labels = ax.get_legend_handles_labels()
handles = [h[0] if isinstance(h, container.ErrorbarContainer) else h for h in handles]
leg = plt.legend(handles, labels, edgecolor='k', facecolor='none')
leg.get_frame().set_linewidth(0.25)
plt.savefig('build/plot_lin.pdf')
plt.close()
