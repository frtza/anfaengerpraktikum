import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import uncertainties.unumpy as unp
from uncertainties.unumpy import nominal_values as noms, std_devs as stds
from uncertainties import ufloat

# Formeln
def I_kgl(m,R): # Kugel
	return 2/5 * m * R**2
def I_sym(m,R): # Zylinder
	return 1/2 * m * R**2
def I_ort(m,R,h): # Zylinder
	return m * (1/4 * R**2 + 1/12 * h**2)
def I_ges(D,T): # Experiment
	return D * T**2 / (4 * np.pi**2)

# Durchmesser d, Höhe h, Masse m für Kugel kgl, Zylinder zln, Gewichte 1 und 2
d_kgl = np.array([14.7,14.7,14.7]) # cm
m_kgl = np.array([1168.9,1169]) # g

d_zln = np.array([21.9]) # cm
h_zln = np.array([15]) # mm
m_zln = np.array([421.7]) # g

d_1 = np.array([34.9]) # mm
h_1 = np.array([30]) # mm
m_1 = np.array([222.7,222.7]) # g

d_2 = np.array([34.9]) # mm
h_2 = np.array([30]) # mm
m_2 = np.array([222.6,222.6]) # g

# Zylinder der Puppe mit Kopf k, Arme a, Torso t, Beine b
m = ufloat(168,0.1) # g

d_k = np.array([28.6,28.34,26.6,23.5,16.2]) # mm
h_k = np.array([54,53]) # mm

d_a = np.array([14.3,14.9,14.7,11.3,13.9,12.1,9.7,8.4,13.16,13.88]) # mm
h_a = np.array([137,136]) # mm

d_t = np.array([40.2,39.06,38.9,42.2,40.5,24.96,26.84,33.2,34.6,38.22]) # mm
h_t = np.array([99,99]) # mm

d_b = np.array([15,19.18,12.06,16.2,11.9,12.2,16.22,16.2,13.3,17.3]) # mm
h_b = np.array([170,169]) # mm

# Statische Messung stat bei Abstand a, Auslenkung phi, Kraft F
a_stat = 20 # cm
phi_stat, F_stat = np.genfromtxt('data/statisch.txt', unpack=True) # deg, N

# Dynamische Messung dyn bei Abstand a, Auslenkung phi, Periode T
phi_dyn = 90 # deg
a_dyn, f_T_dyn = np.genfromtxt('data/dynamisch.txt', unpack=True) # cm, s
T_dyn = f_T_dyn / 5 # s

# Perioden T für Kugel kgl, Zylinder zln bei Auslenkung von phi_dyn
f_T_kgl = np.genfromtxt('data/kugel.txt', unpack=True) # s
f_T_zln = np.genfromtxt('data/zylinder.txt', unpack=True) # s
T_kgl = f_T_kgl / 5 # s
T_zln = f_T_zln / 5 # s

# Puppe aus Zylindern in Posen A und B bei Auslenkung von phi 90, 120
phi_90, phi_120 = np.array([90]), np.array([120]) # deg, deg
f_T_A_90, f_T_A_120, f_T_B_90, f_T_B_120 = np.genfromtxt('data/puppe.txt', unpack=True) # s, s, s, s
T_A_90, T_A_120, T_B_90, T_B_120 = f_T_A_90 / 5, f_T_A_120 / 5, f_T_B_90 / 5, f_T_B_120 / 5 # s, s, s, s

# Mittlere Messwerte für Körper
d_kgl = ufloat(np.mean(d_kgl)*10,np.std(d_kgl)*10) # mm
m_kgl = ufloat(np.mean(m_kgl),0.1) # g

d_zln = ufloat(np.mean(d_zln)*10,np.std(d_zln)*10) # mm
h_zln = ufloat(np.mean(h_zln),np.std(h_zln)) # mm
m_zln = ufloat(np.mean(m_zln),0.1) # g

d_1 = ufloat(np.mean(d_1),np.std(d_1)) # mm
h_1 = ufloat(np.mean(h_1),np.std(h_1)) # mm
m_1 = ufloat(np.mean(m_1),0.1) # g

d_2 = ufloat(np.mean(d_2),np.std(d_2)) # mm
h_2 = ufloat(np.mean(h_2),np.std(h_2)) # mm
m_2 = ufloat(np.mean(m_2),0.1) # g

# Tabelle mit Dimensionen der Testmassen
table_header = r'''	\begin{tabular}
		{r
		 S[table-format=3.1]
		 @{${}\pm{}$}
		 S[table-format=1.1]
		 S[table-format=2.1]
		 S[table-format=2.0]}
		\toprule
		{$ $} &
		\multicolumn{2}{c}{$m \mathbin{/} \unit{\gram}$} &
		{$d \mathbin{/} \unit{\milli\meter}$} &
		{$h \mathbin{/} \unit{\milli\meter}$} \\
		\cmidrule[\lightrulewidth]{2-5}
'''
table_footer = r'''		\bottomrule
	\end{tabular}
'''
with open('build/dim_tab_test.tex', 'w') as f:
	f.write(table_header)
	f.write(r'		\textsc{1. Testmasse} & ')
	f.write(f'{m_1.n:.1f} & {m_1.s:.1f} & {d_1.n:.1f} & {h_1.n:.0f} ')
	f.write(r'\\')
	f.write('\n')
	f.write(r'		\textsc{2. Testmasse} & ')
	f.write(f'{m_2.n:.1f} & {m_2.s:.1f} & {d_2.n:.1f} & {h_2.n:.0f} ')
	f.write(r'\\')
	f.write('\n')
	f.write(table_footer)

# Tabelle mit Dimensionen der Körper
table_header = r'''	\begin{tabular}
		{r
		 S[table-format=4.1]
		 @{${}\pm{}$}
		 S[table-format=1.1]
		 S[table-format=3.0]
		 S[table-format=2.0]}
		\toprule
		{$ $} &
		\multicolumn{2}{c}{$m \mathbin{/} \unit{\gram}$} &
		{$d \mathbin{/} \unit{\milli\meter}$} &
		{$h \mathbin{/} \unit{\milli\meter}$} \\
		\cmidrule[\lightrulewidth]{2-5}
'''
table_footer = r'''		\bottomrule
	\end{tabular}
'''
with open('build/dim_tab_körper.tex', 'w') as f:
	f.write(table_header)
	f.write(r'		\textsc{Holzkugel} & ')
	f.write(f'{m_kgl.n:.1f} & {m_kgl.s:.1f} & {d_kgl.n:.0f} & {{-}} ')
	f.write(r'\\')
	f.write('\n')
	f.write(r'		\textsc{Holzzylinder} & ')
	f.write(f'{m_zln.n:.1f} & {m_zln.s:.1f} & {d_zln.n:.0f} & {h_zln.n:.0f} ')
	f.write(r'\\')
	f.write('\n')
	f.write(table_footer)

# Tabelle mit Messungen der Puppe
table_header = r'''	\begin{tabular}
		{S[table-format=2.2]
		 S[table-format=2.0]
		 S[table-format=2.2]
		 S[table-format=3.0]
		 S[table-format=2.2]
		 S[table-format=2.0]
		 S[table-format=2.2]
		 S[table-format=3.0]}
		\toprule
		\multicolumn{2}{c}{\textsc{Kopf}} &
		\multicolumn{2}{c}{\textsc{Arme}} &
		\multicolumn{2}{c}{\textsc{Torso}} &
		\multicolumn{2}{c}{\textsc{Beine}} \\
		\cmidrule(lr){1-2}\cmidrule(lr){3-4}\cmidrule(lr){5-6}\cmidrule(lr){7-8}
		{$d \mathbin{/} \unit{\milli\meter}$} &
		{$h \mathbin{/} \unit{\milli\meter}$} &
		{$d \mathbin{/} \unit{\milli\meter}$} &
		{$h \mathbin{/} \unit{\milli\meter}$} &
		{$d \mathbin{/} \unit{\milli\meter}$} &
		{$h \mathbin{/} \unit{\milli\meter}$} &
		{$d \mathbin{/} \unit{\milli\meter}$} &
		{$h \mathbin{/} \unit{\milli\meter}$} \\
		\midrule
'''
row_template = r'		{0:2.2f} & {1:2.0f} & {2:2.2f} & {3:3.0f} & {4:2.2f} & {5:2.0f} & {6:2.2f} & {7:3.0f} \\'
with open('build/mess_tab_puppe.tex', 'w') as f:
    f.write(table_header)
    for row in zip(d_k[:2], h_k[:2], d_a[:2], h_a[:2], d_t[:2], h_t[:2], d_b[:2], h_b[:2]):
        f.write(row_template.format(*row))
        f.write('\n')
row_template = r'		{0:2.2f} & {{-}} & {1:2.2f} & {{-}} & {2:2.2f} & {{-}} & {3:2.2f} & {{-}} \\'
with open('build/mess_tab_puppe.tex', 'a') as f:
    for row in zip(d_k[2:5], d_a[2:5], d_t[2:5], d_b[2:5]):
        f.write(row_template.format(*row))
        f.write('\n')
row_template = r'		{{-}} & {{-}} & {0:2.2f} & {{-}} & {1:2.2f} & {{-}} & {2:2.2f} & {{-}} \\'
with open('build/mess_tab_puppe.tex', 'a') as f:
    for row in zip(d_a[5:], d_t[5:], d_b[5:]):
        f.write(row_template.format(*row))
        f.write('\n')
    f.write(table_footer)

# Zylindervolumen
def V_Z(R, h):
	return np.pi * R**2 * h

# Mittlere Werte für Puppe
d_k = ufloat(np.mean(d_k),np.std(d_k)) # mm
h_k = ufloat(np.mean(h_k),np.std(h_k)) # mm
V_k = V_Z(d_k/2, h_k) / 1000 # cm**3

d_a = ufloat(np.mean(d_a),np.std(d_a)) # mm
h_a = ufloat(np.mean(h_a),np.std(h_a)) # mm
V_a = V_Z(d_a/2, h_a) / 1000 # cm**3

d_t = ufloat(np.mean(d_t),np.std(d_t)) # mm
h_t = ufloat(np.mean(h_t),np.std(h_t)) # mm
V_t = V_Z(d_t/2, h_t) / 1000 # cm**3

d_b = ufloat(np.mean(d_b),np.std(d_b)) # mm
h_b = ufloat(np.mean(h_b),np.std(h_b)) # mm
V_b = V_Z(d_b/2, h_b) / 1000 # cm**3

V = V_k + V_a + V_t + V_b # cm**3
rho = m / V # g / cm**3

m_k = rho * V_k # g
m_a = rho * V_a # g
m_t = rho * V_t # g
m_b = rho * V_b # g

# Tabelle mit Dimensionen der Puppe
with open('build/m_puppe.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{m.n:.1f}({m.s:.1f})')
	f.write(r'}{\gram}')
with open('build/V_puppe.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{V.n:.2f}({V.s:.2f})')
	f.write(r'}{\centi\meter\cubed}')
with open('build/rho_puppe.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{rho.n:.2f}({rho.s:.2f})')
	f.write(r'}{\gram\per\centi\meter\cubed}')
table_header = r'''	\sisetup{retain-zero-uncertainty}
	\begin{tabular}
		{r
		 S[table-format=2.2]
		 @{${}\pm{}$}
		 S[table-format=1.2]
		 S[table-format=3.1]
		 @{${}\pm{}$}
		 S[table-format=1.1]
		 S[table-format=3.2]
		 @{${}\pm{}$}
		 S[table-format=2.2]
		 S[table-format=2.2]
		 @{${}\pm{}$}
		 S[table-format=2.2]}
		\toprule
		{$ $} &
		\multicolumn{2}{c}{$d \mathbin{/} \unit{\milli\meter}$} &
		\multicolumn{2}{c}{$h \mathbin{/} \unit{\milli\meter}$} &
		\multicolumn{2}{c}{$V \mathbin{/} \unit{\centi\meter\cubed}$} &
		\multicolumn{2}{c}{$m \mathbin{/} \unit{\gram}$} \\
		\cmidrule[\lightrulewidth]{2-9}
'''
with open('build/dim_tab_puppe.tex', 'w') as f:
	f.write(table_header)
	f.write(r'		\textsc{Kopf} & ')
	f.write(f'{d_k.n:.2f} & {d_k.s:.2f} & {h_k.n:.1f} & {h_k.s:.1f} & {V_k.n:.2f} & {V_k.s:.2f} & {m_k.n:.2f} & {m_k.s:.2f} ')
	f.write(r'\\')
	f.write('\n')
	f.write(r'		\textsc{Arme} & ')
	f.write(f'{d_a.n:.2f} & {d_a.s:.2f} & {h_a.n:.1f} & {h_a.s:.1f} & {V_a.n:.2f} & {V_a.s:.2f} & {m_a.n:.2f} & {m_a.s:.2f} ')
	f.write(r'\\')
	f.write('\n')
	f.write(r'		\textsc{Torso} & ')
	f.write(f'{d_t.n:.2f} & {d_t.s:.2f} & {h_t.n:.1f} & {h_t.s:.1f} & {V_t.n:.2f} & {V_t.s:.2f} & {m_t.n:.2f} & {m_t.s:.2f} ')
	f.write(r'\\')
	f.write('\n')
	f.write(r'		\textsc{Beine} & ')
	f.write(f'{d_b.n:.2f} & {d_b.s:.2f} & {h_b.n:.1f} & {h_b.s:.1f} & {V_b.n:.2f} & {V_b.s:.2f} & {m_b.n:.2f} & {m_b.s:.2f} ')
	f.write(r'\\')
	f.write('\n')
	f.write(table_footer)

# Tabelle statische Messung Winkelrichtgröße
phi_stat_deg = phi_stat
phi_stat_rad = np.pi * phi_stat / 180
D = F_stat * a_stat / phi_stat_rad # N cm
with open('build/a_statisch.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{a_stat:.0f}')
	f.write(r'}{\centi\meter}')
table_header = r'''	\begin{tabular}
		{S[table-format=3.0]
		 c
		 S[table-format=1.3]
		 S[table-format=1.2]}
		\toprule
		{$\varphi \mathbin{/} \unit{\degree}$} &
		{$\varphi \mathbin{/} \unit{\radian}$} &
		{$F \mathbin{/} \unit{\newton}$} &
		{$D \mathbin{/} \qty{e-2}{\newton\meter}$} \\
		\midrule
'''
row_template = r'		{0:3.0f} & \num{{{1:1.2f}}}$\hspace{{0.15ex}}\pi$ & {2:1.3f} & {3:3.2f} \\[0.15ex]'
with open('build/mess_tab_statisch.tex', 'w') as f:
	f.write(table_header)
	for row in zip(phi_stat, phi_stat / 180, F_stat, D):
		f.write(row_template.format(*row))
		f.write('\n')
	f.write(table_footer)
D = ufloat(np.mean(D),np.std(D)) # N cm
with open('build/D.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{D.n:.2f}({D.s:.2f})e-2')
	f.write(r'}{\newton\meter}')

# Tabelle dynamische Messung Trägheitsmoment
with open('build/phi_dynamisch.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{phi_dyn:.0f}')
	f.write(r'}{\degree}')
table_header = r'''	\begin{tabular}
		{S[table-format=2.0]
		 S[table-format=3.0]
		 S[table-format=1.2]
		 S[table-format=2.2]}
		\toprule
		{$a \mathbin{/} \unit{\centi\meter}$} &
		{$a^2\hspace{-0.15ex} \mathbin{/} \unit{\centi\meter\squared}$} &
		{$T\hspace{-0.05ex} \mathbin{/} \unit{\second}$} &
		{$T^{\hspace{0.15ex}2}\hspace{-0.15ex} \mathbin{/} \unit{\second\squared}$} \\
		\midrule
'''
row_template = r'		{0:2.0f} & {1:3.0f} & {2:1.2f} & {3:2.2f} \\'
with open('build/mess_tab_dynamisch.tex', 'w') as f:
	f.write(table_header)
	for row in zip(a_dyn, a_dyn**2, T_dyn, T_dyn**2):
		f.write(row_template.format(*row))
		f.write('\n')
	f.write(table_footer)

# Plot der dynamischen Messergebnisse
par, cov = np.polyfit((a_dyn/100)**2, T_dyn**2, deg=1, cov=True)
err = np.sqrt(np.diag(cov))
with open('build/p_dynamisch.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{par[0]:.0f}({err[0]:.0f})')
	f.write(r'}{\second\squared\per\meter\squared}')
with open('build/q_dynamisch.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{par[1]:.2f}({err[1]:.2f})')
	f.write(r'}{\second\squared}')
ax = plt.subplot(111)
ax.xaxis.set_minor_locator(AutoMinorLocator(5))
ax.yaxis.set_minor_locator(AutoMinorLocator(5))
ax.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
x = np.linspace(0,0.08,2)
plt.fill_between(x, (par[0] + 3*err[0])*x + (par[1] + 3*err[1]), (par[0] - 3*err[0])*x + (par[1] - 3*err[1]),
				 alpha=0.4, facecolor='olivedrab')
plt.plot(x, par[0]*x + par[1], '-', lw=1, c='olivedrab', label='Ausgleichsgerade')
plt.plot((a_dyn/100)**2, T_dyn**2, 'kx', ms=3.21, label='Messdaten')
plt.xlabel(r'$a^2\hspace{-0.15ex} \mathbin{/} \unit{\meter\squared}$')
plt.ylabel(r'$T^{\hspace{0.15ex}2}\hspace{-0.15ex} \mathbin{/} \unit{\second\squared}$')
leg = plt.legend(edgecolor='k', facecolor='none')
leg.get_frame().set_linewidth(0.25)
plt.xlim(0.00202,0.071979999999999995)
plt.ylim(5.200887714512961,50.21541755235693)
plt.savefig('build/plot.pdf')

# Berechnung Trägheitsmoment des Stabes
m_stab = ufloat(95.8,0.1) # g
with open('build/m_stab.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{m_stab.n:.1f}({m_stab.s:.1f})')
	f.write(r'}{\gram}')
h_stab = ufloat(60,0) # cm
with open('build/h_stab.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{h_stab.n:.0f}({h_stab.s:.0f})')
	f.write(r'}{\centi\meter}')
I_stab = 1/12 * m_stab * (h_stab/10)**2 # g dm**2
with open('build/I_stab.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{I_stab.n:.1f}({I_stab.s:.1f})e-5')
	f.write(r'}{\kilo\gram\meter\squared}')

# Berechnung Trägheitsmoment Testmassen
m = (m_1 + m_2) / 2 # g
with open('build/m.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{m.n:.1f}({m.s:.1f})')
	f.write(r'}{\gram}')
R = (d_1 + d_2) / 4 # mm
with open('build/R.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{R.n:.1f}({R.s:.1f})')
	f.write(r'}{\milli\meter}')
h = (h_1 + h_2) / 2 # mm
with open('build/h.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{h.n:.0f}({h.s:.0f})')
	f.write(r'}{\milli\meter}')
I_0 = I_ort(m, R/10, h/10)
with open('build/I_0.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{I_0.n:.1f}({I_0.s:.1f})e-7')
	f.write(r'}{\kilo\gram\meter\squared}')

# Berechnung Eigenträgheitsmoment
q = ufloat(par[1],err[1]) # s**2
I_D = 1/(4 * np.pi**2) * D * q * 1000 - 2 * I_0/100 # g dm**2
with open('build/I_D.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{I_D.n:.1f}({I_D.s:.1f})e-5')
	f.write(r'}{\kilo\gram\meter\squared}')
I_D_korr = I_D - I_stab # g dm**2
with open('build/I_D_korr.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{I_D_korr.n:.1f}({I_D_korr.s:.1f})e-5')
	f.write(r'}{\kilo\gram\meter\squared}')


# Tabelle zu Schwingzeiten der Körper
table_header = r'''	\begin{tabular}
		{S[table-format=1.2]
		 S[table-format=1.2]
		 S[table-format=2.2]
		 S[table-format=1.3]}
		\toprule	
		\multicolumn{4}{c}{$T\hspace{-0.05ex} \mathbin{/} \unit{\second}$} \\
		\midrule
		\multicolumn{2}{c}{\textsc{Holzkugel}} &
		\multicolumn{2}{c}{\textsc{Holzzylinder}} \\
		\cmidrule(lr){1-2}\cmidrule(lr){3-4}
'''
row_template = r'		{0:1.2f} & {1:1.2f} & {2:1.2f} & {3:1.2f} \\'
with open('build/period_tab_körper.tex', 'w') as f:
	f.write(table_header)
	for row in zip(T_kgl[:5], T_kgl[5:], T_zln[:5], T_zln[5:]):
		f.write(row_template.format(*row))
		f.write('\n')
	f.write(table_footer)

# Mittlere Perioden der Körper
T_kgl = ufloat(np.mean(T_kgl),np.std(T_kgl))
with open('build/T_kgl.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{T_kgl.n:.2f}({T_kgl.s:.2f})')
	f.write(r'}{\second}')
T_zln = ufloat(np.mean(T_zln),np.std(T_zln))
with open('build/T_zln.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{T_zln.n:.2f}({T_zln.s:.2f})')
	f.write(r'}{\second}')

# Theoretische Trägheitsmomente der Körper
I_kgl_theo = I_kgl(m_kgl, d_kgl/200)
with open('build/I_kgl_theo.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{I_kgl_theo.n:.2f}({I_kgl_theo.s:.2f})e-5')
	f.write(r'}{\kilo\gram\meter\squared}')
I_zln_theo = I_sym(m_zln, d_zln/200)
with open('build/I_zln_theo.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{I_zln_theo.n:.2f}({I_zln_theo.s:.2f})e-5')
	f.write(r'}{\kilo\gram\meter\squared}')

# Experimentelle Trägheitsmomente der Körper
I_kgl_exp = I_ges(1000*D, T_kgl) - I_D
with open('build/I_kgl_exp.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{I_kgl_exp.n:.2f}({I_kgl_exp.s:.2f})e-5')
	f.write(r'}{\kilo\gram\meter\squared}')
I_zln_exp = I_ges(1000*D, T_zln) - I_D
with open('build/I_zln_exp.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{I_zln_exp.n:.2f}({I_zln_exp.s:.2f})e-5')
	f.write(r'}{\kilo\gram\meter\squared}')
I_kgl_exp_korr = I_ges(1000*D, T_kgl) - I_D_korr
with open('build/I_kgl_exp_korr.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{I_kgl_exp_korr.n:.2f}({I_kgl_exp_korr.s:.2f})e-5')
	f.write(r'}{\kilo\gram\meter\squared}')
I_zln_exp_korr = I_ges(1000*D, T_zln) - I_D_korr
with open('build/I_zln_exp_korr.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{I_zln_exp_korr.n:.2f}({I_zln_exp_korr.s:.2f})e-5')
	f.write(r'}{\kilo\gram\meter\squared}')

# Tabelle zu Schwingzeiten der Puppe
table_header = r'''	\begin{tabular}
		{S[table-format=1.2]
		 S[table-format=1.2]
		 S[table-format=1.2]
		 S[table-format=1.2]}
		\toprule	
		\multicolumn{4}{c}{$T\hspace{-0.05ex} \mathbin{/} \unit{\second}$} \\
		\midrule
		\multicolumn{2}{c}{\textsc{1. Körperhaltung}} &
		\multicolumn{2}{c}{\textsc{2. Körperhaltung}} \\
		\cmidrule(lr){1-2}\cmidrule(lr){3-4}
		{$\varphi = \qty{90}{\degree}$} &
		{$\varphi = \qty{120}{\degree}$} &
		{$\varphi = \qty{90}{\degree}$} &
		{$\varphi = \qty{120}{\degree}$} \\
		\cmidrule(lr){1-1}\cmidrule(lr){2-2}\cmidrule(lr){3-3}\cmidrule(lr){4-4}
'''
row_template = r'		{0:1.2f} & {1:1.2f} & {2:1.2f} & {3:1.2f} \\'
with open('build/period_tab_puppe.tex', 'w') as f:
	f.write(table_header)
	for row in zip(T_A_90, T_A_120, T_B_90, T_B_120):
		f.write(row_template.format(*row))
		f.write('\n')
	f.write(table_footer)

# Gemittelte Perioden der Puppe
T_1 = np.concatenate((T_A_90,T_A_120))
T_1 = ufloat(np.mean(T_1),np.std(T_1))
with open('build/T_1.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{T_1.n:.2f}({T_1.s:.2f})')
	f.write(r'}{\second}')
T_2 = np.concatenate((T_B_90,T_B_120))
T_2 = ufloat(np.mean(T_2),np.std(T_2))
with open('build/T_2.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{T_2.n:.2f}({T_2.s:.2f})')
	f.write(r'}{\second}')

# Theoretisches Trägheitsmoment in Haltung 1
I_k = I_sym(m_k, d_k/200)
I_a = I_ort(m_a, d_a/200, h_a/100) + m_a * (h_a/200 + d_t/200)**2
I_t = I_sym(m_t, d_t/200)
I_b = I_sym(m_b, d_b/200) + m_b * (d_b/200)**2
I_1_theo = I_k + I_a + I_t + I_b
with open('build/I_1_theo.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{I_1_theo.n:.2f}({I_1_theo.s:.2f})e-5')
	f.write(r'}{\kilo\gram\meter\squared}')

# Experimentelles Trägheitsmoment in Haltung 1
I_1_exp = I_ges(1000*D, T_1) - I_D
with open('build/I_1_exp.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{I_1_exp.n:.2f}({I_1_exp.s:.2f})e-5')
	f.write(r'}{\kilo\gram\meter\squared}')
I_1_exp_korr = I_ges(1000*D, T_1) - I_D_korr
with open('build/I_1_exp_korr.tex', 'w') as f:
	f.write(r'(\num{')
	f.write(f'{I_1_exp_korr.n:.2f}')
	f.write(r'}\pm\num{')
	f.write(f'{I_1_exp_korr.s:.2f}')
	f.write(r'}) \cdot \qty{e-5}{\kilo\gram\meter\squared}')

# Theoretisches Trägheitsmoment in Haltung 2
I_k = I_sym(m_k, d_k/200)
I_a = I_ort(m_a, d_a/200, h_a/100) + m_a * (h_a/200 + d_t/200)**2
I_t = I_sym(m_t, d_t/200)
I_b = I_sym(m_b, d_b/200) + m_b * ((d_b/200)**2 + (h_b/200)**2)
I_2_theo = I_k + I_a + I_t + I_b
with open('build/I_2_theo.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{I_2_theo.n:.2f}({I_2_theo.s:.2f})e-5')
	f.write(r'}{\kilo\gram\meter\squared}')

# Experimentelles Trägheitsmoment in Haltung 2
I_2_exp = I_ges(1000*D, T_2) - I_D
with open('build/I_2_exp.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{I_2_exp.n:.2f}({I_2_exp.s:.2f})e-5')
	f.write(r'}{\kilo\gram\meter\squared}')
I_2_exp_korr = I_ges(1000*D, T_2) - I_D_korr
with open('build/I_2_exp_korr.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{I_2_exp_korr.n:.2f}({I_2_exp_korr.s:.2f})e-5')
	f.write(r'}{\kilo\gram\meter\squared}')

# Relative Abweichungen der Puppe
rel_kgl = I_kgl_exp_korr.n / I_kgl_theo.n * 100
with open('build/rel_kgl.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{rel_kgl:.1f}')
	f.write(r'}{\percent}')
rel_zln = I_zln_exp_korr.n / I_zln_theo.n * 100
with open('build/rel_zln.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{rel_zln:.1f}')
	f.write(r'}{\percent}')

# Relative Abweichungen der Puppe
rel_1 = I_1_exp_korr.n / I_1_theo.n * 100
with open('build/rel_1.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{rel_1:.1f}')
	f.write(r'}{\percent}')
rel_1_max = I_1_exp_korr.n + I_1_exp_korr.s / I_1_theo.n * 100
with open('build/rel_1_max.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{rel_1_max:.1f}')
	f.write(r'}{\percent}')
rel_2 = I_2_exp_korr / I_2_theo * 100
with open('build/rel_2.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{rel_2.n:.1f}')
	f.write(r'}{\percent}')
rel_theo = I_1_theo.n / I_2_theo.n * 100
with open('build/rel_theo.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{rel_theo:.1f}')
	f.write(r'}{\percent}')
rel_exp = I_1_exp_korr.n / I_2_exp_korr.n * 100
with open('build/rel_exp.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{rel_exp:.1f}')
	f.write(r'}{\percent}')
rel_exp_max = I_1_exp_korr.n + I_1_exp_korr.s / I_2_exp_korr.n * 100
with open('build/rel_exp_max.tex', 'w') as f:
	f.write(r'\qty{')
	f.write(f'{rel_exp_max:.1f}')
	f.write(r'}{\percent}')

