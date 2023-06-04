import matplotlib.pyplot as plt
import numpy as np
import uncertainties.unumpy as unp
from matplotlib.ticker import AutoMinorLocator
from scipy.optimize import curve_fit

#################################################################################################################
# Zeitabhängige Entladungskurve

t, U_t = np.genfromtxt('data/t-U_data.txt', unpack=True) 
U_t_0 = 3.2

par_t_U, cov_t_U = np.polyfit(t, np.log(U_t/U_t_0), deg=1, cov=True)
err_t_U = np.sqrt(np.diag(cov_t_U))

val_t_U = unp.uarray(par_t_U, err_t_U)

x = np.linspace(-0.02,0.70,10000)

plt.figure(figsize=(5.78,2.57))
ax = plt.subplot(1,2,1)
ax.xaxis.set_minor_locator(AutoMinorLocator(3))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
f, = plt.plot(x, par_t_U[0] * x + par_t_U[1], '-', c='olivedrab')
p, = plt.plot(t, np.log(U_t/U_t_0), 'kx', ms=3.21)
plt.xlabel(r'$t \mathbin{/} \unit{\milli\second}$')
plt.ylabel(r'$\ln(U\hspace{-0.2ex}/U_{\! 0})$')
leg = plt.legend([p, f], ['Messdaten', 'Regression'], edgecolor='k', facecolor='none')
leg.get_frame().set_linewidth(0.25)

ax = plt.subplot(1,2,2)
ax.xaxis.set_minor_locator(AutoMinorLocator(3))
ax.yaxis.set_minor_locator(AutoMinorLocator(3))
f, = plt.plot(x, U_t_0 * np.exp(par_t_U[0] * x + par_t_U[1]), '-', c='olivedrab')
p, = plt.plot(t, U_t, 'kx', ms=3.21)
plt.xlabel(r'$t \mathbin{/} \unit{\milli\second}$')
plt.ylabel(r'$U \hspace{-0.2ex} \mathbin{/} \unit{\volt}$')
leg = plt.legend([p, f], ['Messdaten', 'Regression'], edgecolor='k', facecolor='none')
leg.get_frame().set_linewidth(0.25)

plt.savefig('build/t-U_plot.pdf')
plt.close()

with open('build/t-U_0.tex', 'w') as file:
	file.write(r'\qty{')
	file.write(f'{U_t_0:.1f}')
	file.write(r'}{\volt}')
with open('build/t-U_v.tex', 'w') as file:
	file.write(r'\qty{')
	file.write(f'{val_t_U[0].n:.2f}({val_t_U[0].s:.2f})')
	file.write(r'}{\per\milli\second}')
with open('build/t-U_w.tex', 'w') as file:
	file.write(r'\num{')
	file.write(f'{val_t_U[1].n:.2f}({val_t_U[1].s:.2f})')
	file.write(r'}')
with open('build/t-U_RC.tex', 'w') as file:
	file.write(r'\qty{')
	file.write(f'{(-1000/val_t_U[0]).n:.1f}({(-1000/val_t_U[0]).s:.1f})')
	file.write(r'}{\micro\second}')

upper_par_t_U, upper_cov_t_U = np.polyfit(1.25*t, np.log(U_t/U_t_0), deg=1, cov=True)
upper_err_t_U = np.sqrt(np.diag(upper_cov_t_U))

upper_val_t_U = unp.uarray(upper_par_t_U, upper_err_t_U)

with open('build/t-U_RC_upper.tex', 'w') as file:
	file.write(r'\qty{')
	file.write(f'{(-1000/upper_val_t_U[0]).n:.1f}({(-1000/upper_val_t_U[0]).s:.1f})')
	file.write(r'}{\micro\second}')

lower_par_t_U, lower_cov_t_U = np.polyfit(0.75*t, np.log(U_t/U_t_0), deg=1, cov=True)
lower_err_t_U = np.sqrt(np.diag(lower_cov_t_U))

lower_val_t_U = unp.uarray(lower_par_t_U, lower_err_t_U)

with open('build/t-U_RC_lower.tex', 'w') as file:
	file.write(r'\qty{')
	file.write(f'{(-1000/lower_val_t_U[0]).n:.1f}({(-1000/lower_val_t_U[0]).s:.1f})')
	file.write(r'}{\micro\second}')

# Zeitabhängige Entladungskurve
#################################################################################################################
# Frequenzabhängiger Amplitudenverlauf

f_U, U_f = np.genfromtxt('data/f-U_data.txt', unpack=True) 
U_f_0 = U_f[0]

def A(f,u):
	return np.sqrt(1/(1 + (u*2*np.pi*f)**2))

par_f_U, cov_f_U = curve_fit(A, f_U, U_f/U_f_0)
err_f_U = np.sqrt(np.diag(cov_f_U))

val_f_U = unp.uarray(par_f_U, err_f_U)

x = np.linspace(9,11000,10000)

plt.figure(figsize=(5.78,3.25))
ax = plt.subplot(1,1,1)
ax.yaxis.set_minor_locator(AutoMinorLocator(4))
f, = plt.plot(x, A(x, *par_f_U), '-', c='olivedrab')
p, = plt.plot(f_U, U_f/U_f_0, 'kx', ms=3.21)
plt.xlabel(r'$\nu \mathbin{/} \unit{\hertz}$')
plt.ylabel(r'$U\hspace{-0.2ex}/U_{\! 0}$')
plt.xscale('log')
leg = plt.legend([p, f], ['Messdaten', 'Regression'], edgecolor='k', facecolor='none')
leg.get_frame().set_linewidth(0.25)

plt.savefig('build/f-U_plot.pdf')
plt.close()

with open('build/f-U_0.tex', 'w') as file:
	file.write(r'\qty{')
	file.write(f'{U_f_0:.0f}')
	file.write(r'}{\milli\volt}')
with open('build/f-U_RC.tex', 'w') as file:
	file.write(r'\qty{')
	file.write(f'{1000000*np.abs(val_f_U[0].n):.1f}({1000000*val_f_U[0].s:.1f})')
	file.write(r'}{\micro\second}')
with open('build/f-U_RC_corr.tex', 'w') as file:
	file.write(r'\qty{')
	file.write(f'{11/11.6*1000000*np.abs(val_f_U[0].n):.1f}({11/11.6*1000000*val_f_U[0].s:.1f})')
	file.write(r'}{\micro\second}')

# Frequenzabhängiger Amplitudenverlauf
#################################################################################################################
# Frequenzabhängiger Phasenverlauf

f_phi, a, b = np.genfromtxt('data/f-phi_data.txt', unpack=True)
phi_f = a/b
while (np.max(phi_f[:-1]) > 0.25):
	phi_f[:-1][phi_f[:-1] > 0.25] = phi_f[:-1][phi_f[:-1] > 0.25] - 0.25

def phase(f,u):
	return np.arctan(-2*np.pi*f*u)

par_f_phi, cov_f_phi = curve_fit(phase, f_phi, 2*np.pi*phi_f)
err_f_phi = np.sqrt(np.diag(cov_f_phi))

val_f_phi = unp.uarray(par_f_phi, err_f_phi)

x = np.linspace(5,15000,10000)

ax = plt.subplot(1,1,1)
ax.locator_params(nbins=10, axis='y')
f, = plt.plot(x, 180/np.pi*phase(x, *par_f_phi), '-', c='olivedrab')
p, = plt.plot(f_phi, 360*phi_f, 'kx', ms=3.21)
plt.xlabel(r'$\nu \mathbin{/} \unit{\hertz}$')
plt.ylabel(r'$\varphi \hspace{-0.2ex} \mathbin{/} \unit{\degree}$')
plt.xscale('log')
leg = plt.legend([p, f], ['Messdaten', 'Regression'], edgecolor='k', facecolor='none')
leg.get_frame().set_linewidth(0.25)

plt.savefig('build/f-phi_plot.pdf')
plt.close()

with open('build/f-phi_RC.tex', 'w') as file:
	file.write(r'\qty{')
	file.write(f'{1000000*np.abs(val_f_phi[0].n):.1f}({1000000*val_f_phi[0].s:.1f})')
	file.write(r'}{\micro\second}')
with open('build/f-phi_RC_corr.tex', 'w') as file:
	file.write(r'\qty{')
	file.write(f'{11/11.6*1000000*np.abs(val_f_phi[0].n):.1f}({11/11.6*1000000*val_f_phi[0].s:.1f})')
	file.write(r'}{\micro\second}')

# Frequenzabhängiger Phasenverlauf
#################################################################################################################
# Polarplot

phi = np.linspace(0,np.pi/2,10000)

fig = plt.figure()
ax = fig.add_subplot(1,1,1, polar=True)
ax.tick_params(labelleft=True, labelright=False, labeltop=False, labelbottom=True)
ax.grid(c='k', linewidth=0.25)
ax.set_theta_zero_location('E')
ax.set_theta_direction(1)
ax.set_thetamin(0)
ax.set_thetamax(90)
ax.set_rorigin(0)
ax.set_rmin(0)
ax.set_rmax(1)
f, = plt.plot(phi, np.cos(phi), '-', c='olivedrab')
p, = plt.plot(np.take(2*np.pi*phi_f, [0,1,2,3,4,6,8]), np.take(U_f/U_f_0, [9,10,11,12,16,19,24]), 'kx', ms=3.21)
plt.xlabel(r'$U \hspace{-0.2ex}\vphantom{\dfrac{\dfrac{A^1}{}}{}} / U_{\! 0}$')
leg = ax.legend([p,f],['Messdaten','Theoriekurve'],edgecolor='k',facecolor='none',bbox_to_anchor=(1.175,1.1))
leg.get_frame().set_linewidth(0.25)

plt.savefig('build/polar_plot.pdf')
plt.close()

# Polarplot
#################################################################################################################
# Tabelle Entladungskurve

table_header = r'''     \begin{tabular}
		{S[table-format=1.2]
		 S[table-format=1.2]
		 S[table-format=1.2]
		 S[table-format=1.2]
		 @{\hspace{5.55ex}}
		 S[table-format=1.2]
		 S[table-format=1.2]
		 S[table-format=1.2]
		 S[table-format=1.2]}
		\toprule
		{$t \mathbin{/} \unit{\milli\second}$} &
		{$U \hspace{-0.15ex} \mathbin{/} \unit{\volt}$} &
		{$U \hspace{-0.2ex} / U_{\!0}$} &
		{$\ln(U \hspace{-0.2ex} / U_{\!0})$} &
		{$t \mathbin{/} \unit{\milli\second}$} &
		{$U \hspace{-0.15ex} \mathbin{/} \unit{\volt}$} &
		{$U \hspace{-0.2ex} / U_{\!0}$} &
		{$\ln(U \hspace{-0.2ex} / U_{\!0})$} \\
		\midrule
'''
table_footer = r'''		\bottomrule
	\end{tabular}
'''
row_template = r'		{0:1.2f} & {1:1.2f} & {2:1.2f} & {3:1.2f} & {4:1.2f} & {5:1.2f} & {6:1.2f} & {7:1.2f} \\'

with open('build/t-U_table.tex', 'w') as f:
    f.write(table_header)
    for row in zip(t[:9], U_t[:9], U_t[:9]/U_t_0, np.log(U_t[:9]/U_t_0), t[9:], U_t[9:], U_t[9:]/U_t_0, np.log(U_t[9:]/U_t_0)):
        f.write(row_template.format(*row))
        f.write('\n')
    f.write(table_footer)

# Tabelle Entladungskurve
#################################################################################################################
# Tabelle Amplitudenkurve

table_header = r'''     \begin{tabular}
		{S[table-format=3.0]
		 S[table-format=3.0]
		 S[table-format=1.3]
		 @{\hspace{6.66ex}}
		 S[table-format=5.0]
		 S[table-format=3.0]
		 S[table-format=1.3]}
		\toprule
		{$\nu \mathbin{/} \unit{\hertz}$} &
		{$U \mathbin{/} \unit{\milli\volt}$} &
		{$U \hspace{-0.2ex} / U_{\!0}$} &
		{$\nu \mathbin{/} \unit{\hertz}$} &
		{$U \mathbin{/} \unit{\milli\volt}$} &
		{$U \hspace{-0.2ex} / U_{\!0}$} \\
		\midrule
'''
table_footer = r'''		\bottomrule
	\end{tabular}
'''
row_template = r'		{0:3.0f} & {1:3.0f} & {2:1.3f} & {3:5.0f} & {4:3.0f} & {5:1.3f} \\'

with open('build/f-U_table.tex', 'w') as f:
    f.write(table_header)
    for row in zip(f_U[:14], U_f[:14], U_f[:14]/U_f_0, f_U[14:], U_f[14:], U_f[14:]/U_f_0):
        f.write(row_template.format(*row))
        f.write('\n')
    f.write(table_footer)

# Tabelle Amplitudenkurve
#################################################################################################################
# Tabelle Phasenkurve

table_header = r'''     \begin{tabular}
		{S[table-format=5.0]
		 S[table-format=1.3]
		 S[table-format=1.3]
		 S[table-format=1.3]
		 S[table-format=1.3]
		 S[table-format=2.1]
		 S[table-format=1.3]}
		\toprule
		{$\nu \mathbin{/} \unit{\hertz}$} &
		{$\tilde a \mathbin{/} \unit{\milli\second}$} &
		{$\,\tilde{\!b} \mathbin{/} \unit{\milli\second}$} &
		{$\tilde{a} / \,\tilde{\!b}$} & {$a / b$} &
		{$\varphi \hspace{-0.2ex} \mathbin{/} \unit{\degree}$} &
		{$\varphi \hspace{-0.2ex} \mathbin{/} \unit{\radian}$} \\
		\midrule
'''
table_footer = r'''		\bottomrule
	\end{tabular}
'''
row_template = r'		{0:5.0f} & {1:1.3f} & {2:1.3f} & {3:1.3f} & {4:1.3f} & {5:2.1f} & {6:1.3f} \\'

with open('build/f-phi_table.tex', 'w') as f:
    f.write(table_header)
    for row in zip(f_phi, a, b, a/b, phi_f, 360*phi_f, 2*np.pi*phi_f):
        f.write(row_template.format(*row))
        f.write('\n')
    f.write(table_footer)

# Tabelle Phasenkurve
#################################################################################################################
