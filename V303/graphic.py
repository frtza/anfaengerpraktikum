import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator

plt.figure(figsize = (3.65, 5.01875))

gs = gridspec.GridSpec(4, 1)
gs.update(hspace = 0)

def format_func(value, tick_number):
    N = int(np.round(2 * value / np.pi))
    if N == 0:
        return "0"
    elif N == 1:
        return r"$\pi/2$"
    elif N == 2:
        return r"$\vphantom{0}\pi$"
    elif N % 2 > 0:
        return r"${0}\pi/2$".format(N)
    else:
        return r"${0}\pi$".format(N // 2)

def set_style(ax):
	plt.axis('on')
	plt.tick_params(
		which='both',
		direction='in',
		bottom=True,
		top=True,
		left=True,
		right=True,
		labelbottom=False,
		labeltop=False,
		labelleft=True,
		labelright=False)
	plt.xlim(0, 6*np.pi)
	plt.ylim(-1.5, 1.5)
	plt.axhline(y=0, color='k', linestyle='--', dashes=(4.25, 2.25), linewidth=0.75)
	ax.xaxis.set_minor_locator(AutoMinorLocator(2))
	ax.yaxis.set_minor_locator(AutoMinorLocator(2))
	ax.xaxis.set_major_locator(plt.MultipleLocator(np.pi))
	ax.xaxis.set_minor_locator(plt.MultipleLocator(np.pi/2))
	ax.xaxis.set_major_formatter(plt.FuncFormatter(format_func))
	x_ticks = ax.xaxis.get_ticklines()[:]
	x_ticks[2].set_visible(False)
	x_ticks[3].set_visible(False)
	x_ticks[14].set_visible(False)
	x_ticks[15].set_visible(False)
	y_ticks = ax.yaxis.get_ticklines()[:]
	y_ticks[4].set_visible(False)
	y_ticks[5].set_visible(False)
	leg = plt.legend(frameon=False, loc='center', labelcolor='linecolor', handlelength=0, handletextpad=0,
		bbox_to_anchor=(0.075,0.35))
	leg.legendHandles[0].set_visible(False)

ax0 = plt.subplot(gs[0])
plt.plot(np.linspace(0, 6*np.pi, 100000), np.sin(np.linspace(0, 6*np.pi, 100000)),
	linewidth=1.125, color='steelblue', label=r'$U_{\! \text{sig}}$')
ax0.set_yticklabels(['', r'$-U_{\! 0}$', '0', r'$U_{\! 0}$'])
set_style(ax0)

ax1 = plt.subplot(gs[1])
plt.plot(np.pi * np.array([0,1,1,2,2,3,3,4,4,5,5,6]), [1,1,-1,-1,1,1,-1,-1,1,1,-1,-1],
	linewidth=1.125, color='salmon', label=r'$U_{\! \text{ref}}$')
set_style(ax1)

ax2 = plt.subplot(gs[2])
plt.plot(np.linspace(0, 6*np.pi, 100000), np.abs(np.sin(np.linspace(0, 6*np.pi, 100000))),
	linewidth=1.125, color='#edc824', label=r'$\, U_{\! \text{mix}}$')
ax2.set_yticklabels(['', r'$-U_{\! 0}$', '0', r'$U_{\! 0}$'])
set_style(ax2)

ax3 = plt.subplot(gs[3])
plt.plot([0, 6*np.pi], [0.75, 0.75], linewidth=1.125, color='olivedrab', label=r'$U_{\! \text{out}}$')
ax3.set_yticklabels(['', r'$-U_{\! 0}$', '0', r'$U_{\! 0}$'])
set_style(ax3)
plt.tick_params(labelbottom=True)

plt.savefig('build/graphic.pdf', bbox_inches='tight')
