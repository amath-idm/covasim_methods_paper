'''
Generate Fig. 1 (calibration) and Fig. 2 (projections) for the report.
'''

import pylab as pl
import sciris as sc
import covasim as cv
import run_scenarios as rs
from plot_utils import plotter as pu_plotter, make_figure, format_axs



#%% Options
do_save = 1
plot_fig1 = 1
plot_fig2 = 1
use_mean = 1 # Use the mean for the 'best' line to get deaths right, but use median otherwise to get bounds right
pl.rc('figure', dpi=180)
fig_dpi = 200
plot_start = '2020-11-01'
plot_end = '2020-12-31'


#%% Loading

T = sc.tic()
msims = cv.load(rs.scenarios_file)
plot_msims = 0
if plot_msims:
    for key,msim in msims.items():
        print(f'Plotting {key}...')
        fig = pl.figure(key)
        cv.maximize(fig=fig)
        msim.plot(fig=fig, to_plot='overview')
refsim = msims[0].base_sim # Pull out any sim

r = sc.objdict({k:msims[k].base_sim.results for k in msims.keys()})
l = sc.objdict({k:msims[k].label for k in msims.keys()})


#%% Plotting

scargs = dict(color='k', marker='d', linestyle='none', alpha=0.4, markersize=5, lw=1.5)
plargs = dict(lw=2, alpha=1.0)
fillargs = dict(alpha=0.2)

def plotter(ax, res, key, **kwargs):
    return pu_plotter(res, refsim, fillargs=fillargs, plargs=plargs, ax=ax, key=key, **kwargs)


fig, axs = make_figure('Fig. 3: Scenarios', n=3)
ax1, ax2, ax3 = axs

# Panel 1
ax1.set_title('Impact of current testing/tracing programs')
plotter(ax1, r.status_quo, 'new_infections', label=l.status_quo, ylabel='New daily infections')
for k in ['zero_delay', 'no_ttq']:
    plotter(ax1, r[k], 'new_infections', label=l[k])

# Panel 2
ax2.set_title('Impact of increased testing and tracing')
plotter(ax2, r.status_quo, 'new_infections', label=l.status_quo, ylabel='New daily infections')
for k in ['increased_testing', 'increased_tracing', 'increased_ttq']:
    plotter(ax2, r[k], 'new_infections', label=l[k])

# # Panel 3
ax3.set_title('Impact of mobility restrictions')
plotter(ax3, r.status_quo, 'new_infections', label=l.status_quo, ylabel='New daily infections')
for k in ['mob_restrict', 'mob_plus_ttq', 'stay_home']:
    plotter(ax3, r[k], 'new_infections', label=l[k])

# Tidying
format_axs(refsim, axs)
pl.show()
if do_save:
    cv.savefig('results/fig3_scenarios.png', dpi=fig_dpi)


print('Done.')
sc.toc(T)