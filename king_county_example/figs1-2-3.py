'''
Generate Fig. 1 (calibration) and Fig. 2 (projections) for the report.
'''

import pandas as pd
import pylab as pl
import sciris as sc
import covasim as cv
import create_sim as cs
import run_baseline as rb
import run_scenarios as rs
from plot_utils import plotter as pu_plotter, make_figure, format_axs, roll, sub
print(cs.__file__)

#%% Options
do_save = 1
plot_fig1 = 1
plot_fig2 = 1
use_mean = 1 # Use the mean for the 'best' line to get deaths right, but use median otherwise to get bounds right
pl.rc('figure', dpi=180)
fig_dpi = 200
calib_end = '2020-11-15'
plot_end = '2020-12-31'

# Figure configuration
# font_size = 22
# font_family = 'Proxima Nova'
# pl.rcParams['font.size'] = font_size
# pl.rcParams['font.family'] = font_family


#%% Loading

T = sc.tic()
msim = cv.load(rb.baseline_file)
msims = cv.load(rs.scenarios_file)
refsim = msim.base_sim
res = msim.results
calib_end = refsim.day(calib_end)
if use_mean:
    m2 = cv.load(rb.baseline_file) # Reload and compute mean
    m2.mean()
    for k,v in m2.results.items():
        try:
            msim.results[k].values = m2.results[k].values
        except:
            pass
    msims2 =  cv.load(rs.scenarios_file) # Reload and compute mean
    for mk,mv2 in msims2.items():
        mv2.mean()
        for k,v in mv2.results.items():
            try:
                msims[mk].results[k].values = mv2.results[k].values
            except:
                pass


r = sc.objdict({k:msims[k].base_sim.results for k in msims.keys()})
l = sc.objdict({k:msims[k].label for k in msims.keys()})
import run_scenarios as rs

# Load data
tsdf = pd.read_excel(cs.epi_data_file)
ctdf = pd.read_excel(cs.contact_tracing_file, sheet_name='final')
ctdf['date'] = pd.to_datetime(ctdf['date']).dt.date
tsdf['date'] = pd.to_datetime(tsdf['date']).dt.date
ctdf = ctdf.dropna() # These went from None to NaT at some point


#%% Plotting

scargs = dict(color='k', marker='d', linestyle='none', alpha=0.4, markersize=5, lw=1.5)
plargs = dict(lw=2, alpha=1.0)
fillargs = dict(alpha=0.2)

def plotter(*args, **kwargs):
    return pu_plotter(res, refsim, fillargs, plargs, *args, **kwargs)

def plotter2(ax, res, key, **kwargs):
    return pu_plotter(res, refsim, fillargs=fillargs, plargs=plargs, ax=ax, key=key, **kwargs)

#%% Plotting

fig, axs = make_figure('Fig. 11: Illustrative example', n=5, figsize=(9,12))
ax1, ax2, ax3, ax4, ax5 = axs

# Panel 1
ax1.plot(refsim.day(sub(tsdf['date'].values)), sub(roll(tsdf['data_tests'])), **scargs, label='Data')
plotter(ax1, 'new_tests', ylabel='Tests')

# Panel 2
ax2.plot(refsim.day(sub(tsdf['date'].values)), sub(roll(tsdf['data_diagnoses'])), **scargs, label='Data')
plotter2(ax2, r.mob_restrict, 'new_diagnoses', ylabel='Diagnoses')

# Panel 3
ax3.plot(refsim.day(ctdf['date'].values), ctdf['estimated_daily'], **scargs, label='Data')
plotter2(ax3, r.mob_restrict, 'new_quarantined', ylabel='Contacts traced')

# Panel 4
ax4.plot(refsim.day(sub(tsdf['date'].values)), sub(roll(tsdf['data_deaths'])), **scargs, label='Data')
plotter2(ax4, r.mob_restrict, 'new_deaths', ylabel='Deaths')
ax4.set_ylim([0,20])

# Panel 5
plotter2(ax5, r.mob_restrict, 'new_infections', label=l.mob_restrict)
plotter2(ax5, r.status_quo, 'new_infections', label='No additional countermeasures', ylabel='New infections')
plotter2(ax5, r.mob_plus_ttq, 'new_infections', label=l.mob_plus_ttq, ylabel='New infections')


# Plot vertical lines
for ax in axs: #[ax4]:
    ax.set_xlim([0, refsim.day(plot_end)])
    ax.axvline(refsim.day(calib_end), ymax=10, c='k', lw=1, alpha=0.5, linestyle='--')

# Tidying
format_axs(refsim, axs)

# Reorder last legend
h,l = ax5.get_legend_handles_labels()
o = [1,0,2]
ax5.legend([h[i] for i in o], [l[i] for i in o], frameon=False)

pl.show()
if do_save:
    cv.savefig('results/fig11_illustrative_example_v2.png', dpi=fig_dpi)



print('Done.')
sc.toc(T)