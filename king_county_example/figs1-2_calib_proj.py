'''
Generate Fig. 1 (calibration) and Fig. 2 (projections) for the report.
'''

import numpy as np
import pandas as pd
import pylab as pl
import sciris as sc
import covasim as cv
import create_sim as cs
import run_baseline as rb
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

# Figure configuration
# font_size = 22
# font_family = 'Proxima Nova'
# pl.rcParams['font.size'] = font_size
# pl.rcParams['font.family'] = font_family


#%% Loading

T = sc.tic()
msim = cv.load(rb.baseline_file)
plot_msim = 0
if plot_msim:
    msim.plot(to_plot='overview')
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

# Load data
tsdf = pd.read_excel(cs.epi_data_file)
ctdf = pd.read_excel(cs.contact_tracing_file, sheet_name='final')
ctdf['date'] = pd.to_datetime(ctdf['date']).dt.date
tsdf['date'] = pd.to_datetime(tsdf['date']).dt.date


#%% Plotting

scargs = dict(color='k', marker='d', linestyle='none', alpha=0.4, markersize=5, lw=1.5)
plargs = dict(lw=2, alpha=1.0)
fillargs = dict(alpha=0.2)

def plotter(*args, **kwargs):
    return pu_plotter(res, refsim, fillargs, plargs, *args, **kwargs)


#%% Fig. 1

if plot_fig1:

    fig, axs = make_figure('Fig. 1: Calibration')
    ax1, ax2, ax3, ax4 = axs

    # Panel 1
    ax1.set_title('Number of daily tests')
    ax1.plot(refsim.day(sub(tsdf['date'].values[:calib_end])), sub(roll(tsdf['new_tests'][:calib_end])), **scargs, label='Data')
    plotter(ax1, 'new_tests', ylabel='Tests', end=calib_end)

    # Panel 2
    ax2.set_title('Number of daily diagnoses and contacts traced')
    ax2.plot(refsim.day(sub(tsdf['date'].values)), sub(roll(tsdf['new_diagnoses'])), **scargs, label='Data (diagnoses)')
    ax2.plot(refsim.day(ctdf['date'].values), ctdf['estimated_daily'], **sc.mergedicts(scargs, dict(color=[1,0,0])), label='Data (contacts traced)')
    plotter(ax2, 'new_diagnoses', label='Model (diagnoses)', end=calib_end)
    plotter(ax2, 'new_quarantined', label='Model (contacts traced)', ylabel='Diagnoses/people traced', end=calib_end)

    # Panel 3
    ax3.set_title('Test positivity rate (%)')
    ax3.plot(refsim.day(sub(tsdf['date'].values)), sub(roll(tsdf['test_yield']*100)), **scargs, label='Data')
    plotter(ax3, 'test_yield', label='Model', pct=True, ylabel='Yield (%)', end=calib_end)
    ax3.set_ylim([0,20])

    # Panel 4
    ax4.set_title('Number of daily deaths')
    ax4.plot(refsim.day(sub(tsdf['date'].values)), sub(roll(tsdf['new_deaths'])), **scargs, label='Data')
    plotter(ax4, 'new_deaths', label='Model', ylabel='Deaths', end=calib_end)
    ax4.set_ylim([0,20])

    # Tidying
    format_axs(refsim, axs)
    pl.show()
    if do_save:
        cv.savefig('results/fig1_calibration.png', dpi=fig_dpi)



#%% Fig. 2

if plot_fig2:

    fig, axs = make_figure('Fig. 2: Projections')
    ax1, ax2, ax3, ax4 = axs

    # Panel 1
    ax1.set_title('Number of cumulative infections')
    plotter(ax1, 'cum_infections', ylabel='Cumulative infections')

    # Panel 2
    ax2.set_title('Number of daily new infections')
    plotter(ax2, 'new_infections', ylabel='Infections')

    # Panel 3
    ax3.set_title('Effective reproduction number')
    plotter(ax3, 'r_eff', ylabel=r'$R_{eff}$')

    # Panel 4
    ax4.set_title('Case detection rate (%)')
    factor = 100
    date = refsim.day(res['date'])
    new_diag = res['new_diagnoses'].values*factor # Use only best estimates
    best = new_diag/np.maximum(res['new_infections'].values, 1)
    low  = new_diag/np.maximum(res['new_infections'].high, 1) # Low = high since denominator
    high = new_diag/np.maximum(res['new_infections'].low, 1) # Avoid divide by 0
    ax4.fill_between(date, roll(low), roll(high), **fillargs)
    ax4.plot(date, roll(best), **plargs, label='Model')
    ax4.set_ylabel('Case detection rate (%)')

    # Tidying
    format_axs(refsim, axs)
    ax3.set_ylim([0.5, 2.5])
    ax3.axhline(1, c='k', lw=1, linestyle=':')
    ax3.legend(frameon=False, loc='upper right')
    pl.show()
    if do_save:
        cv.savefig('results/fig2_projections.png', dpi=fig_dpi)


print('Done.')
sc.toc(T)