'''
Plot sims with expliict parameter values. For implicit, see run_scenarios.py.
'''

import numpy as np
import sciris as sc
import covasim as cv
import create_sim as cs
import pandas as pd
print(cs.__file__)

do_save = True
baseline_file = 'results/kc-baseline-nov17.msim'

T = sc.tic()

indices    = np.arange(8)
ncpus      = min(len(indices), 12)
from_cache = 0

# Load data
tsdf = pd.read_excel(cs.epi_data_file)
ctdf = pd.read_excel(cs.contact_tracing_file, sheet_name='final')
ctdf['date'] = pd.to_datetime(ctdf['date']).dt.date
tsdf['date'] = pd.to_datetime(tsdf['date']).dt.date


projpars = sc.objdict(
        trace_time   = {k:1   for k in 'hswcl'},
        i_factor     = {k:0.1 for k in 'hswcl'},
        q_factor     = {k:0.2 for k in 'hswcl'},
        rel_severe   = 1.0, # Change to 0.5 to match hospitalizations, but disrupts optimization
        death_change = 0.3,
        test_delay   = 1,
        trace_probs  = sc.mergedicts({k:0.6 for k in 'hsl'}, {'w':0.05}, {'c': 0.0}),
)

kwargs = dict(projpars=projpars, from_cache=from_cache, do_shrink=0)


if __name__ == '__main__':

    # Settings
    if ncpus>1:
        sims = sc.parallelize(cs.run_sim, iterarg=indices, ncpus=ncpus, kwargs=kwargs)
    else:
        sims = [cs.run_sim(index, **kwargs) for index in indices]

    msim = cv.MultiSim(sims=sims)
    msim.median()
    sim = msim.base_sim
    res = msim.results
    plot_msim = True
    if plot_msim:
        msim.plot(to_plot='overview', fig_args={'figsize':(38,19)})

    if do_save:
        msim.save(baseline_file)

    sc.toc(T)
