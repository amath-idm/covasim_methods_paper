'''
Plot sims with expliict parameter values. For implicit, see run_scenarios.py.

Restrictions: https://medium.com/wagovernor/inslee-announces-statewide-restrictions-for-four-weeks-c0b7da87d34e
'''

import numpy as np
import sciris as sc
import covasim as cv
import create_sim as cs
import pandas as pd
print(cs.__file__)

do_save = True
scenarios_file = 'results/kc-scenarios-nov17.msims'


T = sc.tic()

indices    = np.arange(8)
ncpus      = min(len(indices), 12)
from_cache = 1

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

scenpars = sc.objdict(
    which = None, # Populated below
    proj_end = '2020-11-14',
    scen_start = '2020-11-15',
)

kwargs = sc.objdict(projpars=projpars, scenpars=scenpars, from_cache=from_cache, do_shrink=0, verbose=0)


if __name__ == '__main__':

    whichlist = sc.objdict(

        # Fig. 3A
        status_quo = 'Status quo',
        no_ttq     = 'No testing/tracing',
        zero_delay = 'Zero-delay testing + tracing',

        # Fig. 3B
        increased_testing = 'Increased testing',
        increased_tracing = 'Increased tracing',
        increased_ttq    = 'Increasted testing + tracing',

        # Fig. 3C
        stay_home     = '"Stay Home, Stay Healthy" restrictions (50% reduction)',
        mob_restrict = 'Nov. 16 restrictions (est. 15% reduction)',
        mob_plus_ttq = 'Nov. 16 restrictions and increased testing + tracing',
    )

    msims = sc.odict()
    for i,which,label in whichlist.enumitems():
        sc.heading(f'Running scenario {i+1}/{len(whichlist)}, {which} -- {label}')
        kwargs.scenpars.which = which
        if ncpus>1:
            sims = sc.parallelize(cs.run_scen, iterarg=indices, ncpus=ncpus, kwargs=kwargs)
        else:
            sims = [cs.run_scen(index, **kwargs) for index in indices]
        msim = cv.MultiSim(sims=sims)
        msim.label = label
        msim.median()
        msim.shrink()
        msims[which] = msim

    if do_save:
        cv.save(scenarios_file, msims)

    print('Done.')
    sc.toc(T)