import covasim as cv
import sciris as sc
import pylab as pl


tstart = sc.tic()

plot_sim    = 0
verbose     = 0
animate     = 0
animate_all = 0

isday = 15
ieday = None
pop_type = 'random'
lkeys = dict(
    random='a',
    hybrid='hswc'
    )
tp = cv.test_prob(start_day=isday, end_day=ieday, symp_prob=0.10, asymp_prob=0.0, symp_quar_prob=1.0, asymp_quar_prob=1.0, test_delay=0)
ct = cv.contact_tracing(trace_probs={k:0.8 for k in lkeys[pop_type]},
                          trace_time={k:3 for k in lkeys[pop_type]},
                          start_day=isday)

pars = dict(
    pop_size = 20000,
    pop_infected = 200,
    pop_type = 'random',
    n_days = 90,
    contacts = 10,
    beta = 0.05,
    rand_seed = 3248,
    )

# Scenario metaparameters
metapars = dict(
    n_runs    = 10, # Number of parallel runs; change to 3 for quick, 11 for real
)

# Define the scenarios
scenarios = {'baseline': {
              'name':'No interventions',
              'pars': {
                  'interventions': None,
                  }
              },
            'testing': {
              'name':'Testing only',
              'pars': {
                  'interventions': tp
                  }
              },
            'ttq': {
              'name':'Testing + tracing',
              'pars': {
                  'interventions': [tp, ct]
                  }
              },
             }

# Run the scenarios
scens = cv.Scenarios(basepars=pars, metapars=metapars, scenarios=scenarios)
scens.run()
scens.plot()
pl.savefig('scenarios-example.png')