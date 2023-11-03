'''
Animation of transmission tree plotting for no interventions, testing only,
and TTQ.
'''

import covasim as cv
import sciris as sc
import numpy as np
from ete3 import Tree, TreeStyle, NodeStyle


tstart = sc.tic()
cv.git_info('transtree_example.info')

plot_sim    = 0
verbose     = 0
animate     = 0
animate_all = 0

isday = 5
ieday = None
pop_type = 'random'
lkeys = dict(
    random='a',
    hybrid='hswc'
    )
tp = cv.test_prob(start_day=isday, end_day=ieday, symp_prob=0.10, asymp_prob=0.0, symp_quar_prob=1.0, asymp_quar_prob=1.0, test_delay=0)
ct = cv.contact_tracing(trace_probs={k:0.5 for k in lkeys[pop_type]},
                          trace_time={k:5 for k in lkeys[pop_type]},
                          start_day=isday)

contacts = dict(
    random=10,
    hybrid=dict(h=2, s=2, w=2, c=2),
    )
beta = dict(
    random=0.5/contacts['random'],
    hybrid=0.1,
    )
pop_infected = dict(
    random=1,
    hybrid=10,
    )
pars = dict(
    pop_size = 300,
    pop_infected = pop_infected[pop_type],
    pop_type = pop_type,
    n_days = 90,
    contacts = contacts[pop_type],
    beta = beta[pop_type],
    rand_seed = 3248,
    )

labels = ['No interventions', 'Testing only', 'Test + trace']
sims = sc.objdict()
sims.base  = cv.Sim(pars) # Baseline
sims.test  = cv.Sim(pars, interventions=tp) # Testing only
sims.trace = cv.Sim(pars, interventions=[tp, ct]) # Testing + contact tracing

tts = sc.objdict()
for key,sim in sims.items():
    sim.run()
    sim.people.make_detailed_transtree()
    tts[key] = sim.people.transtree.detailed
    if plot_sim:
        to_plot = cv.get_sim_plots()
        to_plot['Total counts']  = ['cum_infections', 'cum_diagnoses', 'cum_quarantined', 'n_quarantined']
        sim.plot(to_plot=to_plot)

""" Swap in ete3 tree drawing in the plotting loop below """

keys = ['base', 'test', 'trace']
for key in [keys[0]]:# sims.keys():

    tt = tts[key]

    tt.sort(key=lambda x: x['date'] if x is not None else 0)  # so source is before target in transmission pairs

    t = Tree()
    nodes = dict()

    for i, entry in enumerate(tt):

        if entry:
            source = entry['source']
            target = entry['target']
            target_date = entry['date']
            if source:
                nodes[target] = nodes[source].add_child(name=target)  # transmitted case
                nodes[target].add_features(date_exposed=target_date)

                # set branch length to time between source and target infection date
                nodes[target].dist = target_date - nodes[source].date_exposed

            else:
                nodes[target] = t.add_child(name=target)  # index case
                nodes[target].add_features(date_exposed=-0.5, dist=5.5)

            date_t = entry.t.date_tested
            date_d = entry.t.date_diagnosed
            date_q = entry.t.date_known_contact
            debug_s = 't=%d : (%s->%d) :' % (target_date, source, target)

            if not np.isnan(date_t):
                nodes[target].add_features(date_tested=date_t)
                debug_s += ' date_tested=%d' % date_t
            if not np.isnan(date_d):
                nodes[target].add_features(date_diagnosed=date_d)
                debug_s += ' date_diagnosed=%d' % date_d
            if not np.isnan(date_q):
                nodes[target].add_features(date_quarantined=date_q)
                debug_s += ' date_known_contact=%d' % date_q
            print(debug_s)

    ts = TreeStyle()
    ts.show_leaf_name = False
    ts.scale = [4,36][0] # Use 4 for columns layout, 36 for rows layout

    for n in t.traverse():

        nstyle = NodeStyle()
        nstyle['size'] = 20
        nstyle['fgcolor'] = 'FireBrick' # '#b22222'
        n.set_style(nstyle)

        if 'date_quarantined' in getattr(n, 'features'):
            n.img_style['fgcolor'] = 'CornflowerBlue'
        elif 'date_diagnosed' in getattr(n, 'features'):
            n.img_style['fgcolor'] = 'LightGreen'
        elif 'date_tested' in getattr(n, 'features'):
            n.img_style['fgcolor'] = 'DarkOrange' # '#ff8c00'

    t.show(tree_style=ts)

sc.toc(tstart)
