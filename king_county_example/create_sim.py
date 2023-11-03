'''
Like create_sim.py, but use new data to project into the future.
'''

# Standard packages
import os
import numpy as np
import pandas as pd
import sciris as sc
import covasim as cv


cv.check_save_version('2.0.3', die=False)

# Define the input files -- difference
epi_data_file  = './data/KC_covasim_2020nov17.xlsx'
safegraph_file = './data/KC_weeklyinteractions_2020nov13.xlsx'
contact_tracing_file = './data/contact_tracing.xlsx' # Not used but for imports
popfile_stem   = './inputs/kc_rnr_seed'
jsonfile       = './data/opt_merged_sep20_sg1.json' # Calibrated runs -- difference
json = sc.loadjson(jsonfile) # Difference

# Define default values for calculating the scenarios -- >> difference
d_calcs = sc.objdict()
d_calcs.routine_test     = 0.08 # Baseline daily symptomatic testing
d_calcs.testing_factor   = 2.5 # Scale-up compared to current testing
d_calcs.symp_asymp_ratio = 80 # How much less asymptomatic vs symptomatic testing there is
d_calcs.max_test         = 0.5 # Maximum (routine) testing probability to consider
d_calcs.max_delay        = 7 # Maximum testing or tracing delays to consider
d_calcs.quar_test        = 0.90 # Baseline quarantine testing
d_calcs.iso_factor       = 0.1 # How much people isolate
d_calcs.iq_ratio         = 2 # How much less people quarantine than isolate
d_calcs.tr_prob          = 0.66 # Tracing probability
d_calcs.t_delay          = 1 # Time for test
d_calcs.tr_delay         = 2 # Time for trace
# d_calcs.reopen           = 1.0 # Amount by which to reopen to

# Define actual scenario parameters based on d_calcs values
d_pars = sc.objdict()
d_pars.symp_prob        = d_calcs.routine_test*d_calcs.testing_factor # 0.2 by default
d_pars.asymp_prob       = d_calcs.routine_test*d_calcs.testing_factor/d_calcs.symp_asymp_ratio # 0.0025 by default
d_pars.symp_quar_prob   = d_calcs.quar_test
d_pars.asymp_quar_prob  = d_calcs.quar_test
d_pars.test_delay       = d_calcs.t_delay
d_pars.trace_probs      = sc.mergedicts({k:d_calcs.tr_prob for k in 'hswl'}, {'c':0})
d_pars.trace_time       = {k:d_calcs.tr_delay for k in 'hswcl'}
d_pars.i_factor         = {k:d_calcs.iso_factor for k in 'hswcl'} # Isolation
d_pars.q_factor         = {k:min(1, d_calcs.iso_factor*d_calcs.iq_ratio) for k in 'hswcl'} # Quarantine
# d_pars.reopen           = d_calcs.reopen
# << difference

# >> difference
def load_pars(index): # Only used by Fig. 4
    entry = json[index]
    pars = sc.objdict(entry['pars'])
    pars.rand_seed = int(entry['index'])
    print(f'Loading parameters from trial {index}, mismatch {entry["mismatch"]}...')
    sc.pp(pars)
    return pars
# << difference


# Generate the population filename
def get_popfile(pars):
    n_popfiles = 5
    popfile = popfile_stem + str(pars['rand_seed']%n_popfiles) + '.ppl'

    # Check that the population file exists
    if not os.path.exists(popfile):
        errormsg = f'WARNING: could not find population file {popfile}! Please regenerate first'
        raise FileNotFoundError(errormsg)

    return popfile


def check_contacts(sim, check=False, verbose=True):
    if check:
        contacts = {}
        for lkey in ['h','w','s','c']:
            contacts[lkey] = len(sim.people.contacts[lkey])
        if not hasattr(sim, 'n_contacts'):
            sim.n_contacts = sc.odict()
        sim.n_contacts[sim.date(sim.t)] = contacts
        if verbose:
            print(f'>>> On day {sim.t}, there were {contacts} contacts')
    return


def make_safegraph(sim):
    ''' Create interventions representing SafeGraph data '''

    # Load data.values
    fn = safegraph_file
    df = pd.read_excel(fn)
    week = df['week']
    s = df['p.tot.schools'].values
    w = df['merged'].values # Difference
    c = sc.dcp(w) # Not different enough to warrant different values

    # Do processing
    days = sim.day([pd.Timestamp(v) for v in week.values]) # loads in nanoseconds otherwise lol
    last_day = days[-1]+1
    i_days = np.arange(days[0], last_day)
    s = np.interp(i_days, days, s)
    w = np.interp(i_days, days, w)
    c = np.interp(i_days, days, c)
    days = i_days

    # Create interventions
    interventions = [
        cv.clip_edges(days=days, changes=s, layers='s', label='clip_s'),
        cv.clip_edges(days=days, changes=w, layers='w', label='clip_w'),
        cv.clip_edges(days=days, changes=c, layers='c', label='clip_c'),
        ]

    return interventions


# >> Difference <<

def remove_ltcf_community(sim, debug=False):
    ''' Ensure LTCF residents don't have community contacts '''
    over_65 = sc.findinds(sim.people.age>65)
    llayer = sim.people.contacts['l']
    clayer = sim.people.contacts['c']
    in_ltcf = np.union1d(llayer['p1'], llayer['p2'])
    over_65_in_ltcf = np.intersect1d(over_65, in_ltcf)
    p1inds = sc.findinds(np.isin(clayer['p1'], over_65_in_ltcf))
    p2inds = sc.findinds(np.isin(clayer['p2'], over_65_in_ltcf))
    popinds = np.union1d(p1inds, p2inds)
    if debug:
        print(f'Over 65: {len(over_65)}')
        print(f'Over 65 LTCF: {len(over_65_in_ltcf)}')
        print(f'Pop inds: {len(popinds)}')
        print(f'Orig len: {len(clayer)}')
    clayer.pop_inds(popinds)
    if debug:
        print(f'New len: {len(clayer)}')
    return clayer


def test_num_subtarg(sim, sev=100.0, u20=0.5, quar=1.0):
    ''' Subtarget severe people with more testing, and young people with less '''
    inds = np.arange(len(sim.people))
    vals = np.ones(len(sim.people))
    vals[sc.findinds(sim.people.age<20 * ~sim.people.severe)] *= u20 # People who are under 20 and severe test as if they're severe; * is element-wise "and"
    vals[sim.people.true('severe')] *= sev
    vals[sim.people.true('quarantined')] *= quar
    return {'inds':inds, 'vals':vals}


def create_sim(eind, projpars, verbose=True): # Difference
    ''' Create a single simulation for further use '''

    p = load_pars(eind) # Difference

    # Basic parameters and sim creation
    pars = {  'pop_size'      : 225e3,
              'pop_scale'     : 10,
              'pop_type'      : 'synthpops',
              'pop_infected'  : 300,
              'beta'          : p.beta,
              'start_day'     : '2020-01-27',
              'end_day'       : '2020-12-31', # Difference
              'rescale'       : True,
              'rescale_factor': 1.1,
              'verbose'       : p.get('verbose', 0.1*verbose),
              'rand_seed'     : int(p.rand_seed),
              'analyzers'     : cv.daily_stats(days=[]),
              'beta_layer'    : dict(h=3.0, s=0.6, w=0.6, c=0.3, l=1.5),
              'rel_severe_prob' : projpars['rel_severe'],
            }

    # Create and initialize the sim
    if pars['verbose']:
        print(f'Creating sim! seed={p.rand_seed}') # Difference
    sim = cv.Sim(pars, label=f'base_sim_{p.rand_seed}', popfile=get_popfile(pars), load_pop=True, datafile=epi_data_file) # Create this here so can be used for test numbers etc.-- difference

    # >> Difference
    sim.sceninfo = sc.objdict()
    sim.sceninfo.calib_end  = '2020-05-31'
    # << Difference

    # Define testing interventions -- 97% sensitivity from https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7177629/
    test_kwargs = dict(daily_tests=sim.data['new_tests'], test_delay=2, sensitivity=0.97, subtarget=test_num_subtarg)
    tn = cv.test_num(symp_test=p.tn, start_day='2020-01-27', end_day=None, **test_kwargs, label='tn') # Difference
    interventions = [tn]

    # Define beta interventions
    hwc_days  = ['2020-02-24', '2020-03-23'] # Change date here, 04-27 or 05-04
    hwc_days = sim.day(hwc_days)
    b_wc_ch   = [1.0, p.bc_wc1] # Simple interpoation -- difference
    b_h_ch    = [1.0, 1.1] # Optional household

    all_b_days = np.arange(hwc_days[0], hwc_days[-1]+1) # Full time series
    all_ch_wc = np.interp(all_b_days, hwc_days, b_wc_ch) # Linearly interpolate
    all_ch_h = np.interp(all_b_days, hwc_days, b_h_ch) # Linearly interpolate
    interventions += [cv.change_beta(days=all_b_days, changes=all_ch_h, layers='h', label='beta_h')]
    lkeys = ['w','c','s'] # Difference
    for lkey in lkeys: # Assume schools also use masks, etc. so have same non-movement beta change
        cb = cv.change_beta(days=all_b_days, changes=all_ch_wc, layers=lkey, label=f'beta_{lkey}')
        interventions += [cb]

    # LTCF beta change
    b_l_days = ['2020-02-24', '2020-03-23']
    b_l_days = np.arange(sim.day(b_l_days[0]), sim.day(b_l_days[1]))
    b_l_ch   = np.linspace(1.0, p.bc_lf, len(b_l_days))
    interventions += [cv.change_beta(days=b_l_days, changes=b_l_ch, layers='l', label='beta_l')]
    sim.people.contacts['c'] = remove_ltcf_community(sim) # Remove community contacts from LTCF

    # SafeGraph intervention & tidy up -- difference
    interventions += make_safegraph(sim)
    sim['interventions'] = interventions

    # Don't show interventions in plots, there are too many
    for interv in sim['interventions']:
        interv.do_plot = False

    # These are copied from parameters.py -- needed to get 60-65 group right
    sim['prognoses']['age_cutoffs'] = np.array([0,      10,      20,      30,      40,      50,      65,      70,      80,      90]) # Age cutoffs (upper limits)

    return sim


def modify_sim_proj(sim, projpars, label=None, runinfo=None):
    ''' Modify the simulation for the projections '''

    print(f'  Note: modifying projection {label} at day={sim.t}, date={sim.date(sim.t)}, projpars:\n{projpars}')

    sim.sceninfo.scen_start = '2020-06-01'
    first_scen_day = sim.day(sim.sceninfo.scen_start)
    sim['rel_death_prob'] *= projpars.death_change

    # Change iso_factor and quar_factor
    sim.pars['iso_factor']  = sc.dcp(projpars['i_factor'])
    sim.pars['quar_factor'] = sc.dcp(projpars['q_factor'])

    # Implement testing & tracing interventions
    ctpars = {k:projpars[k] for k in ['trace_probs', 'trace_time']}
    tn = sim.get_intervention('tn')
    tn.test_delay = projpars.test_delay
    ct = cv.contact_tracing(start_day=first_scen_day, label='contact_tracing', **ctpars)
    ct.initialize(sim)
    sim['interventions'] += [ct]

    # Final tidying
    sim.label = label
    sim.runinfo = runinfo
    sim.sceninfo.projpars = projpars

    return sim


def modify_sim_scen(sim, scenpars, label=None, runinfo=None, verbose=False):
    ''' Modify the simulation for the scenarios '''

    print(f'  Note: modifying scenario {label} {scenpars.which} at day={sim.t}, date={sim.date(sim.t)}')

    tn = sim.get_intervention(cv.test_num)
    ct = sim.get_intervention(cv.contact_tracing)
    ce_w = sim.get_intervention('clip_w')
    ce_c = sim.get_intervention('clip_c')
    scen_start = sim.day(scenpars.scen_start)

    def change_mobility(sim, amount):
        for ce in [ce_w, ce_c]: # Clip edges interventions
            valid_days = sc.findinds(ce.days<scen_start)
            ce.days = np.append(ce.days[valid_days], scen_start) # NB, repr of intervention will be wrong with direct modification!
            ce.changes = np.append(ce.changes[valid_days], amount)
        return

    scens_applied = []

    if scenpars.which == 'status_quo':
        scens_applied.append('No change')

    if scenpars.which == 'no_ttq':
        scens_applied.append('No TTQ')
        tn.daily_tests[scen_start:] = 0
        ct.trace_probs = {k:0.0 for k in 'hswcl'}

    if scenpars.which == 'zero_delay':
        scens_applied.append('Zero delay')
        tn.test_delay = 0
        ct.trace_time = {k:0.0 for k in 'hswcl'}

    if scenpars.which in ['increased_testing', 'increased_ttq', 'mob_plus_ttq']:
        scens_applied.append('Increased testing')
        factor = 1.5 # Factor increase in testing rates
        tests = np.array(tn.daily_tests[scen_start:], dtype=float) # Convert to float and back
        tests *= factor
        tests = np.array(tests, dtype=int)
        tn.daily_tests[scen_start:] = tests

    if scenpars.which in ['increased_tracing', 'increased_ttq', 'mob_plus_ttq']:
        scens_applied.append('Increased tracing')
        hsl = 0.90 # Probability of home tracing
        w   = 0.50 # Probability of work tracing
        c   = 0.0 # Probability of community tracing
        ct.trace_probs = sc.mergedicts({k:hsl for k in 'hsl'}, {'w':w, 'c':c})

    if scenpars.which == 'stay_home':
        scens_applied.append('Stay-at-home')
        amount = 0.50 # Minimum value during stay home, stay healthy, corrected from 0.33 to 0.50 for reduced compliance
        change_mobility(sim, amount)

    if scenpars.which in ['mob_restrict', 'mob_plus_ttq']:
        scens_applied.append('Mobility restrictions')
        amount = 0.85*0.85 # Assume restaurants are half open, so half of a 30% reduction (=85%) of current level (also 85%)
        change_mobility(sim, amount)

    if not len(scens_applied):
        raise Exception(f'Scenario "{scenpars.which}" not understood')

    # Final tidying
    sim['verbose'] = 0.1*verbose # Update here since loaded from file
    sim.label = label + ' -- ' + scenpars.which
    sim.sceninfo.scenpars = scenpars
    sim.sceninfo.scens_applied = scens_applied

    return sim


def load_from_cache(filename, index, from_cache, verbose=False):
    sim = None
    sim_loaded = 0
    if from_cache and os.path.exists(filename):
        tries = 0
        while not sim_loaded and tries<3:
            try:
                print(f'run_sim(): Loading sim from cache ({filename})...')
                sim = cv.load(filename)
                assert isinstance(sim, cv.Sim) # Sometimes unpickling fails
                tn = sim.get_intervention('tn')
                # if tn.subtarget == sc.sc_fileio.Failed:
                tn.subtarget = test_num_subtarg
                if verbose:
                    print('FYI, replacing subtargeting function -- all good')
                sim_loaded = 1
            except Exception as E:
                print(f'Loading failed on try {tries}! {E}')
                string = '\n\n'.join([f'Index {index}', filename, f'Try {tries}', str(E), sc.traceback()])
                fn = f'simcache_exception_index{index}_try{tries}.err'
                sc.savetext(fn, string)
                tries += 1
                sc.timedsleep(0.5)
    return sim, sim_loaded


def run_sim(index, projpars=None, label=None, runinfo=None, do_shrink=True, from_cache=True, verbose=True):
    ''' Load, modify, and run the simulation '''

    if projpars is None:
        print('WARNING, loading default parameters!')
        projpars = sc.dcp(d_pars)

    if label is None:
        label = f'full_sim_trial{index}_nolabel'

    # Run sim or load up until scenarios start
    filename = f'./simcache/projection_trial{index}.sim'
    sim, sim_loaded = load_from_cache(filename, index, from_cache)

    if not sim_loaded:
        print(f'run_sim(): Creating new sim (and saving to {filename})...')
        sim = create_sim(index, projpars=projpars, verbose=verbose)
        sim.run(until=sim.sceninfo.calib_end)
        sim.save(filename, keep_people=True)

    # Modify the sim with these scenario values
    sim = modify_sim_proj(sim, projpars=projpars, label=label, runinfo=runinfo)

    # Rerun till the end and optionally shrink
    sim.run(restore_pars=False) # So we can see the changes made
    if do_shrink:
        sim.shrink()

    return sim


def run_scen(index, projpars=None, label=None, runinfo=None, do_shrink=True, from_cache=True, verbose=False, scenpars=None):
    ''' Load, modify, and run the simulation '''

    if projpars is None:
        print('WARNING, loading default parameters!')
        projpars = sc.dcp(d_pars)

    if label is None:
        label = f'scenario_trial{index}'

    # Run sim or load up until scenarios start
    filename = f'./simcache/scenario_trial{index}.sim'
    sim, sim_loaded = load_from_cache(filename, index, from_cache)

    if not sim_loaded:
        print(f'run_scen(): Creating new sim (and saving to {filename})...')
        sim = create_sim(index, projpars=projpars, verbose=verbose)
        sim.run(until=sim.sceninfo.calib_end)
        sim = modify_sim_proj(sim, projpars=projpars, label=label, runinfo=runinfo)
        sim.run(until=scenpars.proj_end) # So we can see the changes made
        sim.save(filename, keep_people=True)

    # Modify the sim with these scenario values
    sim = modify_sim_scen(sim, scenpars=scenpars, label=label, runinfo=runinfo, verbose=verbose)

    # Rerun till the end and optionally shrink
    sim.run(restore_pars=False) # So we can see the changes made
    if do_shrink:
        sim.shrink()


    return sim


if __name__ == '__main__': # Differences -- see also explore_sims.py

    raise NotImplementedError('See explore_sim script')
