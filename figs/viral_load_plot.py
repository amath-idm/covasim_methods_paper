"""
Instantiate a population and infect all on day 0, watch what happens next and make plots.
"""

import pandas as pd
import numpy as np
import covasim as cv
import matplotlib.pyplot as plt
import seaborn as sns
import sciris as sc

sns.set(font_scale=2)
sns.set_style('white')#, {'axis.grid': True})
# plt.rcParams['font.family'] = 'Proxima Nova'

do_save = True
# cv.check_save_version('1.5.1', die=True)


def create_sim(TRACE_PROB=None,  TEST_PROB=None, TRACE_TIME=None, TEST_DELAY=None):

    # initialize simulation
    sim = cv.Sim(version='2.0.3')

    # PARAMETERS
    pars = {}
    pars.update({'pop_size': 5e4}) # start with a small pool << actual population

    # diagnosed individuals maintain same beta
    pars.update({'pop_infected': pars['pop_size']}) # Infect ALL!
    pars.update({'beta': 0}) # No transmission
    pars.update({'n_days': 30}) # 40d is long enough for everything to play out

    # update parameters
    sim.update_pars(pars=pars)

    return sim

# nonideal, but avouds having to figure out how to import utils from covasim
# and allows changes for here.
def compute_viral_load(t,     time_start, time_recovered, time_dead,  frac_time,    load_ratio,    high_cap):
    '''
    Calculate relative transmissibility for time t. Includes time varying
    viral load, pre/asymptomatic factor, diagonsis factor, etc.

    Args:
        t: (int) timestep
        time_start: (float[]) individuals' infectious date
        time_recovered: (float[]) individuals' recovered date
        time_dead: (float[]) individuals' death date
        frac_time: (float) frac of time in high load
        load_ratio: (float) ratio for high to low viral load
        high_cap: (float) cap on the number of days with high viral load

    Returns:
        load (float): viral load
    '''

    # Get the end date from recover or death
    n = len(time_dead)
    time_stop = np.ones(n)*time_recovered # This is needed to make a copy
    inds = ~np.isnan(time_dead)
    time_stop[inds] = time_dead[inds]

    # Calculate which individuals with be high past the cap and when it should happen
    infect_days_total = time_stop-time_start
    trans_day = frac_time*infect_days_total
    inds = trans_day > high_cap
    cap_frac = high_cap/infect_days_total[inds]

    # Get corrected time to switch from high to low
    trans_point = np.ones(n)*frac_time
    trans_point[inds] = cap_frac

    # Calculate load
    load = np.ones(n) # allocate an array of ones with the correct dtype
    early = (t-time_start)/infect_days_total < trans_point # are we in the early or late phase
    load = (load_ratio * early + load * ~early)/(load+frac_time*(load_ratio-load)) # calculate load

    # set to zero if not infectious
    inds = (t-time_start) < 0
    load[inds] = 0
    inds = (t-time_start) >= infect_days_total
    load[inds] = 0

    return load

def get_viral_loads(sim, t_start, t_end):
    viralLoad = np.zeros((len(sim.people.date_exposed),t_end-t_start))
    par1 = sim.pars['beta_dist']['par1']
    par2 = sim.pars['beta_dist']['par2']
    mean  = np.log(par1**2 / np.sqrt(par2 + par1**2)) # Computes the mean of the underlying normal distribution
    sigma = np.sqrt(np.log(par2/par1**2 + 1)) # Computes sigma for the underlying normal distribution
    rel_trans = np.random.lognormal(mean=mean, sigma=sigma, size=len(sim.people.date_exposed))

    for i in np.arange(t_end-t_start):
        viralLoad[:,i] = compute_viral_load(t_start+i, sim.people.date_infectious,\
                 sim.people.date_recovered, sim.people.date_dead,\
                 sim.pars['viral_dist']['frac_time'], sim.pars['viral_dist']['load_ratio'],\
                 sim.pars['viral_dist']['high_cap'])*rel_trans
    return viralLoad

if __name__ == "__main__":

    # directory for storing results
    do_run = True

    fn = 'sim_inds.sim'
    if do_run:
        sim = create_sim()
        sim.run()
        sim.save(fn)
    else:
        sim = cv.Sim.load(fn)


    t_start = 0
    t_end = 35
    interval = t_end - t_start
    viralLoad = get_viral_loads(sim, t_start, t_end)

    cols = ['age', 'date_exposed', 'date_infectious', 'date_symptomatic',
            'date_severe', 'date_critical', 'date_recovered', 'date_dead']
    raw = {c: getattr(sim.people, c) for c in cols}
    dat = pd.DataFrame(raw, columns=cols)
    dat['AgeBin'] = pd.cut(dat['age'], bins=list(range(0,101,10)))
    dat['symptomatic'] = ~dat['date_symptomatic'].isna()
    dat['severe'] = ~dat['date_severe'].isna()
    dat['critical'] = ~dat['date_critical'].isna()
    dat['dead'] = ~dat['date_dead'].isna()

    #Group first person in each age bin
#    age_bins = dat['AgeBin'].unique()
    idx = []
#    for age in age_bins:
#        idx.append((dat['AgeBin']==age).argmax(axis=0))
#    idx = [5, 27, 1, 0, 17, 13, 49995, 7, 18, 19, 91]

    #Group a fixed number of samples from each age bin
#    num_age_bin = 2
#    index = dat.groupby('AgeBin').apply(lambda s: s.sample(min(len(s), num_age_bin),\
#                                        replace=False).index)
#    for i in np.arange(len(index)):
#        for j in np.arange(len(index.values[i])):
#            idx.append(index.values[i][j])

    #Random sample
    np.random.seed(12) # Was 5
    idx = dat.sample(30, replace=False).index
    idx = idx.values
    #Make sure we have at least one death
    #idx_temp = (~np.isnan(dat['date_dead'])).argmax(axis=0)
    # idx_temp = 207
    # idx = np.append(idx, dat.index[idx_temp])

    sort_idx = np.argsort(-dat['age'][idx])
    idx = np.array(idx)[sort_idx]

    # Set the colormap
    try:
        import cmasher as cm
        import matplotlib.colors as colors

        def truncate_colormap(cmap, minval=0.0, maxval=1.0, n=256):
            ''' From https://stackoverflow.com/questions/18926031/how-to-extract-a-subset-of-a-colormap-as-a-new-colormap-in-matplotlib '''
            new_cmap = colors.LinearSegmentedColormap.from_list(
                'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
                cmap(np.linspace(minval, maxval, n)))
            return new_cmap

        cmap = cm.ember.reversed()
        cmap = truncate_colormap(cmap, minval=0.0, maxval=0.5)
    except Exception as E:
        cmap = 'Wistia'
        print(f'Could not load cmasher colormaps ({str(E)}), using {cmap}')

    dy = 0.004
    use_zero_color = False # Optionally use a different color for zero viral load
    dat_reduced = dat.loc[idx,:]
    viralLoad_reduced = viralLoad[idx,:]
    viralLoad_reduced = viralLoad_reduced/np.amax(viralLoad_reduced[:,:])
    vl_norm = np.log10(viralLoad_reduced+0.1)
    vl_colors = sc.arraycolors(vl_norm, cmap=cmap)
    fig = plt.figure(figsize=(12, 14))
    ax = plt.axes([0.10,0.08, 0.9, 0.9])
    for i in np.arange(len(idx)):
        point_size = 140
        line_size = 5
        im = plt.scatter(np.zeros(interval), np.zeros(interval), c=vl_norm[i,:], cmap=cmap, s=0) # For the colorbar
        x = np.arange(interval)
        y = np.ones(interval)*.01*i
        plt.scatter(x, y, c=vl_colors[i,:], s=point_size, edgecolors='k', linewidths=0)
        if use_zero_color:
            zero_color = [[0.8,0.8,0.8]]
            zero_inds = viralLoad_reduced[i,:]==0
            plt.scatter(x[zero_inds], y[zero_inds], c=zero_color, s=point_size, edgecolors='k', linewidths=0)
        t_symp = dat_reduced['date_symptomatic'].array[i]
        if ~np.isnan(t_symp):
            symp_line, = plt.plot([t_symp-.5, t_symp-.5], [0.01*i-dy,0.01*i+dy], linewidth=line_size, color='c')
        t_sev = dat_reduced['date_severe'].array[i]
        if ~np.isnan(t_sev):
            sev_line, = plt.plot([t_sev-.5, t_sev-.5], [0.01*i-dy,0.01*i+dy], linewidth=line_size, color='b')
        t_crit = dat_reduced['date_critical'].array[i]
        if ~np.isnan(t_crit):
            crit_line, = plt.plot([t_crit-.5, t_crit-.5], [0.01*i-dy,0.01*i+dy], linewidth=line_size, color='r')
        t_dead = dat_reduced['date_dead'].array[i]
        if ~np.isnan(t_dead):
            inter_temp = max(0,int(t_end-t_dead))
            if inter_temp == 0:
                x = []
            else:
                x = np.arange(t_dead,t_end)
            dead_line, = plt.plot([t_dead-.5, t_dead-.5],\
                 [0.01*i-dy,0.01*i+dy], linewidth=line_size, color='k')
            plt.scatter(x, np.ones(inter_temp)*.01*i,\
                        c='w', s=point_size)
        t_rec = dat_reduced['date_recovered'].array[i]
        if ~np.isnan(t_rec):
            inter_temp = max(0,int(t_end-t_rec))
            if inter_temp == 0:
                x = []
            else:
                x = np.arange(t_rec,t_end)
            rec_line, = plt.plot([t_rec-.5, t_rec-.5],\
                 [0.01*i-dy,0.01*i+dy], linewidth=line_size, color=[0.0,0.8,0.0])
            plt.scatter(x, np.ones(inter_temp)*.01*i,\
                        c='w', s=point_size+2)
    plt.xlim(t_start-0.5, t_end)
    plt.ylim(-dy,0.01*len(idx)+dy)
    plt.yticks(0.01*np.arange(len(idx)), labels=np.int32(np.array(round(dat_reduced['age']))))
    plt.xlabel('Day post exposure')
    plt.ylabel('Age of individual (years)', labelpad=10.0)
    cbar = plt.colorbar(im,ticks=[])
    plt.text(1.163*t_end,0.135,'Viral load',rotation=270)
    plt.legend((symp_line, sev_line, crit_line, rec_line, dead_line),\
               ('Symptomatic', 'Severe', 'Critical', 'Recovery', 'Death'),\
               #bbox_to_anchor=(1,0.6))
               loc='upper right')
    plt.show()

    if do_save:
        cv.savefig('viral-load-by-age_v3.png', dpi=150)
        print('Plot saved.')