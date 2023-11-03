'''
Create figure for benchmarking the model
'''

import numpy   as np
import pylab   as pl
import sciris  as sc
import covasim as cv
import matplotlib as mpl

start_at_1000 = True
if start_at_1000:
    pop_sizes = np.array([1e3, 2e3, 5e3, 10e3, 20e3, 50e3, 100e3, 200e3, 500e3, 1000e3])
else:
    pop_sizes = np.array([100, 200, 500, 1e3, 2e3, 5e3, 10e3, 20e3, 50e3, 100e3, 200e3, 500e3, 1000e3])
n_sizes = len(pop_sizes)

pars = dict(
    pop_type='random', # Was 'hybrid',
    n_days=100,
)

#%% Test speed

rerun_speed = False

if rerun_speed:

    results = {k:[] for k in pop_sizes}
    T0 = sc.tic()
    repeats = 3
    for p in pop_sizes:
        for r in range(repeats):
            sim = cv.Sim(pars, pop_size=p, verbose=0)
            sim.initialize()
            T = sc.tic()
            sim.run()
            elapsed = sc.toc(T, output=True)
            results[p].append(elapsed)
            print(f'Ran {p}: {elapsed:0.2f} s ({sc.toc(T0, output=True):0.2f} s); infections: {sim.summary["cum_infections"]}')

    sc.toc(T0)

else:
    results = {
        # Original
        # 100: [0.27550673484802246, 0.2355809211730957, 0.2365407943725586, 0.2450428009033203, 0.2367246150970459],
        # 200: [0.24162030220031738, 0.24155545234680176, 0.23961853981018066, 0.24763059616088867, 0.3022189140319824],
        # 500: [0.2786753177642822, 0.2647714614868164, 0.2662160396575928, 0.2647979259490967, 0.2593207359313965],
        # 1000.0: [0.2829773426055908, 0.2863593101501465, 0.27786850929260254, 0.2793593406677246, 0.33785247802734375],
        # 2000.0: [0.3230257034301758, 0.32685089111328125, 0.3400561809539795, 0.3364284038543701, 0.39693570137023926],
        # 5000.0: [0.46991491317749023, 0.5052773952484131, 0.4511244297027588, 0.5047757625579834, 0.44927477836608887],
        # 10000.0: [0.7071657180786133, 0.7393951416015625, 0.7829132080078125, 0.7315793037414551, 0.7270290851593018],
        # 20000.0: [1.173414707183838, 1.114853858947754, 1.1554431915283203, 1.1675190925598145, 1.1691868305206299],
        # 50000.0: [2.5121631622314453, 2.579982280731201, 2.5996506214141846, 2.589463472366333, 2.787954807281494],
        # 100000.0: [5.422731876373291, 5.251848459243774, 5.226514101028442, 5.56469988822937, 5.3730902671813965],
        # 200000.0: [10.182217359542847, 10.160524368286133, 9.880542516708374, 10.04503583908081, 10.12425684928894],
        # 500000.0: [24.50388264656067, 25.800947666168213, 25.745570182800293, 25.43501901626587, 26.076338291168213],
        # 1000000.0: [52.80573606491089, 51.96183943748474, 53.75685095787048]

         # New, hybrid
         #    1000.0: [0.23, 0.24, 0.22],
         #    2000.0: [0.25, 0.27, 0.27],
         #    5000.0: [0.30, 0.31, 0.31],
         #   10000.0: [0.39, 0.39, 0.38],
         #   20000.0: [0.55, 0.63, 0.57],
         #   50000.0: [1.11, 1.25, 1.09],
         #  100000.0: [2.05, 2.22, 2.11],
         #  200000.0: [4.05, 4.34, 4.07],
         #  500000.0: [10.3, 10.6, 10.3],
         # 1000000.0: [18.8, 20.6, 19.0],

         # New, random
            1000.0: [0.08309531211853027, 0.07750940322875977, 0.08118724822998047],
            2000.0: [0.09429717063903809, 0.09399294853210449, 0.09648609161376953],
            5000.0: [0.12806320190429688, 0.1304476261138916, 0.134260892868042],
           10000.0: [0.17401885986328125, 0.17401432991027832, 0.17141389846801758],
           20000.0: [0.27066969871520996, 0.2872636318206787, 0.2722644805908203],
           50000.0: [0.591052770614624, 0.6034688949584961, 0.5902981758117676],
          100000.0: [1.1512539386749268, 1.1550524234771729, 1.1539506912231445],
          200000.0: [2.2536728382110596, 2.2686684131622314, 2.2643091678619385],
          500000.0: [5.822768211364746, 5.899953603744507, 6.142377853393555],
         1000000.0: [14.346067667007446, 13.602745294570923, 14.191432237625122],
    }

# Analyze results

r = sc.objdict()
for k in 'blh':
    r[k] = np.zeros(n_sizes)

for i,p in enumerate(pop_sizes):
    r.b[i] = np.median(np.array(results[p]))
    r.l[i] = np.array(results[p]).min()
    r.h[i] = np.array(results[p]).max()


#%% Test memory -- this is slow so just run manually

rerun_memory = False

if rerun_memory:
    print('Please copy memory usage manually into mresults dict')
    sim = cv.Sim(pars, pop_size=100e3, verbose=0)
    sc.mprofile(sim.run)

else:
    mresults = {
        100: 0.6,
        200: 0.9,
        500: 1.6,
        1000: 2.3,
        2000: 3.0,
        5000: 7.4,
        10000: 10.7,
        20000: 15.4,
        50000: 41.5,
        100000: 101.5,
        200000: 199.7,
        500000: 509.5,
        1000000: 855.9,
        }

m_arr = np.array([mresults[p] for p in pop_sizes])


#%% Test days and infections

rerun_days = False

days = np.arange(0,101,10)
days[0] = 1 # Can't run with 0 days

if rerun_days:

    dpars = dict(
        pop_type='hybrid',
        pop_size=100e3,
    )


    dresults = {k:0 for k in days}
    iresults = {k:0 for k in days}

    T0 = sc.tic()
    for d in days:
        sim = cv.Sim(dpars, n_days=d, verbose=0)
        sim.initialize()
        T = sc.tic()
        sim.run(verbose=0)
        elapsed = sc.toc(T, output=True)
        ires = sim.summary["cum_infections"]
        dresults[d] = elapsed
        iresults[d] = ires
        print(f'Ran {d}: {elapsed:0.2f} s ({sc.toc(T0, output=True):0.2f} s); infections: {ires}')

    sc.toc(T0)

else:
    dresults = {
         1: 0.04146766662597656,
         10: 0.21635770797729492,
         20: 0.4205513000488281,
         30: 0.6406049728393555,
         40: 0.8705778121948242,
         50: 1.0619094371795654,
         60: 1.4124529361724854,
         70: 1.7442219257354736,
         80: 2.076596260070801,
         90: 2.5131449699401855,
         100: 2.9370434284210205
     }
    iresults = {
        1: 10.0,
        10: 35.0,
        20: 134.0,
        30: 445.0,
        40: 1382.0,
        50: 4190.0,
        60: 12088.0,
        70: 29523.0,
        80: 53293.0,
        90: 70305.0,
        100: 77794.0
    }
d_arr = np.array([dresults[d] for d in days])
i_arr = np.array([iresults[d] for d in days])

#%% Plot results

do_plot = True
err_bars = False

if do_plot:

    fig = pl.figure(figsize=(12,9))
    axis_args = dict(left=0.10, bottom=0.10, right=0.95, top=0.97, wspace=0.25, hspace=0.25)
    pl.subplots_adjust(**axis_args)
    pl.rcParams['font.family'] = 'Proxima Nova'
    pl.rcParams['font.size'] = 18
    msize = 12

    ax = pl.subplot(2,1,1)
    if err_bars:
        errcol = [0.5, 0.5, 0.5]  # Error bars are too small to see!
        for i,p in enumerate(pop_sizes):
            pl.loglog([p, p], [r.l[i], r.h[i]], c=errcol)
    pl.plot(pop_sizes, r.b, 'o', markersize=msize, label='Simulation', zorder=100)
    pl.plot(pop_sizes, pop_sizes/70e3, c='k', label='Constant time (7 million person-days / second)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_ylim([0.01, 100])
    sc.commaticks(axis='x')
    ax.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    pl.xlabel('Population size')
    pl.ylabel('CPU time per 100 days (s)')
    pl.legend()

    ax = pl.subplot(2,1,2)
    pl.plot(pop_sizes, m_arr, 'o', markersize=msize, label='Simulation', zorder=100)
    pl.plot(pop_sizes, pop_sizes/1e3, c='k', label='Constant memory (1000 people / MB)')
    ax.set_xscale('log')
    ax.set_yscale('log')
    # ax.set_ylim([1, 1000])
    sc.commaticks(axis='x')
    sc.commaticks(axis='y')
    # ax.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    pl.xlabel('Population size')
    pl.ylabel('Memory usage (MB)')
    pl.legend()

    pl.show()
    cv.savefig('covasim-performance_v2.png', dpi=150)

    # ax = pl.subplot(2,2,2)
    # pl.plot(days, d_arr, 'o', markersize=msize, label='Simulation')
    # # pl.plot(pop_sizes, pop_sizes/1e3, c='k', label='Constant memory (1 KB / person)')
    # # ax.set_ylim([1, 1000])
    # # sc.commaticks(axis='x')
    # # sc.commaticks(axis='y')
    # # ax.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    # # pl.xlabel('Population size')
    # # pl.ylabel('Memory usage (MB)')
    # pl.legend()

    # ax = pl.subplot(2,2,4)
    # pl.plot(i_arr, d_arr, 'o', markersize=msize, label='Simulation')
    # ax.set_xscale('log')
    # ax.set_yscale('log')
    # # pl.plot(pop_sizes, pop_sizes/1e3, c='k', label='Constant memory (1 KB / person)')
    # # ax.set_ylim([1, 1000])
    # # sc.commaticks(axis='x')
    # # sc.commaticks(axis='y')
    # # ax.get_yaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    # # pl.xlabel('Population size')
    # # pl.ylabel('Memory usage (MB)')
    # pl.legend()


