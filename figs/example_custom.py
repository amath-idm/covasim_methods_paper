import covasim as cv

def protect_elderly(sim):
    if sim.t == sim.day('2021-04-01'):
        elderly = sim.people.age>70
        sim.people.rel_sus[elderly] = 0.0

pars = {'pop_size':100e3, 'start_day':'2021-03-01', 'n_days':120}
s1 = cv.Sim(pars, label='Default')
s2 = cv.Sim(pars, label='Protect the elderly', interventions=protect_elderly)
cv.MultiSim([s1, s2]).run().plot(to_plot=['cum_deaths', 'cum_infections'])


# ----
# Save
cv.savefig('example_fig09_custom.png')
