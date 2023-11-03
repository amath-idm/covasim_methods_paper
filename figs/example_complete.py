import covasim as cv

# Define the simulation parameters
pars = dict(
    pop_size     = 200e3,        # Define a population size of 200,000 people
    pop_infected = 75,           # Start with 75 infected individuals
    beta         = 0.012,        # Calibrate overall transmission to this setting
    pop_type     = 'hybrid',     # Use realistic household, school, and work contacts
    start_day    = '2020-02-10', # First day of the simulation
    end_day      = '2020-06-29', # Last day of the simulation
)

# Define the interventions
tr_probs = dict(h=0.9, s=0.7, w=0.7, c=0.3) # Probability that a contact in each layer will be traced
tr_time  = dict(h=0.0, s=1.0, w=1.0, c=3.0) # Time required to trace contacts in each layer (days)
m_days    = ['2020-03-26', '2020-04-10', '2020-05-05'] # Days on which mobility changes
m_changes = [0.7, 0.4, 0.8] # Define the mobility changes on these days
interventions = [
    cv.clip_edges(days='2020-03-26', changes=0.0, layers='s'),        # Close schools
    cv.clip_edges(days=m_days, changes=m_changes, layers=['w', 'c']), # Close/reopen work + community
    cv.test_prob(start_day='2020-05-20', symp_prob=0.1, symp_quar_prob=0.8, test_delay=1), # Testing
    cv.contact_tracing(start_day='2020-05-20', trace_probs=tr_probs, trace_time=tr_time)   # Contact tracing
]

# Create and run the simulation
sim = cv.Sim(pars=pars, interventions=interventions)
sim.run()
sim.plot()


# ----
# Save
cv.savefig('example_fig08_interventions.png')
