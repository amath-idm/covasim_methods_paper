'''
Create an idealized network plot.
'''

#%% Preliminaries

import numpy as np
import sciris as sc
import pylab as pl
import covasim as cv
import svg_data as svg

# Options
do_save = 1
do_show = 1
cv.options.set(dpi=100, font_family='Roboto')
sc.tic()

# Get actual Covasim network -- good seeds: 3-senegal, 1-malawi, 6-malawi, 8-malawi, 9-malawi
n_people = 20
pars = sc.objdict(
    rand_seed = 9,
    pop_size  = n_people,
    pop_type  = 'hybrid',
    location  = 'malawi'
)
sim = cv.Sim(pars)
sim.initialize()
if do_save:
    sim.save('network.sim', keep_people=True)

# Find clusters in the network
c = sc.objdict() # c = clusters
lkeys = sim.layer_keys()
for lk in lkeys:
    contacts = sim.people.contacts[lk]
    inds = np.unique(sc.cat(contacts['p1'], contacts['p2'])).tolist()
    c[lk] = sc.odict()
    while len(inds):
        root = inds.pop(0)
        rk = str(root)
        c[lk][rk] = [root]
        contacts = sim.people.contacts[lk].find_contacts(root).tolist()
        contacts = [i for i in contacts if c not in c[lk][rk]] # Sometimes get self connections
        c[lk][rk].extend(contacts)
        for i in contacts:
            if i in inds:
                inds.remove(i)
    for rk in c[lk]:
        c[lk][rk] = np.unique(c[lk][rk]).tolist()

# Remove workplace duplicates
diff = list(set(c.w[0]) - set(c.w[1]))
c.w[1] = c.w[0] # Swap order
c.w[0] = diff # Reset
c.w[1] = list(set(c.w[1]) - set(c.w[0]))

# Perform checks
print('Network:')
print(c)
n_h = len(c.h)
n_s = 1
n_w = 2
assert len(c.s) == n_s and len(c.w) == n_w, 'Sorry, hard-coded to certain numbers of schools/workplaces!'


#%% Define helper functions

def inds2circ(n, x0=0, y0=0, r=1, aspect=1, rot=0):
    ''' Convert n indices into a circle of coordinates with a specified radius, aspect ratio, and rotation '''
    x = np.zeros(n)
    y = np.zeros(n)
    radians = 2*np.pi*np.arange(n)/n+rot
    x = np.cos(radians)
    y = np.sin(radians)
    x = x*r + x0
    y = y*r/aspect + y0
    return x,y


def plargs(lkey):
    ''' Get plot arguments '''
    args = sc.objdict(c=lcolors[lkey], lw=lw[lkey], alpha=0.15, zorder=-10)
    if lkey == 'c': # Too many lines otherwise
        args.alpha *= 0.5
    return args


def i2k(i):
    ''' Convert an index into a sortable string key '''
    return f'p{i:02d}'


def label(x, y, s):
    ''' Plot a label '''
    x = sc.promotetoarray(x).mean()
    return ax.text(x, y, s, fontsize=16, horizontalalignment='center', bbox=dict(facecolor='w', alpha=0.4, lw=0))


def plotsvg(x=50, y=50, marker=None, size=50, fill='none', color='k', lw=1.0, alpha=1.0):
    ''' Add an SVG image to the plot '''
    return pl.scatter([x], [y], marker=marker, s=size*1000, facecolor=fill, edgecolors=color, linewidths=lw, alpha=alpha)


#%% Set up figure
figx = 16
figy = 16
fig = pl.figure(figsize=(figx,figy))
ax0 = pl.axes([0, 0.54, 1, 0.45]) # Upper panel
ax = pl.axes([0, 0, 1, 0.5]) # Lower panel
ax.set_xlim([3,100]) # Nominally [0,100], but trim a little
ax.set_ylim([10,100])
lcolors = sc.objdict(h='#377eb8', s='#e41a1c', w='#4daf4a', c='#a24e99') # From ttq-app.covasim.org
lw = sc.objdict(sim['beta_layer']) # Line widths
lw[:] *= 3 # Make thicker
ax0.axis('off')
ax.axis('off')
fig.text(0.02, 0.97, 'A', fontsize=40)
fig.text(0.02, 0.47, 'B', fontsize=40)


#%% Top panel -- load and plot image
schematic_fn = 'contact_network_mar29_edit.png'
imdata = pl.imread(schematic_fn)
ax0.imshow(imdata)


#%% Bottom panel

# Plot households
hy0 = 50
xmin = 10
xmax = 90
hx0 = np.linspace(xmin, xmax, n_h)
r = 5 # Radius
δ = 5 # Standard distance unit
rots = [0.5, 0.5, 0.8, 0.6, 0.3] # Rotation of each household circle (in radians)
ppl = sc.odict()
aspect = 1.0 # Adjust default aspect ratio of each household circle
age_labels = False # Whether to plot age labels above each person
for hou,house in c.h.enumvals():
    ages = sim.people.age[house]
    markers = [svg.person for h in house] # Or ['$♂️$' if sim.people.sex[h] else '$♀️$' for h in house]
    hx,hy = inds2circ(len(house), x0=hx0[hou], y0=hy0, r=r, aspect=aspect, rot=rots[hou])
    ihouse = np.arange(len(house))
    label(hx, hy0+2.8*δ, f'Household {hou+1}')
    for i in ihouse:
        idk = i2k(house[i])
        ppl[idk] = sc.objdict()
        ppl[idk].hx = hx[i]
        ppl[idk].hy = hy[i]
    for i in ihouse:
        x = hx[i]
        y = hy[i]+0.06*ages[i] # So lines appear at people's feet
        s = 40*(min(ages[i],30)+7) # Make older people larger
        pl.scatter(x, y, s=s, marker=markers[i], c='k', alpha=0.5) # Or c=sc.gridcolors(n_people)[house[i]]])
        if age_labels:
            pl.text(x-0.4, y+0.001*s+1.3, int(np.round(ages[i])))
    for i1 in ihouse:
        for i2 in ihouse:
            pl.plot([hx[i1], hx[i2]], [hy[i1], hy[i2]], **plargs('h'))

# Find people not belonging to a household cluster (single person households) and remove
missing = []
for i in np.arange(n_people):
    key = i2k(i)
    if key not in ppl:
        ppl[key] = None
        missing.append(i) # If no household contacts
ppl.sort()


# Plot community
cy0 = 85
cx0 = np.linspace(xmin+2.6*δ, xmax-2.8*δ, n_people)
pl.fill_between([xmin+2.0*δ, xmax-1.4*δ], cy0-2*δ, cy0+2.0*δ, alpha=0.3, color=lcolors.c)
plotsvg(y=cy0, marker=svg.community, size=550)
label(cx0, cy0+2.4*δ, 'Community')
for comm in c.c.values():
    icomm = np.arange(len(comm))
    plotted = []
    for i1 in icomm:
        for i2 in icomm:
            p1 = comm[i1]
            p2 = comm[i2]
            p1p2 = (p1,p2)
            if p1p2 not in plotted: # Don't plot the same connection more than once
                plotted.append(p1p2)
                if p1 not in missing:
                    pl.plot([ppl[p1].hx, cx0[p2]], [ppl[p1].hy, cy0], **plargs('c'))


# Plot school
sy0 = 15
sxmin = 45
sxmax = 55
nsch = len(c.s[0])
sx0 = np.linspace(sxmin, sxmax, nsch)
plotsvg(y=sy0+2*δ, marker=svg.school, size=30, fill=lcolors.s, alpha=0.5, lw=1.3)
label(sx0, sy0+3.1*δ, 'School')
for sch in c.s.values():
    isch = np.arange(nsch)
    for i1 in isch:
        for i2 in isch:
            p1 = sch[i1]
            p2 = sch[i2]
            if p1 not in missing:
                pl.plot([ppl[p1].hx, sx0[i2]], [ppl[p1].hy, sy0], **plargs('s'))


# Plot workplaces
wy0 = 17
wxmins = [15, 75]
wxmaxs = [25, 85]
plotsvg(x=19, y=wy0+1*δ, marker=svg.workplace1, size=30, fill=lcolors.w, alpha=0.5, lw=1.3)
plotsvg(x=79, y=wy0+1*δ, marker=svg.workplace2, size=30, fill=lcolors.w, alpha=0.5, lw=1.3)
for wi,work in c.w.enumvals():
    wxmin = wxmins[wi]
    wxmax = wxmaxs[wi]
    nwork = len(work)
    iwork = np.arange(nwork)
    wx0 = np.linspace(wxmin, wxmax, nwork)
    label(wx0, wy0+3.3*δ, f'Workplace {wi+1}')
    for i1 in iwork:
        for i2 in iwork:
            p1 = work[i1]
            p2 = work[i2]
            if p1 not in missing:
                pl.plot([ppl[p1].hx, wx0[i2]], [ppl[p1].hy, wy0], **plargs('w'))


#%% Tidy up
if do_save:
    cv.savefig('network_schematic.png')

if do_show:
    pl.show()

sc.toc()
print('Done.')