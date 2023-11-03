"""
Show contacts within different layers. Used with PR #482 to check that network
properties are preserved.
"""

import covasim as cv
import networkx as nx
import pylab as pl
import numpy as np
import sciris as sc
import synthpops as sp

import matplotlib as mplt
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, LinearSegmentedColormap
from matplotlib.ticker import LogLocator, LogFormatter
import matplotlib.font_manager as font_manager
from mpl_toolkits.axes_grid1 import make_axes_locatable

import os
import cmocean
import cmasher as cmr
import seaborn as sns

# NB -- must be on synthpops commit 311e8f9

if not sp.config.full_data_available:
    pytest.skip("Data not available, tests not possible", allow_module_level=True)

try:
    username = os.path.split(os.path.expanduser('~'))[-1]
    fontdirdict = {
        'dmistry': '/home/dmistry/Dropbox (IDM)/GoogleFonts',
        'cliffk': '/home/cliffk/idm/covid-19/GoogleFonts',
    }
    if username not in fontdirdict:
        fontdirdict[username] = os.path.expanduser(os.path.expanduser('~'), 'Dropbox', 'GoogleFonts')

    font_path = fontdirdict[username]

    fontpath = fontdirdict[username]
    font_style = 'Roboto_Condensed'
    fontstyle_path = os.path.join(fontpath, font_style, font_style.replace('_', '') + '-Light.ttf')
    prop = font_manager.FontProperties(fname=fontstyle_path)
    mplt.rcParams['font.family'] = prop.get_name()
except:
    mplt.rcParams['font.family'] = 'Roboto'



pop_size = [118, 124, 125, 126, 127, 128, 131][3]
pop_type = 'synthpops'
undirected = True

np.random.seed(1)

contacts = dict(
    random = {'a':20},
    hybrid = {'h': 2.0, 's': 4, 'w': 6, 'c': 10},
    synthpops = {'h': 2.0, 's': 4, 'w': 6, 'c': 10},
    )

pars = {
    'pop_size': pop_size, # start with a small pool
    'pop_type': pop_type, # synthpops, hybrid
    'contacts': contacts[pop_type],
}

# Create sim
sim = cv.Sim(pars=pars)
sim.initialize()

fig = plt.figure(figsize=(7,10), dpi=130)
axis_args = dict(left=0.10, bottom=0.05, right=0.70, top=0.95, wspace=0.25, hspace=0.15)
plt.subplots_adjust(**axis_args)
mplt.rcParams['font.family'] = 'Roboto Condensed'
mplt.rcParams['font.size'] = 16

mapping = dict(a='All', h='Household networks', 
    s='School network', 
    w='Workplace networks', c='Community')
discrete = False

# cmap = 'nipy_spectral' # 'gist_rainbow' # 'parula''

cmap1 = plt.get_cmap('cmr.ember')
cmap2 = plt.get_cmap('cmr.lavender_r')
new_cmap_name = 'ember_lavender'

colors1 = cmap1(np.linspace(0.4, 1, 96))  # to truncate darker end of cmap1 change 0 to a value greater than 0, less than 1
colors2 = cmap2(np.linspace(0., 1, 128))  # to truncate darker end of cmap2 change 1 to a value less than 1, greater than 0

# transition_steps = 0  # heat+freeze
transition_steps = 4 # increase if closest ends of the color maps are far apart, values to try: 4, 8, 16
transition = mplt.colors.LinearSegmentedColormap.from_list("transition", [cmap1(1.), cmap2(0)])(np.linspace(0,1,transition_steps))
colors = np.vstack((colors1, transition, colors2))
colors = np.flipud(colors)

new_cmap = mplt.colors.LinearSegmentedColormap.from_list(new_cmap_name, colors)
cmap = new_cmap

# cmap = plt.get_cmap('cmr.pride')
# cmap = mplt.colors.LinearSegmentedColormap.from_list("new_cmap", cmap(np.linspace(0, 0.85, 128)))


age_cutoffs = np.arange(0,101,10) # np.array([0, 4, 6, 18, 22, 30, 45, 65, 80, 90, 100])
if discrete:
    raw_colors = sc.vectocolor(len(age_cutoffs), cmap=cmap)
    colors = []
    for age in sim.people.age:
        ind = sc.findinds(age_cutoffs<=age)[-1]
        colors.append(raw_colors[ind])
    colors = np.array(colors)
else:
    age_map = sim.people.age*0.1+np.sqrt(sim.people.age)
    colors = sc.vectocolor(age_map, cmap=cmap)


# Create the legend
ax = pl.axes([0.85, 0.05, 0.14, 0.93])
ax.axis('off')
for age in age_cutoffs:
    nearest_age = sc.findnearest(sim.people.age, age)
    col = colors[nearest_age]
    plt.plot(np.nan, np.nan, 'o', c=col, label=f'Age {age}')
plt.legend()


keys = list(contacts[pop_type].keys())
keys = keys[:3] # Ignore community
# keys = ['h', 'w'] # Ignore schools

# Find indices
idict = {}
hdfs = {}
for layer in keys:
    hdf = sim.people.contacts[layer].to_df()
    hdfs[layer] = hdf
    idict[layer] = list(set(list(hdf['p1'].unique()) + list(hdf['p2'].unique())))
    if layer == 'h':
        orig_h = idict[layer]
        idict[layer] = list(range(pop_size))

trimmed_s = sc.dcp(idict['s'])
for ind in idict['h']:
    if ind in trimmed_s and ind not in orig_h:
        trimmed_s.remove(ind)
trimmed_h = sc.dcp(idict['h'])
for ind in trimmed_s:
    if ind in trimmed_h:
        trimmed_h.remove(ind)
for ind in idict['w']:
    if ind in trimmed_h:
        trimmed_h.remove(ind)

ndict = dict(h=20, s=50, w=50)
kdict = dict(h=0.7, s=1.0, w=2.0)
mdict = dict(h='^', s='o', w='s')
use_graphviz = True
keys = ['h', 'w']
for i, layer in enumerate(keys):
    ax = pl.subplot(2,1,i+1)
    hdf = hdfs[layer]
    inds = idict[layer]
    color = colors[inds]

    G = nx.DiGraph()
    G.add_nodes_from(inds)
    p1 = hdf['p1']
    p2 = hdf['p2']
    G.add_edges_from(zip(p1,p2))
    if undirected:
        G.add_edges_from(zip(p2,p1))
    print(f'Layer: {layer}, nodes: {G.number_of_nodes()}, edges: {G.number_of_edges()}')

    if use_graphviz:
        pos = nx.nx_pydot.graphviz_layout(G)
    else:
        pos = nx.spring_layout(G, k=kdict[layer], iterations=200)
    nx.draw(G, pos=pos, ax=ax, node_size=ndict[layer], node_shape=mdict[layer], width=0.1, arrows=False, node_color=color)
    if layer=='h': # Warning, assumes that everyone is in a household
        for sublayer in 'hsw':
            sli = idict[sublayer]
            if sublayer == 's':
                sli = trimmed_s
            elif sublayer == 'h':
                sli = trimmed_h
            subG = G = nx.DiGraph()
            subG.add_nodes_from(sli)
            subpos = {i:pos[i] for i in sli}
            nx.draw(subG, pos=subpos, ax=ax, node_size=50, node_shape=mdict[sublayer], width=0.1, arrows=False, node_color=colors[sli])
    ax.set_title(mapping[layer])
# plt.show()
plt.savefig('network_layout_v3.png', dpi=200)
