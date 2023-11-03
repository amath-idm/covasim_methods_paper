'''
Make Fig. 8 and Fig. 9
'''

import pylab as pl
import sciris as sc
import covasim as cv

fig8    = 0
fig9    = 1
do_save = 1

# Load image data
fns = sc.objdict()
fns.fig8a = 'example_fig08_interventions.png'
fns.fig8b = 'example_listing_8b.png'
fns.fig9a = 'example_fig09_custom.png'
fns.fig9b = 'example_listing_9b.png'

d = sc.objdict()
for k,fn in fns.items():
    print(f'Loading {k}...')
    d[k] = pl.imread(fn)


#%% Set up figure 8
if fig8:
    figx = 12.5
    figy = 18
    fig = pl.figure(figsize=(figx,figy))
    ax8a = pl.axes([0, 0.44, 1, 0.55]) # Upper panel
    ax8b = pl.axes([0.03, 0.03, 1, 0.41]) # Lower panel
    ax8a.axis('off')
    ax8b.axis('off')
    fig.text(0.02, 0.95, 'A', fontsize=30)
    fig.text(0.02, 0.42, 'B', fontsize=30)

    # Show images
    ax8a.imshow(d.fig8a)
    ax8b.imshow(d.fig8b)

    fn = 'combined_fig08.png'
    if do_save:
        cv.savefig(fn)
    sc.runcommand(f'trim {fn}')
    pl.show()


#%% Set up figure 9
if fig9:
    figx = 12.5
    figy = 13
    fig = pl.figure(figsize=(figx,figy))
    ax8a = pl.axes([0, 0.25, 1, 0.75]) # Upper panel
    ax8b = pl.axes([-0.04, 0.005, 1, 0.25]) # Lower panel
    ax8a.axis('off')
    ax8b.axis('off')
    fig.text(0.05, 0.95, 'A', fontsize=30)
    fig.text(0.05, 0.23, 'B', fontsize=30)

    # Show images
    ax8a.imshow(d.fig9a)
    ax8b.imshow(d.fig9b)

    fn = 'combined_fig09.png'
    if do_save:
        cv.savefig(fn)
    sc.runcommand(f'trim {fn}')
    pl.show()


print('Done.')