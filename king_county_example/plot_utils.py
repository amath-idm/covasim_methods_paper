'''
Shared utilities for plotting
'''

import pylab as pl
import pandas as pd
import sciris as sc
import numpy as np
import datetime as dt
import matplotlib.ticker as ticker

#%% Functions

def make_figure(fignum=None, n=4, figsize=(9, 10)):
    fig = pl.figure(fignum, figsize=figsize)
    hspace = 0.5 if n==4 else 0.5
    pl.subplots_adjust(bottom=0.03, top=0.97, left=0.10, right=0.98, wspace=0.3, hspace=hspace)
    axs = fig.subplots(n, 1)
    for i,letter in enumerate('ABCDEFG'[:n]):
        pl.figtext(0.01, 0.97-(0.252)*(4/n)*i, letter, fontsize=18)
    return fig, axs


def sub(data, wind=2):
    return data[::wind]


def roll(data):
    if hasattr(data, 'values'):
        data = data.values
    output = pd.Series(data).rolling(window=7).mean()
    return output


def plotter(res, refsim, fillargs, plargs, ax, key, label='Model', pct=False, ylabel=None, end=-1):
    factor = 100 if pct else 1
    date = refsim.day(res['date'])[:end]
    best = res[key].values[:end]*factor
    low  = res[key].low[:end]*factor
    high = res[key].high[:end]*factor
    ax.fill_between(date, roll(low), roll(high), **fillargs)
    ax.plot(date, roll(best), **plargs, label=label)
    if ylabel:
        ax.set_ylabel(ylabel)
    return


def format_axs(refsim, axs):
    @ticker.FuncFormatter
    def date_formatter(x, pos):
        # print(x)
        return (refsim['start_day'] + dt.timedelta(days=x)).strftime('%b-%d')
    for i,ax in enumerate(axs):
        # bbox = (0.22,1.05) # Move legend up a bit
        day_stride = 28
        xmax = ax.get_xlim()[1]
        ax.set_xticks(np.arange(0, xmax, day_stride))
        ax.xaxis.set_major_formatter(date_formatter)
        ax.legend(frameon=False, loc='upper left')
        sc.boxoff(ax=ax)
        sc.setylim(ax=ax)
        if ax.get_ylim()[1] > 1000:
            sc.commaticks(ax=ax)
    return