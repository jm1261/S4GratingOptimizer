import numpy as np
import src.analysis as anal
import matplotlib.pyplot as plt

from pathlib import Path


root = Path().absolute()
out_path = Path(f'{root}/Fano_x0')

wavelength = np.arange(-1, 1.01, 0.01)
fixed_fano_params = [0, 0.5, 1, 1, 0]
fixed_fano = [
    anal.fano_resonance(
        x=x,
        x0=fixed_fano_params[0],
        gamma=fixed_fano_params[1],
        q=fixed_fano_params[2],
        amplitude=fixed_fano_params[3],
        damping=fixed_fano_params[4])
    for x in wavelength]

moving_fano_params = [0, 0.5, 1, 1, 0]
damp_rang = np.arange(-1, 1.01, 0.1)
damp_range = [round(damp, 1) for damp in damp_rang]
moving_fanos = []
for damp in damp_range:
    moving_fano = [
        anal.fano_resonance(
            x=x,
            x0=damp,
            gamma=moving_fano_params[1],
            q=moving_fano_params[2],
            amplitude=moving_fano_params[3],
            damping=moving_fano_params[4])
        for x in wavelength]
    moving_fanos.append(moving_fano)

overlaps = []
figure_merit = []
xs = []
for index, fano in enumerate(moving_fanos):
    fig, ax1= plt.subplots(
        nrows=1,
        ncols=1,
        figsize=[10, 7])
    #overlap = np.trapz(
    #    y=np.array(fano) * np.array(fixed_fano),
    #    x=wavelength)
    overlap_difference = [((fix - mov)**2) for fix, mov in zip(fixed_fano, fano)]
    overlap = sum(overlap_difference)
    overlaps.append(overlap)
    figure_merit.append(overlap)
    xs.append(damp_range[index])
    ax1.plot(wavelength, fano, 'r', lw=2, label='Moving Fano')
    ax1.plot(wavelength, fixed_fano, 'b', lw=2, label='Fixed Fano')
    ax1.plot(wavelength, np.array(fano) * np.array(fixed_fano), 'g', lw=2, label=f'{overlap}')
    ax1.legend(loc=0, prop={'size': 14})
    plt.savefig(
        Path(f'{out_path}/x0_{damp_range[index]}.png'))
    fig.clf()
    plt.cla()
    plt.close(fig)

max_overlap = max(overlaps)
min_figure = min(figure_merit)

overlap_mark = [np.argmax(overlaps)]
figure_mark = [np.argmin(figure_merit)]
fig, ax1 = plt.subplots(
    nrows=1,
    ncols=1,
    figsize=[10, 7])
ax1.plot(xs, figure_merit, markevery=figure_mark, color='b', marker='x', markersize=10, lw=2, label=f'{min_figure}')
ax1.plot(xs, overlaps, markevery=overlap_mark, color='r', marker='x', markersize=10, lw=2, label=f'{max_overlap}')
ax1.legend(loc=0, prop={'size': 14})
plt.savefig(
    Path(f'{out_path}/x0_figure_merit.png'))
fig.clf()
plt.cla()
plt.close(fig)
