import numpy as np
import src.analysis as anal
import scipy.optimize as opt
import matplotlib.pyplot as plt

from pathlib import Path


def get_sim_fano_parameters(initial_guesses,
                            wavelength,
                            intensity,
                            sample_name,
                            plot_figure,
                            out_path=False):
    try:
        popt, pcov = opt.curve_fit(
            anal.fano_resonance,
            wavelength,
            intensity,
            initial_guesses)
        errors = np.sqrt(np.diag(pcov))
    except RuntimeError:
        popt = [0, 0, 0, 0, 0]
        errors = [0, 0, 0, 0, 0]
    return {
        f'{sample_name} S4 Fano Fit Parameters': [
            'Peak', 'Gamma', 'q', 'Ampltidude', 'Damping'],
        f'{sample_name} S4 Fano Fit': [value for value in popt],
        f'{sample_name} S4 Fano Errors': [value for value in errors]}


root = Path().absolute()
out_path = Path(f'{root}/Test')

n_range = np.arange(0.5, 5.1, 0.1)
k_range = np.arange(0, 5.1, 0.1)
period_range = range(300, 501, 20)
ff_range = np.arange(0.1, 0.91, 0.1)
grating_range = range(30, 150, 5)
film_thickness = 150
constants = [21, 1, 0, 1.4542, 0, 0, 1, film_thickness]
wavelength_range = [500, 750, 1]
guess_names = [
    'material_n', 'material_k', 'period', 'fill_factor', 'grating_thickness']
constant_names = [
    'harmonics', 'cover_n', 'cover_k', 'substrate_n', 'substrate_k', 'TE', 'TM',
    'film_thickness']
constants_dict = {
    'S4 Strings': [name for name in constant_names],
    'S4 Values': [value for value in constants]}

initial_guesses = [1.8, 0.01, 350, 0.7, 50]
n_dict = {}
for n in n_range:
    initial_guesses[0] = n
    variables_dict = {
        'S4 Strings': [name for name in guess_names],
        'S4 Guesses': [value for value in initial_guesses]}
    arg_string = anal.parameters_to_arguments(
        constants=constants_dict,
        variables=variables_dict,
        wavelength_range=wavelength_range)
    print(arg_string)
    S4_output = anal.S4_RCWA(
        lua_script='1D_4layer_grating.lua',
        argument_string=arg_string)
    wavelength, _, intensity = anal.read_S4_output(
        process_string=S4_output)
    sim_fano = get_sim_fano_parameters(
        initial_guesses=[wavelength[np.argmax(intensity)], 10, 5, 0.6, 1],
        wavelength=wavelength,
        intensity=intensity,
        sample_name='test',
        plot_figure='False',
        out_path=False)
    n_dict.update({f'{n}': sim_fano})

fig, ax1 = plt.subplots(1, figsize=[10, 7])
variables = []
x0s = []
gammas = []
qs = []
amps = []
damps = []
for variable, dicts in n_dict.items():
    fano_parameters = dicts['test S4 Fano Fit']
    variables.append(float(variable))
    x0s.append(fano_parameters[0])
    gammas.append(fano_parameters[1])
    qs.append(fano_parameters[2])
    amps.append(fano_parameters[3])
    damps.append(fano_parameters[4])

line1, = ax1.plot(variables, x0s, 'b', lw=2, label='x0')
ax1.set_xlabel('n [au]', fontsize=14, fontweight='bold')
ax1.set_ylabel('x0 [nm]', fontsize=14, fontweight='bold')
ax1.yaxis.label.set_color(line1.get_color())
ax1.tick_params(axis='y', colors=line1.get_color(), size=12)
ax1.tick_params(axis='x', colors='black', size=12)

ax1gamma = ax1.twinx()
ax1gamma.spines.right.set_position(('axes', 1))
line2, = ax1gamma.plot(variables, gammas, 'r', lw=2, label='Gamma')
ax1gamma.set_ylabel('gamma', fontsize=14, fontweight='bold')
ax1gamma.yaxis.label.set_color(line2.get_color())
ax1gamma.tick_params(axis='y', colors=line2.get_color(), size=12)

ax1q = ax1.twinx()
ax1q.spines.right.set_position(('axes', 1.2))
line3, = ax1q.plot(variables, qs, 'g', lw=2, label='q')
ax1q.set_ylabel('q', fontsize=14, fontweight='bold')
ax1q.yaxis.label.set_color(line3.get_color())
ax1q.tick_params(axis='y', colors=line3.get_color(), size=12)

ax1a = ax1.twinx()
ax1a.spines.right.set_position(('axes', 1.4))
line4, = ax1a.plot(variables, amps, 'orange', lw=2, label='amps')
ax1a.set_ylabel('amplitude', fontsize=14, fontweight='bold')
ax1a.yaxis.label.set_color(line4.get_color())
ax1a.tick_params(axis='y', colors=line4.get_color(), size=12)

ax1d = ax1.twinx()
ax1d.spines.right.set_position(('axes', 1.6))
line5, = ax1d.plot(variables, damps, 'black', lw=2, label='damps')
ax1d.set_ylabel('damping', fontsize=14, fontweight='bold')
ax1d.yaxis.label.set_color(line5.get_color())
ax1d.tick_params(axis='y', colors=line5.get_color(), size=12)

fig.tight_layout()
plt.savefig(Path(f'{out_path}/n.png'))
plt.cla()
fig.clf()
plt.close(fig)


initial_guesses = [1.8, 0.01, 350, 0.7, 50]
k_dict = {}
for k in k_range:
    initial_guesses[1] = k
    variables_dict = {
        'S4 Strings': [name for name in guess_names],
        'S4 Guesses': [value for value in initial_guesses]}
    arg_string = anal.parameters_to_arguments(
        constants=constants_dict,
        variables=variables_dict,
        wavelength_range=wavelength_range)
    print(arg_string)
    S4_output = anal.S4_RCWA(
        lua_script='1D_4layer_grating.lua',
        argument_string=arg_string)
    wavelength, _, intensity = anal.read_S4_output(
        process_string=S4_output)
    sim_fano = get_sim_fano_parameters(
        initial_guesses=[wavelength[np.argmax(intensity)], 10, 5, 0.6, 1],
        wavelength=wavelength,
        intensity=intensity,
        sample_name='test',
        plot_figure='False',
        out_path=False)
    k_dict.update({f'{k}': sim_fano})

fig, ax1 = plt.subplots(1, figsize=[10, 7])
variables = []
x0s = []
gammas = []
qs = []
amps = []
damps = []
for variable, dicts in k_dict.items():
    fano_parameters = dicts['test S4 Fano Fit']
    variables.append(float(variable))
    x0s.append(fano_parameters[0])
    gammas.append(fano_parameters[1])
    qs.append(fano_parameters[2])
    amps.append(fano_parameters[3])
    damps.append(fano_parameters[4])

line1, = ax1.plot(variables, x0s, 'b', lw=2, label='x0')
ax1.set_xlabel('k [au]', fontsize=14, fontweight='bold')
ax1.set_ylabel('x0 [nm]', fontsize=14, fontweight='bold')
ax1.yaxis.label.set_color(line1.get_color())
ax1.tick_params(axis='y', colors=line1.get_color(), size=12)
ax1.tick_params(axis='x', colors='black', size=12)

ax1gamma = ax1.twinx()
ax1gamma.spines.right.set_position(('axes', 1))
line2, = ax1gamma.plot(variables, gammas, 'r', lw=2, label='Gamma')
ax1gamma.set_ylabel('gamma', fontsize=14, fontweight='bold')
ax1gamma.yaxis.label.set_color(line2.get_color())
ax1gamma.tick_params(axis='y', colors=line2.get_color(), size=12)

ax1q = ax1.twinx()
ax1q.spines.right.set_position(('axes', 1.2))
line3, = ax1q.plot(variables, qs, 'g', lw=2, label='q')
ax1q.set_ylabel('q', fontsize=14, fontweight='bold')
ax1q.yaxis.label.set_color(line3.get_color())
ax1q.tick_params(axis='y', colors=line3.get_color(), size=12)

ax1a = ax1.twinx()
ax1a.spines.right.set_position(('axes', 1.4))
line4, = ax1a.plot(variables, amps, 'orange', lw=2, label='amps')
ax1a.set_ylabel('amplitude', fontsize=14, fontweight='bold')
ax1a.yaxis.label.set_color(line4.get_color())
ax1a.tick_params(axis='y', colors=line4.get_color(), size=12)

ax1d = ax1.twinx()
ax1d.spines.right.set_position(('axes', 1.6))
line5, = ax1d.plot(variables, damps, 'black', lw=2, label='damps')
ax1d.set_ylabel('damping', fontsize=14, fontweight='bold')
ax1d.yaxis.label.set_color(line5.get_color())
ax1d.tick_params(axis='y', colors=line5.get_color(), size=12)

fig.tight_layout()
plt.savefig(Path(f'{out_path}/k.png'))
plt.cla()
fig.clf()
plt.close(fig)

initial_guesses = [1.8, 0.01, 350, 0.7, 50]
p_dict = {}
for p in period_range:
    initial_guesses[2] = p
    variables_dict = {
        'S4 Strings': [name for name in guess_names],
        'S4 Guesses': [value for value in initial_guesses]}
    arg_string = anal.parameters_to_arguments(
        constants=constants_dict,
        variables=variables_dict,
        wavelength_range=wavelength_range)
    print(arg_string)
    S4_output = anal.S4_RCWA(
        lua_script='1D_4layer_grating.lua',
        argument_string=arg_string)
    wavelength, _, intensity = anal.read_S4_output(
        process_string=S4_output)
    sim_fano = get_sim_fano_parameters(
        initial_guesses=[wavelength[np.argmax(intensity)], 10, 5, 0.6, 1],
        wavelength=wavelength,
        intensity=intensity,
        sample_name='test',
        plot_figure='False',
        out_path=False)
    p_dict.update({f'{p}': sim_fano})

fig, ax1 = plt.subplots(1, figsize=[10, 7])
variables = []
x0s = []
gammas = []
qs = []
amps = []
damps = []
for variable, dicts in p_dict.items():
    fano_parameters = dicts['test S4 Fano Fit']
    variables.append(float(variable))
    x0s.append(fano_parameters[0])
    gammas.append(fano_parameters[1])
    qs.append(fano_parameters[2])
    amps.append(fano_parameters[3])
    damps.append(fano_parameters[4])

line1, = ax1.plot(variables, x0s, 'b', lw=2, label='x0')
ax1.set_xlabel('Period [nm]', fontsize=14, fontweight='bold')
ax1.set_ylabel('x0 [nm]', fontsize=14, fontweight='bold')
ax1.yaxis.label.set_color(line1.get_color())
ax1.tick_params(axis='y', colors=line1.get_color(), size=12)
ax1.tick_params(axis='x', colors='black', size=12)

ax1gamma = ax1.twinx()
ax1gamma.spines.right.set_position(('axes', 1))
line2, = ax1gamma.plot(variables, gammas, 'r', lw=2, label='Gamma')
ax1gamma.set_ylabel('gamma', fontsize=14, fontweight='bold')
ax1gamma.yaxis.label.set_color(line2.get_color())
ax1gamma.tick_params(axis='y', colors=line2.get_color(), size=12)

ax1q = ax1.twinx()
ax1q.spines.right.set_position(('axes', 1.2))
line3, = ax1q.plot(variables, qs, 'g', lw=2, label='q')
ax1q.set_ylabel('q', fontsize=14, fontweight='bold')
ax1q.yaxis.label.set_color(line3.get_color())
ax1q.tick_params(axis='y', colors=line3.get_color(), size=12)

ax1a = ax1.twinx()
ax1a.spines.right.set_position(('axes', 1.4))
line4, = ax1a.plot(variables, amps, 'orange', lw=2, label='amps')
ax1a.set_ylabel('amplitude', fontsize=14, fontweight='bold')
ax1a.yaxis.label.set_color(line4.get_color())
ax1a.tick_params(axis='y', colors=line4.get_color(), size=12)

ax1d = ax1.twinx()
ax1d.spines.right.set_position(('axes', 1.6))
line5, = ax1d.plot(variables, damps, 'black', lw=2, label='damps')
ax1d.set_ylabel('damping', fontsize=14, fontweight='bold')
ax1d.yaxis.label.set_color(line5.get_color())
ax1d.tick_params(axis='y', colors=line5.get_color(), size=12)

fig.tight_layout()
plt.savefig(Path(f'{out_path}/period.png'))
plt.cla()
fig.clf()
plt.close(fig)

initial_guesses = [1.8, 0.01, 350, 0.7, 50]
ff_dict = {}
for ff in ff_range:
    initial_guesses[3] = ff
    variables_dict = {
        'S4 Strings': [name for name in guess_names],
        'S4 Guesses': [value for value in initial_guesses]}
    arg_string = anal.parameters_to_arguments(
        constants=constants_dict,
        variables=variables_dict,
        wavelength_range=wavelength_range)
    print(arg_string)
    S4_output = anal.S4_RCWA(
        lua_script='1D_4layer_grating.lua',
        argument_string=arg_string)
    wavelength, _, intensity = anal.read_S4_output(
        process_string=S4_output)
    sim_fano = get_sim_fano_parameters(
        initial_guesses=[wavelength[np.argmax(intensity)], 10, 5, 0.6, 1],
        wavelength=wavelength,
        intensity=intensity,
        sample_name='test',
        plot_figure='False',
        out_path=False)
    ff_dict.update({f'{ff}': sim_fano})

fig, ax1 = plt.subplots(1, figsize=[10, 7])
variables = []
x0s = []
gammas = []
qs = []
amps = []
damps = []
for variable, dicts in ff_dict.items():
    fano_parameters = dicts['test S4 Fano Fit']
    variables.append(float(variable))
    x0s.append(fano_parameters[0])
    gammas.append(fano_parameters[1])
    qs.append(fano_parameters[2])
    amps.append(fano_parameters[3])
    damps.append(fano_parameters[4])

line1, = ax1.plot(variables, x0s, 'b', lw=2, label='x0')
ax1.set_xlabel('fill factor [au]', fontsize=14, fontweight='bold')
ax1.set_ylabel('x0 [nm]', fontsize=14, fontweight='bold')
ax1.yaxis.label.set_color(line1.get_color())
ax1.tick_params(axis='y', colors=line1.get_color(), size=12)
ax1.tick_params(axis='x', colors='black', size=12)

ax1gamma = ax1.twinx()
ax1gamma.spines.right.set_position(('axes', 1))
line2, = ax1gamma.plot(variables, gammas, 'r', lw=2, label='Gamma')
ax1gamma.set_ylabel('gamma', fontsize=14, fontweight='bold')
ax1gamma.yaxis.label.set_color(line2.get_color())
ax1gamma.tick_params(axis='y', colors=line2.get_color(), size=12)

ax1q = ax1.twinx()
ax1q.spines.right.set_position(('axes', 1.2))
line3, = ax1q.plot(variables, qs, 'g', lw=2, label='q')
ax1q.set_ylabel('q', fontsize=14, fontweight='bold')
ax1q.yaxis.label.set_color(line3.get_color())
ax1q.tick_params(axis='y', colors=line3.get_color(), size=12)

ax1a = ax1.twinx()
ax1a.spines.right.set_position(('axes', 1.4))
line4, = ax1a.plot(variables, amps, 'orange', lw=2, label='amps')
ax1a.set_ylabel('amplitude', fontsize=14, fontweight='bold')
ax1a.yaxis.label.set_color(line4.get_color())
ax1a.tick_params(axis='y', colors=line4.get_color(), size=12)

ax1d = ax1.twinx()
ax1d.spines.right.set_position(('axes', 1.6))
line5, = ax1d.plot(variables, damps, 'black', lw=2, label='damps')
ax1d.set_ylabel('damping', fontsize=14, fontweight='bold')
ax1d.yaxis.label.set_color(line5.get_color())
ax1d.tick_params(axis='y', colors=line5.get_color(), size=12)

fig.tight_layout()
plt.savefig(Path(f'{out_path}/ff.png'))
plt.cla()
fig.clf()
plt.close(fig)

initial_guesses = [1.8, 0.01, 350, 0.7, 50]
g_dict = {}
for g in grating_range:
    initial_guesses[4] = g
    variables_dict = {
        'S4 Strings': [name for name in guess_names],
        'S4 Guesses': [value for value in initial_guesses]}
    arg_string = anal.parameters_to_arguments(
        constants=constants_dict,
        variables=variables_dict,
        wavelength_range=wavelength_range)
    print(arg_string)
    S4_output = anal.S4_RCWA(
        lua_script='1D_4layer_grating.lua',
        argument_string=arg_string)
    wavelength, _, intensity = anal.read_S4_output(
        process_string=S4_output)
    sim_fano = get_sim_fano_parameters(
        initial_guesses=[wavelength[np.argmax(intensity)], 10, 5, 0.6, 1],
        wavelength=wavelength,
        intensity=intensity,
        sample_name='test',
        plot_figure='False',
        out_path=False)
    g_dict.update({f'{g}': sim_fano})

fig, ax1 = plt.subplots(1, figsize=[10, 7])
variables = []
x0s = []
gammas = []
qs = []
amps = []
damps = []
for variable, dicts in g_dict.items():
    fano_parameters = dicts['test S4 Fano Fit']
    variables.append(float(variable))
    x0s.append(fano_parameters[0])
    gammas.append(fano_parameters[1])
    qs.append(fano_parameters[2])
    amps.append(fano_parameters[3])
    damps.append(fano_parameters[4])

line1, = ax1.plot(variables, x0s, 'b', lw=2, label='x0')
ax1.set_xlabel('Grating Thickness [nm]', fontsize=14, fontweight='bold')
ax1.set_ylabel('x0 [nm]', fontsize=14, fontweight='bold')
ax1.yaxis.label.set_color(line1.get_color())
ax1.tick_params(axis='y', colors=line1.get_color(), size=12)
ax1.tick_params(axis='x', colors='black', size=12)

ax1gamma = ax1.twinx()
ax1gamma.spines.right.set_position(('axes', 1))
line2, = ax1gamma.plot(variables, gammas, 'r', lw=2, label='Gamma')
ax1gamma.set_ylabel('gamma', fontsize=14, fontweight='bold')
ax1gamma.yaxis.label.set_color(line2.get_color())
ax1gamma.tick_params(axis='y', colors=line2.get_color(), size=12)

ax1q = ax1.twinx()
ax1q.spines.right.set_position(('axes', 1.2))
line3, = ax1q.plot(variables, qs, 'g', lw=2, label='q')
ax1q.set_ylabel('q', fontsize=14, fontweight='bold')
ax1q.yaxis.label.set_color(line3.get_color())
ax1q.tick_params(axis='y', colors=line3.get_color(), size=12)

ax1a = ax1.twinx()
ax1a.spines.right.set_position(('axes', 1.4))
line4, = ax1a.plot(variables, amps, 'orange', lw=2, label='amps')
ax1a.set_ylabel('amplitude', fontsize=14, fontweight='bold')
ax1a.yaxis.label.set_color(line4.get_color())
ax1a.tick_params(axis='y', colors=line4.get_color(), size=12)

ax1d = ax1.twinx()
ax1d.spines.right.set_position(('axes', 1.6))
line5, = ax1d.plot(variables, damps, 'black', lw=2, label='damps')
ax1d.set_ylabel('damping', fontsize=14, fontweight='bold')
ax1d.yaxis.label.set_color(line5.get_color())
ax1d.tick_params(axis='y', colors=line5.get_color(), size=12)

fig.tight_layout()
plt.savefig(Path(f'{out_path}/grating.png'))
plt.cla()
fig.clf()
plt.close(fig)
