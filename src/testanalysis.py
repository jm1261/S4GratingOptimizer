import os
import subprocess
import numpy as np
import scipy.optimize as opt

from io import StringIO
from pathlib import Path
from src.plotting import fanofitplot, fomplot
from src.fileIO import save_json_dicts, load_json
from src.fileIO import S4_args, load_fomlog, log_figure_of_merit


def no_damping_fano(x,
                    x0,
                    amplitude,
                    assymetry,
                    gamma,
                    offset):
    '''
    No damping fano resonance equation from literature.
    Args:
        x: <float> wavelength/energy value at which to analyse the fano
        x0: <float> peak wavelength/energy value of the fano resonance
        amplitude: <float> fano resonance amplitude
        assymetry: <float> assymetry factor
        gamma: <float> peak width factor
        offset: <float> DC offset
    Returns:
        y: <float> fano resonance value at x
    '''
    numerator = ((assymetry * gamma) + (x - x0)) ** 2
    denominator = (gamma ** 2) + ((x - x0) ** 2)
    y = (amplitude * (numerator / denominator)) + offset
    return y


def fano_resonances(x,
                    x0,
                    gamma,
                    q,
                    amplitude,
                    damping):
    '''
    Fano resonance peak equation.
    Args:
        x: <array> x-axis point, wavelength in nm
        x0: <float> peak x point, wavelength in nm
        gamma: <float> reduced frequency
        q: <float> shape factor
        amplitude: <float> peak amplitude
        damping: <float> damping factor
    Returns:
        y: <array> data array for fano peak
    '''
    omega = 2 * ((x - x0) / gamma)
    numerator = ((q + omega) ** 2) + damping
    denominator = 1 + (omega ** 2)
    y = amplitude * (numerator / denominator)
    return y


def normalise_intensity(intensity):
    '''
    Normalise intensity to minimum intensity value.
    Args:
        intensity: <array> intensity array
    Returns:
        normalised_intensity: <array> normalised intensity array
    '''
    minimum_intensity = min(intensity)
    normalised_intensity = [i - minimum_intensity for i in intensity]
    return normalised_intensity


def simulation_wavelength_range(experimental_wavelength,
                                experimental_intensity):
    '''
    Generate simulation wavelength range to match experimentally measured
    wavelength range.
    Args:
        experimental_wavelength: <array> experimental wavelength array
        experimental_intensity: <array> experimental intensity array
    Returns:
        wavelengths: <array> new experimental wavelength range
        intensity: <array> experimental intensity interpolated at wavelength
                    range
        simulation_wavelength: <array> [start, stop, step]
    '''
    wavelength_range = [
        int(round(min(experimental_wavelength), 0)),
        int(round(max(experimental_wavelength), 0)),
        1]
    wavelengths = np.arange(
        wavelength_range[0],
        wavelength_range[1] + (wavelength_range[2] / 10),
        wavelength_range[2])
    intensity = np.interp(
        x=wavelengths,
        xp=experimental_wavelength,
        fp=experimental_intensity)
    simulation_wavelength = [
        wavelength_range[0] - 25,
        wavelength_range[1] + 25,
        wavelength_range[2]]
    return wavelengths, intensity, simulation_wavelength


def argument_string(names,
                    values):
    '''
    Create argument string for S4 variables. Variable names and parameters must
    be in corresponding order within arrays.
    Args:
        names: <array> array of variable names that match variables in S4
                lua script
        values: <array> array of matching values for argument names
    Returns:
        arg_string: <string> '; ' joined argument string for S4 variables
    '''
    arguments = [
        f'{name} = {value}'
        for name, value in zip(names, values)]
    arg_string = ('; ').join(arguments)
    return arg_string


def parameters_to_arguments(constants,
                            variables,
                            wavelength_range):
    '''
    Turn constants and variables dictionaries to argument names and values.
    Args:
        constants: <dict> simulation constants dictionary
        variables: <dict> simulation variables dictionary
        wavelength_range: <array> [start, stop, step]
    Returns:
        names: <array> argument names
        values: <array> argument values
    '''
    names = []
    values = []
    for index, name in enumerate(constants['S4 Strings']):
        if name == 'peak_wavelength':
            pass
        else:
            names.append(name)
            values.append((constants['S4 Values'])[index])
    for index, name in enumerate(variables['S4 Strings']):
        if name == 'peak_wavelength':
            pass
        else:
            names.append(name)
            values.append((variables['S4 Guesses'])[index])
    wavelength_names = [
        'wavelength_initial',
        'wavelength_final',
        'wavelength_step']
    for name, value in zip(wavelength_names, wavelength_range):
        names.append(name)
        values.append(value)
    arg_string = argument_string(
        names=names,
        values=values)
    return arg_string


def S4_RCWA(lua_script,
            argument_string):
    '''
    Run S4 RCWA simulation with desired argument string and lua script.
    Args:
        lua_script: <string> name of target lua script for S4, must be in same
                    directory as python script
        argument_string: <string> argument string in corrector format with
                            variable and parameter names and values
    Returns:
        process: <string> subprocess PIPE string output from S4
    '''
    command = f'S4 -a "{argument_string}" {lua_script}'
    process = subprocess.run(
        command,
        shell=True,
        stdout=subprocess.PIPE)
    return process


def read_S4_output(process_string):
    '''
    Read S4 string output from subprocess stdout. Is dependent on lua script
    output.
    Args:
        process_string: <string> subprocess PIPE stdout string
    Returns:
        wavelength: <array> wavelength array (nm)
        transmission: <array> transmission array (au)
        reflection: <array> reflection array (au)
    '''
    try:
        wavelength, transmission, reflection = np.genfromtxt(
            fname=StringIO(
                process_string.stdout.decode('utf-8')),
            delimiter='\t',
            unpack=True)
    except:
        wavelength, transmission, reflection = np.genfromtxt(
            fname=StringIO(
                process_string.stdout.decode('utf-8')),
            delimiter=',',
            unpack=True)
    return wavelength, transmission, reflection


def get_exp_fano_parameters(initial_guesses,
                            wavelength,
                            intensity,
                            sample_name,
                            plot_figure,
                            out_path):
    '''
    Use scipy optimize and fano_resonance equations to iterate and find the
    peak in intensity and the corresponding wavelength and error value. Tries to
    fit a fano peak with damping, taken from literature, if that fails will fit
    a non-damping fano taken from literature.
    Args:
        initial_guesses: <array> experimental fano parameters initial guesses
        wavelength: <array> wavelength array in nm
        intensity: <array> intensity array
        sample_name: <string> secondaray sample identifier string
        plot_figure: <string> if "True" will plot normalised spectrum with fano
        out_path: <string> path to save
    Returns:
        results_dictionary: <dict> dictionary containing:
            popt: <array> fano fit parameters
            pcov: <array> fano fit errors
        normalised_intensity: <array> experimental intensity normalised
        damping: <bool> if True, damping fano used
    '''
    normalised_intensity = normalise_intensity(intensity=intensity)
    try:
        popt, pcov = opt.curve_fit(
            fano_resonances,
            wavelength,
            normalised_intensity,
            initial_guesses)
        errors = np.sqrt(np.diag(pcov))
    except RuntimeError:
        print('\nRun Time Error No Peak Found')
        popt = [0, 0, 0, 0, 0]
        errors = [0, 0, 0, 0, 0]
    except ValueError:
        print('\nValue Error No Peak Found')
        popt = [0, 0, 0, 0, 0]
        errors = [0, 0, 0, 0, 0]
    damping = True
    fano_parameters = ['Peak', 'Gamma', 'q', 'Amplitude', 'Damping']
    if sum(popt) == 0:
        try:
            popt, pcov = opt.curve_fit(
                no_damping_fano,
                wavelength,
                normalised_intensity,
                initial_guesses)
            errors = np.sqrt(np.diag(pcov))
        except RuntimeError:
            print('\nRun Time Error No Peak Found')
            popt = [0, 0, 0, 0, 0]
            errors = [0, 0, 0, 0, 0]
        except ValueError:
            print('\nValue Error No Peak Found')
            popt = [0, 0, 0, 0, 0]
            errors = [0, 0, 0, 0, 0]
        if plot_figure == 'True':
            fanofitplot(
                wavelength=wavelength,
                intensity=normalised_intensity,
                fano=[
                    no_damping_fano(
                        x=wav,
                        x0=popt[0],
                        amplitude=popt[1],
                        assymetry=popt[2],
                        gamma=popt[3],
                        offset=popt[4])
                    for wav in wavelength],
                peak_wavelength=popt[0],
                peak_error=errors[0],
                out_path=Path(f'{out_path}/{sample_name}_Experimental.png'))
        damping = False
        fano_parameters = ['Peak', 'Amplitude', 'Assymetry', 'Gamma', 'Offset']
    else:
        if plot_figure == 'True':
            fanofitplot(
                wavelength=wavelength,
                intensity=normalised_intensity,
                fano=[
                    fano_resonances(
                        x=wav,
                        x0=popt[0],
                        gamma=popt[1],
                        q=popt[2],
                        amplitude=popt[3],
                        damping=popt[4])
                    for wav in wavelength],
                peak_wavelength=popt[0],
                peak_error=errors[0],
                out_path=Path(f'{out_path}/{sample_name}_Experimental.png'))
    results_dictionary = {
        f'{sample_name} Fano Fit Parameters': [fp for fp in fano_parameters],
        f'{sample_name} Fano Fit': [value for value in popt],
        f'{sample_name} Fano Errors': [err for err in errors]}
    return results_dictionary, normalised_intensity, damping


def get_S4_fano(wavelength,
                intensity,
                initial_guesses,
                sample_name,
                plot_figure,
                out_path,
                damping=False):
    '''
    Use scipy optimize and fano_resonance equations to iterate and find the
    peak in intensity and the corresponding wavelength and error value. Tries to
    fit a fano peak with damping, taken from literature, if that fails will fit
    a non-damping fano taken from literature.
    Args:
        wavelength: <array> wavelength array in nm
        intensity: <array> intensity array
        initial_guesses: <array> initial guesses for fano resonance equation
        sample_name: <string> sample name identifier string
        plot_figure: <string> if "True" will plot normalised spectrum with fano
        out_path: <string> path to save
        damping: <bool> if True, uses fano resonance equation with damping
    Returns:
        results_dictionary: <dict> dictionary containing:
            popt: <array> S4 fano fit parameters
            pcov: <array> S4 fano fit errors
    '''
    if damping:
        try:
            popt, pcov = opt.curve_fit(
                fano_resonances,
                wavelength,
                intensity,
                initial_guesses)
            errors = np.sqrt(np.diag(pcov))
        except RuntimeError:
            print('\nRun Time Error \nNo Peak Found')
            popt = [0, 0, 0, 0, 0]
            errors = [0, 0, 0, 0, 0]
        except ValueError:
            print('\nValue Error \nNo Peak Found')
            popt = [0, 0, 0, 0, 0]
            errors = [0, 0, 0, 0, 0]
        fano_parameters = ['Peak', 'Gamma', 'q', 'Amplitude', 'Damping']
        if plot_figure == 'True':
            fanofitplot(
                wavelength=wavelength,
                intensity=intensity,
                fano=[
                    fano_resonances(
                        x=wav,
                        x0=popt[0],
                        gamma=popt[1],
                        q=popt[2],
                        amplitude=popt[3],
                        damping=popt[4])
                    for wav in wavelength],
                peak_wavelength=popt[0],
                peak_error=errors[0],
                out_path=Path(f'{out_path}/{sample_name}_Simulation.png'))
    else:
        try:
            popt, pcov = opt.curve_fit(
                no_damping_fano,
                wavelength,
                intensity,
                initial_guesses)
            errors = np.sqrt(np.diag(pcov))
        except RuntimeError:
            print('\nRun Time Error \nNo Peak Found')
            popt = [0, 0, 0, 0, 0]
            errors = [0, 0, 0, 0, 0]
        except ValueError:
            print('\nValue Error \nNo Peak Found')
            popt = [0, 0, 0, 0, 0]
            errors = [0, 0, 0, 0, 0]
        fano_parameters = ['Peak', 'Amplitude', 'Assymetry', 'Gamma', 'Offset']
        if plot_figure == 'True':
            fanofitplot(
                wavelength=wavelength,
                intensity=intensity,
                fano=[
                    no_damping_fano(
                        x=wav,
                        x0=popt[0],
                        amplitude=popt[1],
                        assymetry=popt[2],
                        gamma=popt[3],
                        offset=popt[4])
                    for wav in wavelength],
                peak_wavelength=popt[0],
                peak_error=errors[0],
                out_path=Path(f'{out_path}/{sample_name}_Simulation.png'))
    return {
        f'{sample_name} S4 Fano Fit Parameters': [fp for fp in fano_parameters],
        f'{sample_name} S4 Fano Fit': [value for value in popt],
        f'{sample_name} S4 Fano Errors': [err for err in errors]}


def compare_fano_parameters(experimental_fano,
                            simulation_fano):
    '''
    Compare experimental fano resonance to simulation fano resonance results.
    Use the square of the difference between values to produce a figure of
    merit.
    Args:
        experimental_fano: <array> experimental fano parameters
        simulation_fano: <array> simulation fano parameters
    Returns:
        figure_merit: <float> figure of merit
        differences: <array> differences between each individual parameter
    '''
    differences = [
        ((exp - sim) ** 2)
        for exp, sim
        in zip(experimental_fano, simulation_fano)]
    figure_merit = sum(differences)
    return figure_merit, differences


def compare_fano_curves(experimental_intensity,
                        simulation_intensity):
    '''
    Compare experimental intensity and simulation intensity overlap. Make sure
    both intensity arrays are the same length.
    Args:
        experimental_intensity: <array> experimental intensity array
        simulation_intensity: <array> simulation intensity array
    Return:
        figure_merit: <float> figure of merit
    '''
    overlap = [
        ((exp - sim) ** 2)
        for exp, sim in zip(experimental_intensity, simulation_intensity)]
    figure_merit = sum(overlap)
    return figure_merit


def compare_experimental_simulation(sample_name,
                                    experimental_fano,
                                    simulation_fano,
                                    experimental_intensity,
                                    simulation_intensity):
    '''
    Compare experimental and simulation fano fits and curves to produce a single
    figure of merit.
    Args:
        sample_name: <string> sample name identifier string
        experimental_fano: <array> experimental fano parameters
        simulation_fano: <array> simulation fano parameters
        experimental_intensity: <array> experimental normalised intensity
        simulation_intensity: <array> simulation normalised intensity
    Returns:
        results_dict: <dict> results dictionary containing:
            figure of merit,
            fano figure of merit
            overlap figure of merit
            differences array
    '''
    if sum(simulation_fano) == 0:
        figure_merit = (500 ** 2)
        figure_merit_fano = simulation_fano
        figure_merit_overlap = (500 ** 2)
        differences = 0
    else:
        figure_merit_fano, differences = compare_fano_parameters(
            experimental_fano=experimental_fano,
            simulation_fano=simulation_fano)
        figure_merit_overlap = compare_fano_curves(
            experimental_intensity=experimental_intensity,
            simulation_intensity=simulation_intensity)
        figure_merit = figure_merit_fano + figure_merit_overlap
    return {
        f'{sample_name} Figure Of Merit': figure_merit,
        f'{sample_name} Fano Error': figure_merit_fano,
        f'{sample_name} Overlap Error': figure_merit_overlap,
        f'{sample_name} Differences': differences}


def create_variables_dictionary(iteration_variable,
                                iteration_variable_name,
                                iteration_constant_names,
                                iteration_constant_values,
                                index):
    '''
    Create variables dictionary for individual variable optimization, takes
    the iteration variable and variable name and adds back into the variables
    dictionary with the currently constant variables and names.
    Args:
        iteration_variable: <float> iteration variable value
        iteration_variable_name: <string> iteration variable name
        iteration_constant_names: <array> variable names for constant variables
        iteration_constant_values: <array> variable values for constant
                                    variables
        index: <int> index for S4 variables currently optimizing
    Returns:
        variables: <dict> reconstructed variables dict
    '''
    strings = []
    values = []
    for i, name in enumerate(iteration_constant_names):
        if i == index:
            strings.append(iteration_variable_name)
            values.append(iteration_variable)
            strings.append(name)
            values.append(iteration_constant_values[i])
        else:
            strings.append(name)
            values.append(iteration_constant_values[i])
    if index == len(iteration_constant_names):
        strings.append(iteration_variable_name)
        values.append(iteration_variable)
    variables = {
        'S4 Strings': [name for name in strings],
        'S4 Guesses': [value for value in values]}
    print(f'variables = {variables}')
    return variables


def trims_sim_intensity(simulation_wavelength,
                        simulation_intensity,
                        experimental_wavelength):
    '''
    Re-trim simulation intensity to match experimental intensity length.
    Args:
        simulation_wavelength: <array> simulation wavelength array
        simulation_intensity: <array> simulation intensity array
        experimental_wavelength: <array> experimental wavelength array
    Returns:
        trimmed_intensity: <array> re-trimmed simulation intensity array
    '''
    wavelengths_min = np.abs(
        simulation_wavelength - min(experimental_wavelength))
    wavelengths_max = np.abs(
        simulation_wavelength - max(experimental_wavelength))
    wavelength_argmin = np.argmin(wavelengths_min)
    wavelength_argmax = np.argmin(wavelengths_max)
    trimmed_intensity = simulation_intensity[
        wavelength_argmin: wavelength_argmax + 1]
    return trimmed_intensity


def figure_of_merit(variables,
                    constants,
                    sample_name,
                    batch_name,
                    peak_parameters,
                    lua_script,
                    plot_figure,
                    out_path,
                    log_fom):
    '''
    Calculate figure of merit between experimental and simulation fano resonance
    parameters and intensity curves while optimizing variable(s). Builds arg
    string for S4 simulation, processes S4 output, and calculates fano resonance
    parameters. Log iterations if desired.
    Args:
        variables: <dict> simulation variables dictionary
        constants: <dict> simulation constants dictionary
        sample_name: <string> sample name identifier string
        batch_name: <string> batch name identifier string
        peak_parameters: <dict> peak parameters dictionary from measurements
                        containing experimental wavelength, fano fit, and errors
        lua_script: <string> lua script file
        plot_figure: <string> if "True" will plot simulation data
        out_path: <string> path to save results
        log_fom: <string> if "True" will log iteration data
    Returns:
        figure_merit: <float> individual figure of merit float for optimizer
    '''
    exp_wav, exp_intensity, sim_wav = simulation_wavelength_range(
        experimental_wavelength=peak_parameters[f'{sample_name} Wavelength'],
        experimental_intensity=peak_parameters[f'{sample_name} Intensity'])
    argument_string = parameters_to_arguments(
        constants=constants,
        variables=variables,
        wavelength_range=sim_wav)
    S4_output = S4_RCWA(
        lua_script=lua_script,
        argument_string=argument_string)
    if constants['Reflection'] == 'True':
        sim_wavelength, _, sim_intensity = read_S4_output(
            process_string=S4_output)
    else:
        sim_wavelength, sim_intensity, _ = read_S4_output(
            process_string=S4_output)
    exp_fano, exp_norm_intensity, damping = get_exp_fano_parameters(
        initial_guesses=peak_parameters[f'{sample_name} Fano Fit'],
        wavelength=exp_wav,
        intensity=exp_intensity,
        sample_name=sample_name,
        plot_figure=plot_figure,
        out_path=out_path)
    sim_fano = get_S4_fano(
        wavelength=sim_wavelength,
        intensity=sim_intensity,
        initial_guesses=exp_fano[f'{sample_name} Fano Fit'],
        sample_name=sample_name,
        plot_figure=plot_figure,
        out_path=out_path,
        damping=damping)
    trimmed_intensity = trims_sim_intensity(
        simulation_wavelength=sim_wavelength,
        simulation_intensity=sim_intensity,
        experimental_wavelength=exp_wav)
    figure_of_merit = compare_experimental_simulation(
        sample_name=sample_name,
        experimental_fano=exp_fano[f'{sample_name} Fano Fit'],
        simulation_fano=sim_fano[f'{sample_name} S4 Fano Fit'],
        experimental_intensity=exp_norm_intensity,
        simulation_intensity=trimmed_intensity)
    results_dict = dict(
        exp_fano,
        **sim_fano,
        **figure_of_merit)
    results_dict.update({f'{sample_name} Variables': variables})
    results_dict.update({f'{sample_name} Constants': constants})
    results_dict.update(
        {f'{sample_name} Wavelength Range': sim_wav})
    save_json_dicts(
        out_path=Path(f'{out_path}/Working_Result.json'),
        dictionary=results_dict)
    if log_fom == 'True':
        log_figure_of_merit(
            variables=variables['S4 Guesses'],
            figure_of_merit=figure_of_merit[f'{sample_name} Figure Of Merit'],
            simulation_fano=sim_fano[f'{sample_name} S4 Fano Fit'],
            experimental_fano=exp_fano[f'{sample_name} Fano Fit'],
            simulation_intensity=trimmed_intensity,
            experimental_intensity=exp_norm_intensity,
            sample_name=sample_name,
            batch_name=f'{batch_name}',
            out_path=out_path)
    return figure_of_merit[f'{sample_name} Figure Of Merit']


def compute_individual_figure_merit(iteration_variable,
                                    iteration_variable_name,
                                    iteration_constant_names,
                                    iteration_constant_values,
                                    index,
                                    constants,
                                    sample_name,
                                    batch_name,
                                    peak_parameters,
                                    lua_script,
                                    plot_figure,
                                    out_path,
                                    log_fom):
    '''
    Calculate figure of merit between experimental and simulation fano resonance
    parameters and intensity curves while optimizing 1 individual variable.
    Args:
        iteration_variable: <float> iteration variable value
        iteration_variable_name: <string> iteration variable name
        iteration_constant_names: <array> variable names for constant variables
        iteration_constant_values: <array> variable values for constant
                                    variables
        index: <int> index for S4 variables currently optimizing
        constants: <dict> simulation constants dictionary
        sample_name: <string> sample name identifier string
        batch_name: <string> batch name identifier string
        peak_parameters: <dict> peak parameters dictionary from measurements
                        containing experimental wavelength, fano fit, and errors
        lua_script: <string> lua script file
        plot_figure: <string> if "True" will plot simulation data
        out_path: <string> path to save results
        log_fom: <string> if "True" will log iteration data
    Returns:
        figure_merit: <float> individual figure of merit float for optimizer
    '''
    variables = create_variables_dictionary(
        iteration_variable=iteration_variable[0],
        iteration_variable_name=iteration_variable_name,
        iteration_constant_names=iteration_constant_names,
        iteration_constant_values=iteration_constant_values,
        index=index)
    figure_merit = figure_of_merit(
        variables=variables,
        constants=constants,
        sample_name=sample_name,
        batch_name=batch_name,
        peak_parameters=peak_parameters,
        lua_script=lua_script,
        plot_figure=plot_figure,
        out_path=out_path,
        log_fom=log_fom)
    return figure_merit


def compute_figure_of_merit(variables,
                            variables_names,
                            constants,
                            sample_name,
                            batch_name,
                            peak_parameters,
                            lua_script,
                            plot_figure,
                            out_path,
                            log_fom):
    '''
    Calculate figure of merit between experimental and simulation fano resonance
    parameters and intensity curves while optimizing variables.
    Args:
        variables: <array> iteration variable values
        variables_names: <array> iteration variable names
        constants: <dict> simulation constants dictionary
        sample_name: <string> sample name identifier string
        batch_name: <string> batch name identifier string
        peak_parameters: <dict> peak parameters dictionary from measurements
                        containing experimental wavelength, fano fit, and errors
        lua_script: <string> lua script file
        plot_figure: <string> if "True" will plot simulation data
        out_path: <string> path to save results
        log_fom: <string> if "True" will log iteration data
    Returns:
        figure_merit: <float> individual figure of merit float for optimizer
    '''
    variables_dict = {
        'S4 Strings': [name for name in variables_names],
        'S4 Guesses': [guess for guess in variables]}
    print(f'variables = {variables_dict}')
    figure_merit = figure_of_merit(
        variables=variables_dict,
        constants=constants,
        sample_name=sample_name,
        batch_name=batch_name,
        peak_parameters=peak_parameters,
        lua_script=lua_script,
        plot_figure=plot_figure,
        out_path=out_path,
        log_fom=log_fom)
    return figure_merit


def opt_individual_variable(variables,
                            constants,
                            sample_name,
                            batch_name,
                            peak_parameters,
                            lua_script,
                            plot_figure,
                            out_path,
                            log_fom,
                            tolerance):
    '''
    Optimize the input variables piecewise for a given tolerance using the
    compute_individual_figure_merit function.
    Args:
        variables: <dict> simulation variables dictionary
        constants: <dict> simulation constants dictionary
        sample_name: <string> sample name identifier string
        batch_name: <string> batch name identifier string
        peak_parameters: <dict> peak parameters dictionary from measurements
                        containing experimental wavelength, fano fit, and errors
        lua_script: <string> lua script file
        plot_figure: <string> if "True" will plot simulation data
        out_path: <string> path to save results
        log_fom: <string> if "True" will log iteration data
        tolerance: <float> figure of merit tolerance, minimum shift in figure of
                    merit before optimizer considers it's done (scientific
                    format, i.e. 1E-2)
    Returns:
        results_dict: <dictionary>
            updated variables dictionary containing S4 strings, optimized values
            , error bounds, and precision for further optimization or result
            output
        iterations: <float> iteration count
    '''
    variable_names = variables['S4 Strings']
    variable_guesses = variables['S4 Guesses']
    variable_bounds = tuple([tuple(p) for p in variables['S4 Error Bounds']])
    variable_precision = variables['S4 Precision']
    optimizer_iterations = []
    optimizer_errors = []
    for index, variable in enumerate(variable_guesses):
        print(f'Optimizing: {variable_names[index]}, tolerance = {tolerance}')
        iteration_constant_values = [
            value for value in variable_guesses if value != variable]
        iteration_constant_names = [
            name for name in variable_names if name != variable_names[index]]
        iteration_variable_name = variable_names[index]
        optimizer_results = opt.minimize(
            fun=compute_individual_figure_merit,
            x0=[variable],
            args=(
                iteration_variable_name,
                iteration_constant_names,
                iteration_constant_values,
                index,
                constants,
                sample_name,
                batch_name,
                peak_parameters,
                lua_script,
                plot_figure,
                out_path,
                log_fom),
            bounds=(variable_bounds[index], ),
            method='L-BFGS-B',
            options={
                'ftol': tolerance,
                'gtol': tolerance,
                'eps': variable_precision[index],
                'maxiter': 1000})
        optimizer_hessian = np.diag(optimizer_results.hess_inv.todense())
        hessian_error = ([
            tolerance * optimizer_results.fun * hess
            for hess in optimizer_hessian])
        variable_guesses[index] = (optimizer_results.x)[0]
        optimizer_iterations.append(optimizer_results.nit)
        optimizer_errors.append(hessian_error)
        print(
            f'Optimized: {variable_names[index]}, '
            f'fom = {optimizer_results.fun}\n')
    iterations = sum(optimizer_iterations)
    results_dict = {
        'S4 Strings': [name for name in variable_names],
        'S4 Guesses': [guess for guess in variable_guesses],
        'S4 Error Bounds': [pair for pair in variables['S4 Error Bounds']],
        'S4 Precision': [precision for precision in variable_precision],
        'FOM': optimizer_results.fun}
    return results_dict, iterations


def optimize_gratings_variables(variables,
                                constants,
                                sample_name,
                                batch_name,
                                peak_parameters,
                                lua_script,
                                plot_figure,
                                out_path,
                                log_fom,
                                tolerance,
                                iterations=0):
    '''
    Optimize S4 fano parameters using scipy optimize minimize to vary simulation
    variables and calculate a minimum figure of merit with respect to the
    experimentally measured parameters. Return all parameters to one results
    dictionary.
    Args:
        variables: <dict> simulation variables dictionary
        constants: <dict> simulation constants dictionary
        sample_name: <string> sample name identifier string
        batch_name: <string> batch name identifier string
        peak_parameters: <dict> experimental peak dictionary
        lua_script: <string> lua script file
        plot_figure: <string> if "True" will plot simulation output
        out_path: <string> path to save results
        log_fom: <string> if "True" log iteration variables
        tolerance: <float> figure of merit tolerance, minimum shift in figure of
                    merit before optimizer considers it's done (scientific
                    format, i.e. 1E-2)
        iterations: <int> number of previous iterations, if applicable
    Returns:
        results_dictionary: <dict> grating results containing:
            variable results from optimizer
            simulation constants
            simulation fano results
    '''
    print(f'Optimizing {batch_name} {sample_name}')
    variable_names = variables['S4 Strings']
    variable_guesses = variables['S4 Guesses']
    variable_bounds = tuple([tuple(p) for p in variables['S4 Error Bounds']])
    optimizer_results = opt.minimize(
        fun=compute_figure_of_merit,
        x0=variable_guesses,
        args=(
            variable_names,
            constants,
            sample_name,
            batch_name,
            peak_parameters,
            lua_script,
            plot_figure,
            out_path,
            log_fom),
        bounds=variable_bounds,
        method='L-BFGS-B',
        options={
            'ftol': tolerance,
            'gtol': tolerance,
            'maxiter': 2000})
    optimizer_hessian = np.diag(optimizer_results.hess_inv.todense())
    optimizer_errors = ([
        tolerance * optimizer_results.fun * hess
        for hess in optimizer_hessian])
    simulation_path = Path(f'{out_path}/Working_Result.json')
    simulation_fano = load_json(file_path=simulation_path)
    variable_results = {
        f'{sample_name} Optimizer Results': [x for x in optimizer_results.x],
        f'{sample_name} Optimizer Errors': [x for x in optimizer_errors],
        f'{sample_name} Iterations': optimizer_results.nit + iterations}
    results_dict = dict(
        variable_results,
        **simulation_fano)
    os.remove(simulation_path)
    print(f'Optimized {batch_name} {sample_name}\n')
    return results_dict


def optimize_S4_grating(parameters_path,
                        measured_paths,
                        batch_name,
                        grating_name,
                        peak_parameters,
                        lua_script,
                        plot_figure,
                        out_path,
                        log_fom):
    '''
    Optimize S4 grating with given variables, arguments, constants, and names.
    Args:
        parameters_path: <string> path to S4_parameters.json file
        measured_paths: <d9ict> paths to measured parameter json files, can
                        be null
        batch_name: <string> batch name identifier string
        grating_name: <string> grating name identifier string
        peak_parameters: <dict> dictionary containing fano fit parameters,
                        fano fit errors, wavelength, and normalised intensity
                        with key == grating name, see readme file
        lua_script: <string> name of lua script for grating
        plot_figure: <string> 'True' or 'False', will plot all output files if
                    'True'
        out_path: <string> path to save out results and plots
        log_fom: <string> 'True' or 'False', will log figure of merit
                optimization if 'True'
    Returns:
        grating_results: <dict> results dictionary
    '''
    arguments, constants, variables = S4_args(
        parameters_path=parameters_path,
        file_paths=measured_paths,
        batch_name=f'{batch_name}',
        grating_name=f'{grating_name}')
    figure_merit = []
    iteration_count = 0
    if arguments['All Arguments'] == 'True':
        optimized_variables = {}
        while figure_merit[0] > 1:
            variables_dict1, iterations1 = opt_individual_variable(
                variables=variables,
                constants=constants,
                sample_name=f'{grating_name}',
                batch_name=f'{batch_name}',
                peak_parameters=peak_parameters,
                lua_script=lua_script,
                plot_figure=plot_figure,
                out_path=out_path,
                log_fom=log_fom,
                tolerance=1E-2)
            figure_merit[0] = variables_dict1['FOM']
            variables_dict2, iterations2 = opt_individual_variable(
                variables=variables_dict1,
                constants=constants,
                sample_name=f'{grating_name}',
                batch_name=f'{batch_name}',
                peak_parameters=peak_parameters,
                lua_script=lua_script,
                plot_figure=plot_figure,
                out_path=out_path,
                log_fom=log_fom,
                tolerance=1E-4)
            figure_merit[0] = variables_dict2['FOM']
            variables_dict3, iterations3 = opt_individual_variable(
                variables=variables_dict2,
                constants=constants,
                sample_name=f'{grating_name}',
                batch_name=f'{batch_name}',
                peak_parameters=peak_parameters,
                lua_script=lua_script,
                plot_figure=plot_figure,
                out_path=out_path,
                log_fom=log_fom,
                tolerance=1E-6)
            figure_merit[0] = variables_dict3['FOM']
        grating_results = optimize_gratings_variables(
            variables=variables_dict3,
            constants=constants,
            sample_name=f'{grating_name}',
            batch_name=f'{batch_name}',
            peak_parameters=peak_parameters,
            lua_script=lua_script,
            plot_figure=plot_figure,
            out_path=out_path,
            log_fom=log_fom,
            tolerance=1E-6,
            iterations=iterations1 + iterations2 + iterations3)
        print(grating_results)
        if log_fom == 'True':
            var, fom, fanos, sim_fano, exp_fano, simint, expint = load_fomlog(
                directory_path=out_path,
                sample_name=f'{grating_name}',
                batch_name=f'{batch_name}',
                length_variables=len(variables['S4 Strings']))
            fomplot(
                wavelengths=grating_results[f'{grating_name} Wavelength Range'],
                variables=var,
                figure_of_merit=fom,
                fano_parameters=fanos,
                simulation_fano=sim_fano,
                experimental_fano=exp_fano,
                simulation_intensity=simint,
                experimental_intensity=expint,
                variable_names=variables['S4 Strings'],
                fano_names=grating_results[
                    f'{grating_name} S4 Fano Fit Parameters'],
                out_path=Path(f'{out_path}/{batch_name}_{grating_name}.png'))
    else:
        grating_results = {
            f'{grating_name} Missing Parameters':
            [f'{argument}' for argument in arguments["Missing Arguments"]]}
        print(grating_results)
    return grating_results
