import os
import subprocess
import numpy as np
import scipy.optimize as opt

from io import StringIO
from pathlib import Path
from src.plotting import fanofitplot, fomplot
from src.fileIO import save_json_dicts, load_json
from src.fileIO import S4_args, load_fomlog, log_figure_of_merit


def no_damping_fano(x, x0, gamma, q, amplitude):
    '''
    Fano resonance peak equation.
    Args:
        x: <array> x-axis point, wavelength in nm
        x0: <float> peak x point, wavelength in nm
        gamma: <float> reduced frequency
        q: <float> shape factor
        amplitude: <float> peak amplitude
    Returns:
        y: <array> data array for fano peak
    '''
    omega = 2 * ((x - x0) / gamma)
    numerator = ((q + omega) ** 2)
    denominator = 1 + (omega ** 2)
    y = amplitude * (numerator / denominator)
    return y


def fano_resonance(x, x0, gamma, q, amplitude, damping):
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
    Normalise intensity to an average value.
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
        wavelength_range: <array> [start, stop, step]
        intensity: <array> experimental_intensity interpolated at wavelength
                            range
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
    return wavelengths, intensity


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
                            out_path=False):
    '''
    Use scipy optimise and fano_resonance to iterate and find the peak in
    intensity and the corresponding wavelength and error value. Tries to fit a
    fano peak with 5 parameters (excluding damping), if not possible then will
    fit with 5 parameters (including damping).
    Args:
        initial_guesses: <array> experimental fano parameters initial guesses
        wavelength: <array> wavelength array in nm
        intensity: <array> intensity array
        sample_name: <string> secondary sample identifier string
        plot_figure: <string> if "True" will plot normalised spectrum with fano
        out_path: <string> path to save, False by default
    Returns:
        results_dictionary: <dict> dictionary containing:
            popt: <array> fano fit parameters:
                peak, gamma, q, amplitude, (damping)
            pcov: <array> fano fit errors
                peak, gamma, q, amplitude, (damping)
        experimental_intensity: <array> experimental intensity array
        damping: <bool> if True, damping in fano equation
    '''
    normalised_intensity = normalise_intensity(intensity=intensity)
    try:
        popt4, pcov4 = opt.curve_fit(
            no_damping_fano,
            wavelength,
            normalised_intensity,
            initial_guesses[0: -1])
        errors4 = np.sqrt(np.diag(pcov4))
    except RuntimeError:
        print('\nRun Time Error \nNo Peak Found')
        popt4 = [0, 0, 0, 0]
        errors4 = [0, 0, 0, 0]
    except ValueError:
        print('\nValue Error \nNo Peak Found')
        popt4 = [0, 0, 0, 0]
        errors4 = [0, 0, 0, 0]
    damping = False
    if sum(popt4) == 0:
        try:
            popt5, pcov5 = opt.curve_fit(
                fano_resonance,
                wavelength,
                intensity,
                initial_guesses)
            errors5 = np.sqrt(np.diag(pcov5))
        except RuntimeError:
            print('\nRun Time Error \nNo Peak Found')
            popt5 = [0, 0, 0, 0, 0]
            errors5 = [0, 0, 0, 0, 0]
        except ValueError:
            print('\nValue Error \nNo Peak Found')
            popt5 = [0, 0, 0, 0, 0]
            errors5 = [0, 0, 0, 0, 0]
        optimized = popt5
        error = errors5
        peak_wavelength = optimized[0]
        peak_wavelength_error = error[0]
        if plot_figure == 'True':
            fanofitplot(
                wavelength=wavelength,
                intensity=intensity,
                fano=[
                    fano_resonance(
                        x=wav,
                        x0=optimized[0],
                        gamma=optimized[1],
                        q=optimized[2],
                        amplitude=optimized[3],
                        damping=optimized[4])
                    for wav in wavelength],
                peak_wavelength=peak_wavelength,
                peak_error=peak_wavelength_error,
                out_path=Path(f'{out_path}/{sample_name}_Experimental.png'))
        return_intensity = intensity
        damping = True
    else:
        optimized = popt4
        error = errors4
        peak_wavelength = optimized[0]
        peak_wavelength_error = error[0]
        if plot_figure == 'True':
            fanofitplot(
                wavelength=wavelength,
                intensity=intensity,
                fano=[
                    no_damping_fano(
                        x=wav,
                        x0=optimized[0],
                        gamma=optimized[1],
                        q=optimized[2],
                        amplitude=optimized[3])
                    for wav in wavelength],
                peak_wavelength=peak_wavelength,
                peak_error=peak_wavelength_error,
                out_path=Path(f'{out_path}/{sample_name}_Experimental.png'))
        return_intensity = normalised_intensity
    results_dictionary = {
        f'{sample_name} Fano Fit Parameters': [
            'Peak', 'Gamma', 'q', 'Ampltidude', 'Damping'],
        f'{sample_name} Fano Fit': [value for value in optimized],
        f'{sample_name} Fano Errors': [err for err in error]}
    return results_dictionary, return_intensity, damping


def get_S4_fano(wavelength,
                intensity,
                sample_name,
                experimental_fano,
                plot_figure,
                damping,
                out_path=False):
    '''
    Use scipy optimise and fano_resonance to iterate and find the peak in
    intensity and the corresponding wavelength and error value for an S4
    simulation.
    Args:
        wavelength: <array> wavelength array in nm
        intensity: <array> intensity array
        sample_name: <string> secondary sample identifier string
        experimental_fano: <array> experimental fano fit parameters
        plot_figure: <string> if "True" will plot normalised spectrum with fano
        out_path: <string> path to save, False by default
        damping: <bool> if True, includes damping in fano fit
    Returns:
        peak_wavelength: <dict> dictionary containing:
            popt: <array> fano fit parameters:
                peak, gamma, q, amplitude, (damping)
            pcov: <array> fano fit errors
                peak, gamma, q, amplitude, (damping)
    '''

    ''' Removed bounds '''

    if damping == 'True':
        #bounds = (
        #    (
        #        max(0, experimental_fano[0] - 50),
        #        max(-50, experimental_fano[1] - 5),
        #        max(-20, experimental_fano[2] - 5),
        #        max(0, experimental_fano[3] - 5),
        #        max(0, experimental_fano[4] - 5)),
        #    (
        #        min(1000, experimental_fano[0] + 50),
        #        min(50, experimental_fano[1] + 5),
        #        min(20, experimental_fano[2] + 5),
        #        min(5, experimental_fano[3] + 5),
        #        min(5, experimental_fano[4] + 5)))
        try:
            #popt, pcov = opt.curve_fit(
            #    fano_resonance,
            #    wavelength,
            #    intensity,
            #    experimental_fano,
            #    bounds=bounds)
            popt, pcov = opt.curve_fit(
                fano_resonance,
                wavelength,
                intensity,
                experimental_fano)
            errors = np.sqrt(np.diag(pcov))
        except RuntimeError:
            print('\nRun Time Error \nNo Peak Found')
            popt = [0, 0, 0, 0, 0]
            errors = [0, 0, 0, 0, 0]
        except ValueError:
            print('\nValue Error \nNo Peak Found')
            popt = [0, 0, 0, 0, 0]
            errors = [0, 0, 0, 0, 0]
        peak_wavelength = popt[0]
        peak_wavelength_error = errors[0]
        if plot_figure == 'True':
            fanofitplot(
            wavelength=wavelength,
            intensity=intensity,
            fano=[
                fano_resonance(
                    x=wav,
                    x0=popt[0],
                    gamma=popt[1],
                    q=popt[2],
                    amplitude=popt[3],
                    damping=popt[4])
                for wav in wavelength],
            peak_wavelength=peak_wavelength,
            peak_error=peak_wavelength_error,
            out_path=Path(f'{out_path}/{sample_name}_Simulation.png'))
    else:
        #bounds = (
        #    (
        #        max(0, experimental_fano[0] - 50),
        #        max(-50, experimental_fano[1] - 5),
        #        max(-20, experimental_fano[2] - 5),
        #        max(0, experimental_fano[3] - 5)),
        #    (
        #        min(1000, experimental_fano[0] + 50),
        #        min(50, experimental_fano[1] + 5),
        #        min(20, experimental_fano[2] + 5),
        #        min(5, experimental_fano[3] + 5)))
        try:
            #popt, pcov = opt.curve_fit(
            #    no_damping_fano,
            #    wavelength,
            #    intensity,
            #    experimental_fano,
            #    bounds=bounds)
            popt, pcov = opt.curve_fit(
                no_damping_fano,
                wavelength,
                intensity,
                experimental_fano)
            errors = np.sqrt(np.diag(pcov))
        except RuntimeError:
            print('\nRun Time Error \nNo Peak Found')
            popt = [0, 0, 0, 0]
            errors = [0, 0, 0, 0]
        except ValueError:
            print('\nValue Error \nNo Peak Found')
            popt = [0, 0, 0, 0]
            errors = [0, 0, 0, 0]
        peak_wavelength = popt[0]
        peak_wavelength_error = errors[0]
        if plot_figure == 'True':
            fanofitplot(
                wavelength=wavelength,
                intensity=intensity,
                fano=[
                    no_damping_fano(
                        x=wav,
                        x0=popt[0],
                        gamma=popt[1],
                        q=popt[2],
                        amplitude=popt[3])
                    for wav in wavelength],
                peak_wavelength=peak_wavelength,
                peak_error=peak_wavelength_error,
                out_path=Path(f'{out_path}/{sample_name}_Simulation.png'))
    return {
        f'{sample_name} S4 Fano Fit Parameters': [
            'Peak', 'Gamma', 'q', 'Ampltidude', 'Damping'],
        f'{sample_name} S4 Fano Fit': [value for value in popt],
        f'{sample_name} S4 Fano Errors': [value for value in errors]}


def compare_weighted_fano_parameter(experimental_fano,
                                    simulation_fano,
                                    weights):
    '''
    Compare experimental fano resonance to simulation fano resonance results.
    Use weighted analysis and the fano-resonance parameters to find a figure of
    merit, where ideal is 0. Parameters are weighted in accordance with relevant
    peak behvaiour. Fano parameters are ['Peak', 'Gamma', 'q', 'Ampltidude',
    'Damping'].
    Args:
        experimental_fano: <array> experimental fano parameters
        simulation_fano: <array> simulation fano parameters
        weights: <array> weights of importance of fano parameters
    Returns:
        figure_of_merit: <float> sum of weighted differences between the
                        experimental and simulation parameters
    '''
    ''' changed back to abs and removed squared factor '''
    #weighted_differences = [
    #    (weights[index] * np.abs((exp - simulation_fano[index]) / exp))
    #    for index, exp in enumerate(experimental_fano)]
    ''' Changed here to just differences '''
    weighted_differences = [
        ((exp - sim) ** 2)
        for exp, sim in zip(experimental_fano, simulation_fano)]
    figure_merit = sum(weighted_differences)
    return figure_merit, weighted_differences


def compare_fano_curves(experimental_intensity,
                        simulation_intensity):
    '''
    Compare experimental intensity and simulation intensity overlap integral.
    Args:
        wavelength: <array> wavelength at which to evaluate curves
        experimental_intensity: <array> experimental intensity curve
        simulation_intensity: <array> simulation intensity curve
    Returns:
        figure_merit: <float> 1 / overlap integral of two curves
    '''
    #overlap = [
    #    (experimental - simulation) ** 2
    #    for experimental, simulation in
    #    zip(experimental_intensity, simulation_intensity)]
    ''' Changed to just differences '''
    overlap = [
        np.abs(exp - sim) for exp, sim
        in zip(experimental_intensity, simulation_intensity)]
    figure_merit = sum(overlap)
    return figure_merit


def compare_experimental_simulation(sample_name,
                                    experimental_fano,
                                    experimental_fano_errors,
                                    simulation_fano,
                                    simulation_fano_errors,
                                    experimental_intensity,
                                    simulation_intensity,
                                    damping=False):
    '''
    Compare experimental and simulation fano fit curves and intensities to
    produce a figure of merit for optimization.
    Args:
        sample_name: <string> sample name identifier string
        experimental_fano: <array> experimental fano fit parameters
        experimental_fano_errors: <array> experimental fano fit errors
        simulation_fano: <array> simulation fano fit parameters
        simulation_fano_errors: <arrau> simulation fano fit errors
        experimental_intensity: <array> experimental intensity
        simulation_intensity: <array> simulation intensity
        damping: <bool> if True, damping included in fano
    Returns:
        figure_merit: <float> figure of merit for optimization
    '''
    if damping:
        #weights = [1.2, 0.9, 1.0, 1.1, 0.8]
        weights = [1.0, 1.0, 1.0, 1.0, 1.0]
    else:
        #weights = [1.2, 0.8, 0.9, 1.1]
        weights = [1.0, 1.0, 1.0, 1.0]
    if sum(simulation_fano) == 0:
        figure_merit = 500.0
        figure_merit_fano = simulation_fano
        figure_merit_integral = 0
        weighted = 0
    else:
        figure_merit_fano, weighted = compare_weighted_fano_parameter(
            experimental_fano=experimental_fano,
            simulation_fano=simulation_fano,
            weights=weights)
        figure_merit_integral = compare_fano_curves(
            experimental_intensity=experimental_intensity,
            simulation_intensity=simulation_intensity)
        #exp_fano_errors = sum([err ** 2 for err in experimental_fano_errors])
        #sim_fano_errors = sum([err ** 2 for err in simulation_fano_errors])
        figure_merit_values = [
            figure_merit_fano,
            figure_merit_integral]
        #figure_weights = [0.8, 1.2]
        figure_weights = [1.0, 0.0]
        ''' changed to no weights and removed integral '''
        figure_merit = sum([
            figure * weight
            for figure, weight
            in zip(figure_merit_values, figure_weights)])
    return {
        f'{sample_name} Figure Of Merit': figure_merit,
        f'{sample_name} Fano Error': figure_merit_fano,
        f'{sample_name} Overlap Error': figure_merit_integral,
        f'{sample_name} Weighted': weighted}


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
    parameters and intensity curves. Build argument string for S4 simulation,
    process S4 output, and calculate fano resonance parameters. Log iterations
    if desired.
    Args:
        variables: <array> simulation variables
        variable_names: <array> array of variable names
        constants: <dict> simulation constants dictionary
        sample_name: <string> sample name identifier string
        batch_name: <string> batch name identifier string
        peak_parameters: <dict> peak parameters dictionary from measurements
                        containing experimental wavelength, fano fit, and
                        errors
        lua_script: <string> lua script file
        plot_figure: <string> if "True" will plot simulation data
        out_path: <string> path to save results
        log_fom: <string> if "True" will log iteration data
    Returns:
        figure_of_merit: <float> iteration figure of merit
    '''

    ''' Changed here '''
    experimental_wavelength, experimental_intensity = simulation_wavelength_range(
        experimental_wavelength=peak_parameters[f'{sample_name} Wavelength'],
        experimental_intensity=peak_parameters[f'{sample_name} Intensity'])
    wavelength_range = [min(experimental_wavelength) - 50, max(experimental_wavelength) + 50, 1]
    ''' Done '''

    variables_dict = {
        'S4 Strings': [name for name in variables_names],
        'S4 Guesses': [guess for guess in variables]}
    argument_string = parameters_to_arguments(
        constants=constants,
        variables=variables_dict,
        wavelength_range=wavelength_range)
    S4_output = S4_RCWA(
        lua_script=lua_script,
        argument_string=argument_string)
    if constants['Reflection'] == 'True':
        #wavelength, _, sim_intensity = read_S4_output(
        #    process_string=S4_output)
        sim_wavelength, _, sim_intensity = read_S4_output(
            process_string=S4_output)
    else:
        #wavelength, sim_intensity, _ = read_S4_output(
        #    process_string=S4_output)
        sim_wavelength, sim_intensity, _ = read_S4_output(
            process_string=S4_output)
    #exp_fano, exp_intensity, damping = get_exp_fano_parameters(
    #    initial_guesses=peak_parameters[f'{sample_name} Fano Fit'],
    #    wavelength=wavelength,
    #    intensity=experimental_intensity,
    #    sample_name=sample_name,
    #    plot_figure=plot_figure,
    #    out_path=out_path)
    exp_fano, exp_intensity, damping = get_exp_fano_parameters(
        initial_guesses=peak_parameters[f'{sample_name} Fano Fit'],
        wavelength=experimental_wavelength,
        intensity=experimental_intensity,
        sample_name=sample_name,
        plot_figure=plot_figure,
        out_path=out_path)
    sim_fano = get_S4_fano(
        wavelength=sim_wavelength,
        intensity=sim_intensity,
        sample_name=sample_name,
        experimental_fano=exp_fano[f'{sample_name} Fano Fit'],
        plot_figure=plot_figure,
        out_path=out_path,
        damping=f'{damping}')

    ''' new '''
    adjusted_wavelengths_min = np.abs(
        sim_wavelength - min(peak_parameters[f'{sample_name} Wavelength']))
    adjusted_wavelengths_max = np.abs(
        sim_wavelength - max(peak_parameters[f'{sample_name} Wavelength']))
    wavelength_argmin = np.argmin(adjusted_wavelengths_min)
    wavelength_argmax = np.argmin(adjusted_wavelengths_max)
    trimmed_wavelengths = sim_wavelength[wavelength_argmin: wavelength_argmax + 1]
    trimmed_intensity = sim_intensity[wavelength_argmin: wavelength_argmax + 1]
    ''' finished '''

    #figure_of_merit = compare_experimental_simulation(
    #    sample_name=sample_name,
    #    experimental_fano=exp_fano[f'{sample_name} Fano Fit'],
    #    experimental_fano_errors=exp_fano[f'{sample_name} Fano Errors'],
    #    simulation_fano=sim_fano[f'{sample_name} S4 Fano Fit'],
    #    simulation_fano_errors=sim_fano[f'{sample_name} S4 Fano Errors'],
    #    experimental_intensity=exp_intensity,
    #    simulation_intensity=sim_intensity,
    #    damping=damping)
    figure_of_merit = compare_experimental_simulation(
        sample_name=sample_name,
        experimental_fano=exp_fano[f'{sample_name} Fano Fit'],
        experimental_fano_errors=exp_fano[f'{sample_name} Fano Errors'],
        simulation_fano=sim_fano[f'{sample_name} S4 Fano Fit'],
        simulation_fano_errors=sim_fano[f'{sample_name} S4 Fano Errors'],
        experimental_intensity=exp_intensity,
        simulation_intensity=trimmed_intensity,
        damping=damping)
    results_dict = dict(
        exp_fano,
        **sim_fano,
        **figure_of_merit)
    results_dict.update({f'{sample_name} Variables': variables_dict})
    results_dict.update({f'{sample_name} Constants': constants})
    results_dict.update({f'{sample_name} Wavelength Range': wavelength_range})
    results_dict.update({f'{sample_name} Damping': f'{damping}'})
    save_json_dicts(
        out_path=Path(f'{out_path}/Working_Result.json'),
        dictionary=results_dict)
    if log_fom == 'True':

        ''' Changed '''
        #log_figure_of_merit(
        #    variables=variables,
        #    figure_of_merit=figure_of_merit[f'{sample_name} Figure Of Merit'],
        #    simulation_fano=sim_fano[f'{sample_name} S4 Fano Fit'],
        #    experimental_fano=exp_fano[f'{sample_name} Fano Fit'],
        #    simulation_intensity=sim_intensity,
        #    experimental_intensity=exp_intensity,
        #    sample_name=sample_name,
        #    batch_name=f'{batch_name}',
        #    out_path=out_path)
        log_figure_of_merit(
            variables=variables,
            figure_of_merit=figure_of_merit[f'{sample_name} Figure Of Merit'],
            simulation_fano=sim_fano[f'{sample_name} S4 Fano Fit'],
            experimental_fano=exp_fano[f'{sample_name} Fano Fit'],
            simulation_intensity=trimmed_intensity,
            experimental_intensity=exp_intensity,
            sample_name=sample_name,
            batch_name=f'{batch_name}',
            out_path=out_path)
        ''' End '''

    return figure_of_merit[f'{sample_name} Figure Of Merit']


def optimize_simulation(sample_name,
                        batch_name,
                        variables,
                        constants,
                        peak_parameters,
                        lua_script,
                        plot_figure,
                        out_path,
                        log_fom):
    '''
    Optimize S4 fano parameters using scipy optimize minimize to vary simulation
    variables and calculate a minimum figure of merit with respect to the
    experimentally measured parameters. Return all parameters to one results
    dictionary.
    Args:
        sample_name: <string> sample name identifier string
        batch_name: <string> batch name identifier string
        variables: <dict> simulation variables dictionary
        constants: <dict> simulation constants dictionary
        peak_parameters: <dict> experimental peak dictionary
        lua_script: <string> lua script file
        plot_figure: <string> if "True" will plot simulation output
        out_path: <string> path to save results
        log_fom: <string> if "True" log iteration variables
        normalise: <bool> if True, intensity spectrum from experiment normalised
        normalise_sim: <bool> if True, intensity spectrim from sim normalised
    Returns:
        results_dictionary: <dict> grating results containing:
            variable results from optimizer
            simulation constants
            simulation fano results
    '''
    variable_guesses = variables['S4 Guesses']
    variable_names = variables['S4 Strings']
    variable_bounds = tuple([
        tuple(variable_pair)
        for variable_pair in variables['S4 Error Bounds']])
    print(variable_bounds)
    ftol = 1E-8
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
            'ftol': ftol,
            'gtol': 1E-8,
            'maxiter': 2000})
    optimizer_hessian = np.diag(optimizer_results.hess_inv.todense())
    optimizer_errors = [
        ftol * optimizer_results.fun * hess
        for hess in optimizer_hessian]
    simulation_path = Path(f'{out_path}/Working_Result.json')
    simulation_fano = load_json(file_path=simulation_path)
    variable_results = {
        f'{sample_name} Optimizer Results': [x for x in optimizer_results.x],
        f'{sample_name} Optimizer Errors': [x for x in optimizer_errors],
        f'{sample_name} Iterations': optimizer_results.nit}
    results_dict = dict(
        variable_results,
        **simulation_fano)
    os.remove(simulation_path)
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
    if arguments['All Arguments'] == 'True':
        grating_results = optimize_simulation(
            sample_name=f'{grating_name}',
            batch_name=f'{batch_name}',
            variables=variables,
            constants=constants,
            peak_parameters=peak_parameters,
            lua_script=lua_script,
            plot_figure=plot_figure,
            out_path=out_path,
            log_fom=log_fom)
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
                out_path=Path(f'{out_path}/{batch_name}_{grating_name}.png'),
                damping=grating_results[f'{grating_name} Damping'])
    else:
        grating_results = {
            f'{grating_name} Missing Parameters':
            [f'{argument}' for argument in arguments["Missing Arguments"]]}
        print(grating_results)
    return grating_results
