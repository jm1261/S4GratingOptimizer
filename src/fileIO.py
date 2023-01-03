import json
import numpy as np

from pathlib import Path


def load_json(file_path):
    '''
    Extract user variables from json dictionary.
    Args:
        file_path: <string> path to file
    Returns:
        dictionary: <dict> use variables dictionary
    '''
    with open(file_path, 'r') as file:
        return json.load(file)


def convert(o):
    '''
    Check type of data string
    '''
    if isinstance(o, np.generic):
        return o.item()
    raise TypeError


def save_json_dicts(out_path,
                    dictionary):
    '''
    Save dictionary to json file.
    Args:
        out_path: <string> path to file, including file name and extension
        dictionary: <dict> python dictionary to save out
    Returns:
        None
    '''
    with open(out_path, 'w') as outfile:
        json.dump(
            dictionary,
            outfile,
            indent=2,
            default=convert)
        outfile.write('\n')


def load_S4_parameters(file_path):
    '''
    Load S4 simulation parameters from S4_parameters.json file in main directory
    of code repository. Without this file, code will not run. Allows user to
    give the lua strings and values of all simulation constants, lua strings,
    guesses, and error bounds of simulation variables, and lua strings and file
    extension strings of all measured variables or constants.
    Args:
        file_path: <string> path to S4 parameters json file
    Returns:
        constants: <dict> constants dictionary containing reflection bool,
                    constants lua name strings, constants values
        variables: <dict> variables dictionary containing lua name strings,
                    initial guesses, error bounds
        measured: <dict> measured variables dictionary containing constant bool,
                    lua name strings, file name strings
    '''
    simulation_parameters = load_json(file_path=file_path)
    constants = simulation_parameters['Constants']
    variables = simulation_parameters['Variables']
    measured = simulation_parameters['Measured']
    return constants, variables, measured


def get_grating_periods(file_path,
                        batch_name,
                        grating_name):
    '''
    Pull grating period and errors from period json result file. Returns 3x the
    measured error (i.e., 3 sigma error range).
    Args:
        file_path: <string> path to file
        batch_name: <string> batch name identifier string, e.g. A1
        grating_name: <string> grating/second_string identifier key
    Returns:
        period: <float> measured grating period
        period_error <float> measured grating period error
    '''
    period_parameters = load_json(file_path=file_path)
    period, error = (
        (period_parameters[f'{batch_name} Average'])[f'{grating_name}'])
    period_error = 10 * error
    return period, period_error


def get_grating_thicknesses(file_path,
                            grating_name):
    '''
    Pull grating thicknesses and errors from grating_thickness json result file.
    Returns 3x the measured error (i.e., 3 sigma error range).
    Args:
        file_path: <string> path to file
    Returns:
        grating_thickness: <float> measured grating thickness
        grating_error: <float> measured grating thickness error
    '''
    grating_parameters = load_json(file_path=file_path)
    grating_thickness = grating_parameters[f'{grating_name} Step Height']
    error = grating_parameters[f'{grating_name} Step Height Error']
    grating_error = 10 * error
    return grating_thickness, grating_error


def get_film_thickness(file_path):
    '''
    Pull average film thickness and error from film_thickness json result file.
    Returns 3x the measured error (i.e., 3 sigma error range).
    Args:
        file_path: <string> path to file
    Returns:
        film_thickness: <float> measured film thickness
        film_error: <float> measured film thickness error
    '''
    film_parameters = load_json(file_path=file_path)
    film_thickness = film_parameters['Average Result']
    error = film_parameters['Average Error']
    film_error = 10 * error
    return film_thickness, film_error


def get_peak_wavelength(file_path,
                        grating_name):
    '''
    Gets target peak wavelength and error from measured experimental peak json.
    Returns 3x the measured error (i.e., 3 sigma error range).
    Args:
        file_path: <string> path to file
        grating_name: <string> grating string identifier
    Returns:
        peak: <float> measured experimental peak wavelength
        peak_error: <float> measured experimental peak wavelength error
    '''
    peak_parameters = load_json(file_path=file_path)
    peak = peak_parameters[f'{grating_name} Peak Wavelength']
    error = peak_parameters[f'{grating_name} Peak Error']
    peak_error = 3 * error
    return peak, peak_error


def measured_S4_parameter_files(file_paths,
                                batch_name,
                                parameter_string,
                                grating_name):
    '''
    Pull in measured grating parameters such as grating period, grating
    thickness, and film thickness measured and stored in json files.
    Args:
        file_paths: <dict> file paths to measured parameter files dictionary
        batch_name: <string> batch name identifier string
        parameter_string: <string> measured parameter string
        grating_name: <string> grating name identifier string
    Returns:
        value: <float> measured parameter value
        error_bound: <array> measured parameter value error bound [lower, upper]
    '''
    file_path = file_paths[f'{parameter_string} Path']
    if parameter_string == 'Period':
        value, error = get_grating_periods(
            file_path=file_path,
            batch_name=batch_name,
            grating_name=(grating_name.split('_'))[0])
        error_bound = [value - error, value + error]
    elif parameter_string == 'Peak':
        value, error = get_peak_wavelength(
            file_path=file_path,
            grating_name=grating_name)
        error_bound = [value - error, value + error]
    elif parameter_string == 'Grating':
        value, error = get_grating_thicknesses(
            file_path=file_path,
            grating_name=(grating_name.split('_'))[0])
        film, film_error = get_film_thickness(file_path=file_paths['Film Path'])
        error_bound = [value - error, film]
    elif parameter_string == 'Film':
        value, error = get_film_thickness(file_path=file_path)
        error_bound = [value - error, value + error]
    return value, error_bound


def get_S4_polarisation(grating_string):
    '''
    Pull grating polarisation (TE/TM) from file string identifier and return
    polarisation array for S4 processing. Will default to TE operation if string
    is not explicit.
    Args:
        grating_string: <string> grating string from file keys, must be in the
                        format name_TE or name_TM
    Returns:
        string: <array> lua strings for polarisaation
        polarisation: <array> [1, 0] or [0, 1] for TE/TM operation respectively
    '''
    split = grating_string.split('_')
    string = ['TE', 'TM']
    if split[1] == 'TM':
        polarisation = [0, 1]
    else:
        polarisation = [1, 0]
    return string, polarisation


def check_4layer_parameters(constants,
                            variables):
    '''
    Check to see if all required lua argument names for 1D_4layer_grating.lua
    are present in supplied constants and variables, code will not progress
    if arguments are missing.
    Args:
        constants: <dict> simulation constants dictionary
        variables: <dict> simulation variables dictionary
    Returns:
        all_keys[0]: <bool> "True" or "False" depending on all present keys
        missing_arguments: <array> array of missing argument names
    '''
    arguments = [
        'harmonics', 'cover_n', 'cover_k', 'substrate_n', 'substrate_k',
        'material_n', 'material_k', 'fill_factor', 'period',
        'grating_thickness', 'film_thickness', 'TE', 'TM',
        'peak_wavelength']
    key_present = []
    user_arguments = []
    missing_arguments = []
    for argument in constants['S4 Strings']:
        user_arguments.append(argument)
    for argument in variables['S4 Strings']:
        user_arguments.append(argument)
    for argument in arguments:
        if argument in user_arguments:
            key_present.append('True')
        else:
            key_present.append('False')
            missing_arguments.append(argument)
    all_keys = ['False' if 'False' in key_present else 'True']
    return {
        "All Arguments": all_keys[0],
        "Missing Arguments": [missing for missing in missing_arguments]}


def S4_args(parameters_path,
            file_paths,
            batch_name,
            grating_name):
    '''
    Pull grating specific simulation arguments, constants and variables for
    a specific grating.
    Args:
        parameters_path: <string> path to parameters info json
        file_paths: <dict> file paths to measured parameter files dictionary
                    can be null if no measured parameters
        batch_name: <string> batch name identifier string
        grating_name: <string> grating name identifier string
    Returns:
        arguments: <dict> arguments dictionary containing true/false and missing
                    lua arguments
        constants: <dict> constants dictionary containing lua names and values
        variables: <dict> variables dictionary containing lua names and values
    '''
    constants, variables, measured = load_S4_parameters(
        file_path=parameters_path)
    m_constants, m_strings, m_files, m_precision = measured.items()
    for index, statement in enumerate(m_constants[1]):
        try:
            value, error = measured_S4_parameter_files(
                file_paths=file_paths,
                batch_name=batch_name,
                grating_name=grating_name,
                parameter_string=(m_files[1])[index])
            if np.isnan(value):
                pass
            else:
                if statement == 'True':
                    constants['S4 Strings'].append((m_strings[1])[index])
                    constants['S4 Values'].append(value)
                else:
                    variables['S4 Strings'].append((m_strings[1])[index])
                    variables['S4 Guesses'].append(value)
                    variables['S4 Error Bounds'].append(error)
                    variables['S4 Precision'].append((m_precision[1])[index])
        except:
            pass
    polarise_string, polarise_value = get_S4_polarisation(
        grating_string=grating_name)
    for string, value in zip(polarise_string, polarise_value):
        constants['S4 Strings'].append(string)
        constants['S4 Values'].append(value)
    arguments = check_4layer_parameters(
        constants=constants,
        variables=variables)
    return arguments, constants, variables


def log_figure_of_merit(variables,
                        figure_of_merit,
                        simulation_fano,
                        experimental_fano,
                        simulation_intensity,
                        experimental_intensity,
                        sample_name,
                        batch_name,
                        out_path):
    '''
    Log figure of merit parameters for plotting and logging data.
    Args:
        variables: <array> simulation iteration variables
        figure_of_merit: <float> figure of merit for iteration
        simulation_fano: <array> simulation fano parameters for iteration
        experimental_fano: <array> experimental fano parameters
        sample_name: <string> sample name identifier string
        batch_name: <string> batch name identifier string
        out_path: <string> path to save
    Returns:
        None
    '''
    variable_line = [variable for variable in variables]
    variable_line.append(figure_of_merit)
    [variable_line.append(variable) for variable in simulation_fano]
    variable_string = ",".join(map(str, variable_line))
    fano_data = (simulation_fano, experimental_fano)
    intensity_data = (simulation_intensity, experimental_intensity)
    outpath = Path(f'{out_path}/{batch_name}_{sample_name}_FOM.csv')
    with open(outpath, 'a') as outfile:
        outfile.write(f'{variable_string}\n')
    outpath = Path(f'{out_path}/{batch_name}_{sample_name}_Fano.csv')
    with open(outpath, 'w') as outfile:
        [
            outfile.write(",".join(map(str, fano)) + '\n')
            for fano in np.array(fano_data).T]
    outpath = Path(f'{out_path}/{batch_name}_{sample_name}_Intensity.csv')
    with open(outpath, 'w') as outfile:
        [
            outfile.write(",".join(map(str, intensity)) + '\n')
            for intensity in np.array(intensity_data).T]


def load_fomlog(directory_path,
                sample_name,
                batch_name,
                length_variables):
    '''
    Load figure of merit logs for plotting.
    Args:
        directory_path: <string> path to figure of merit logs directory
        sample_name: <string> sample name identifier (same as save out)
        batch_name: <string> batch name identifier string
        length_variables: <int> length of iteration variables
    Returns:
        varbs: <2D array> array of iteration variables over iteration loops
        fom: <array> figure of merit over iteration loops
        fano_parameters: <2D array> array of simulation fano parameters over
                        iteration loops
        sim_fano: <array> best figure of merit simulation fano parameters
        exp_fano: <array> experimental fano parameters
        sim_int: <array> simulation intensity array
        exp_int: <array> experimental intensity array
    '''
    data = np.genfromtxt(
        fname=Path(f'{directory_path}/{batch_name}_{sample_name}_FOM.csv'),
        delimiter=',')
    varbs = [data[:, i] for i in range(0, length_variables, 1)]
    fom = data[:, length_variables]
    fano_params = [
        data[:, i] for i in range(length_variables + 1, len(data[0, :]))]
    sim_fano, exp_fano = np.genfromtxt(
        fname=Path(f'{directory_path}/{batch_name}_{sample_name}_Fano.csv'),
        delimiter=',',
        unpack=True)
    sim_int, exp_int = np.genfromtxt(
        fname=Path(
            f'{directory_path}/{batch_name}_{sample_name}_Intensity.csv'),
        delimiter=',',
        unpack=True)
    return varbs, fom, fano_params, sim_fano, exp_fano, sim_int, exp_int
