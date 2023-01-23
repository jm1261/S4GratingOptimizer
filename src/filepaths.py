import os

from pathlib import Path
from sys import platform
from src.fileIO import load_json
from src.GUI import prompt_for_path


def check_platform():
    '''
    Check operating system.
    Args:
        None
    Returns:
        operating_system: <string> "Windows", "Linux", or "Mac"
    '''
    if platform == 'linux' or platform == 'linux2':
        operating_system = 'Linux'
    elif platform == 'darwin':
        operating_system = 'Mac'
    elif platform == 'win32':
        operating_system = 'Windows'
    return operating_system


def get_directory_paths(root_path):
    '''
    Get target data path and results path from info dictionary file.
    Args:
        root_path: <string> path to root directory
    Returns:
        data_path: <string> path to data directory
        bg_path: <string> path to background directory
        results_path: <string> path to results directory
        info: <dict> information dictionary (info.json)
    '''
    info = load_json(file_path=Path(f'{root_path}/info.json'))
    directory_paths = {}
    for key, value in info.items():
        if 'Path' in key:
            directory_paths.update({key: Path(f'{root_path}/{value}')})
    return info, directory_paths


def extractfile(directory_path,
                file_string):
    '''
    Pull file from directory path.
    Args:
        directory_path: <string> path to file
        file_string: <string> string contained within file name
    Returns:
        array: <array> array of selected files
    '''
    directory_list = sorted(os.listdir(directory_path))
    return [file for file in directory_list if file_string in file]


def get_files_paths(directory_path,
                    file_string):
    '''
    Get target files and directory paths depending on the operating system.
    Args:
        directory_path: <string> path to data directory
        file_string: <string> file extension (e.g. .csv)
    Returns:
        file_paths: <string> path to files
    '''
    operating_system = check_platform()
    if operating_system == 'Linux' or operating_system == 'Mac':
        file_list = extractfile(
            directory_path=directory_path,
            file_string=file_string)
        file_paths = [Path(f'{directory_path}/{file}') for file in file_list]
    elif operating_system == 'Windows':
        file_paths = prompt_for_path(
            default=directory_path,
            title='Select Target File(s)',
            file_path=True,
            file_type=[(f'{file_string}', f'*{file_string}')])
    return file_paths


def get_parent_directory(file_path):
    '''
    Find parent directory name of target file.
    Args:
        file_path: <string> path to file
    Returns:
        parent_directory: <string> parent directory name (not path)
    '''
    operating_system = check_platform()
    dirpath = os.path.dirname(file_path)
    if operating_system == 'Linux':
        dirpathsplit = dirpath.split('/')
    else:
        dirpathsplit = dirpath.split('\\')
    parent_directory = dirpathsplit[-1]
    return parent_directory


def get_filename(file_path):
    '''
    Splits file path to remove directory path and file extensions.
    Args:
        file_path: <string> path to file
    Returns:
        file_name: <string> file name without path or extensions
    '''
    return os.path.splitext(os.path.basename(file_path))[0]


def S4_sample_information(file_path):
    '''
    Pull sample parameters for file name string for various processes.
    Args:
        file_path: <string> path to file
    Returns:
        sample_parameters: <dict>
    '''
    parent_directory = get_parent_directory(file_path=file_path)
    file_name = get_filename(file_path=file_path)
    file_split = file_name.split('_')
    return {
        "Parent Directory": parent_directory,
        f'{parent_directory} File Name': file_name,
        f'{parent_directory} File Path': f'{file_path}',
        f'{parent_directory} Primary String': file_split[0],
        f'{parent_directory} Secondary String': file_split[1]}


def sample_information(file_path):
    '''
    Pull sample parameters based on which type of file is being analysed.
    Args:
        file_path: <string> path to file
    Returns:
        sample_parameters: <dict>
    '''
    parent_directory = get_parent_directory(file_path=file_path)
    if parent_directory == 'S4':
        sample_parameters = S4_sample_information(file_path=file_path)
    else:
        sample_parameters = {}
    return sample_parameters


def get_S4_batches(file_paths):
    '''
    Find all sample batches in series of file paths and append file paths to
    batch names for loop processing.
    Args:
        file_paths: <array> array of target file paths
    Returns:
        parent: <string> parent directory string
        batches: <dict>
            Batch indicators: <dict> file paths for film thickness, grating
                                thickness, peak wavelength, period measurement
    '''
    batches = {}
    for file in file_paths:
        sample_parameters = sample_information(file_path=file)
        parent = sample_parameters['Parent Directory']
        primary_key = f'{parent} Primary String'
        secondary_key = f'{parent} Secondary String'
        file_type = {f'{sample_parameters[secondary_key]} Path': file}
        if sample_parameters[primary_key] in batches.keys():
            batches[f'{sample_parameters[primary_key]}'].update(file_type)
        else:
            batches.update({f'{sample_parameters[primary_key]}': file_type})
    return parent, batches


def get_experimental_gratings(file_path):
    '''
    Get all batch gratings to processed irrespective of polarisation.
    Args:
        file_path: <string> path to peak parameters json
    Returns:
        peak_parameters: <dict> peak parameters json
        gratings: <array> array of grating strings
    '''
    peak_parameters = load_json(file_path=file_path)
    all_peaks = [p for p in peak_parameters['Spectrum Secondary String']]
    all_gratings = [(p.split('_'))[0] for p in all_peaks]
    gratings = []
    for grating in all_gratings:
        if grating in gratings:
            pass
        else:
            gratings.append(grating)
    return peak_parameters, gratings
