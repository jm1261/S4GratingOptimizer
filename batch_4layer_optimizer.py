import src.fileIO as io
import src.filepaths as fp
import src.analysis as anal

from pathlib import Path


def batch_4_layer_optimizer(batch_name,
                            file_paths,
                            parameters_path,
                            plot_figure,
                            out_path,
                            log_fom):
    '''
    Batch 4 layer grating optimizer.
    Args:
        batch_name: <string> batch name identifier string
        file_paths: <dict> measured file paths dictionary, can be null
        parameters_path: <string> path to S4_parameters.json file
        plot_figure: <string> If "True", will plot output files
        out_path: <string> path to save out figure or log files
        log_fom: <string> if "True" will log figure of merit optimization
    Returns:
        batch_dictionary: <dict> batch results dictionary
    '''
    peak_parameters, gratings = fp.get_experimental_gratings(
        file_path=file_paths['Peak Path'])
    batch_dictionary = {}
    all_paths = [file_paths[key] for key in file_paths.keys()]
    batch_dictionary.update({'File Paths': [f'{path}' for path in all_paths]})
    batch_dictionary.update({'Gratings': gratings})
    for grating in gratings:
        grating_results = anal.optimize_S4_grating(
            parameters_path=parameters_path,
            measured_paths=file_paths,
            batch_name=f'{batch_name}',
            grating_name=f'{grating}',
            peak_parameters=peak_parameters,
            lua_script='1D_4layer_grating.lua',
            plot_figure=plot_figure,
            out_path=out_path,
            log_fom=log_fom)
        batch_dictionary.update({f'{grating}': grating_results})
    return batch_dictionary


if __name__ == '__main__':
    ''' Organisation '''
    root = Path().absolute()
    info, directory_paths = fp.get_directory_paths(root_path=root)
    file_paths = fp.get_files_paths(
        directory_path=directory_paths['S4 Path'],
        file_string='.json')
    parent, batches = fp.get_S4_batches(file_paths=file_paths)

    ''' Loop Experimental Files '''
    for batch, filepaths in batches.items():
        out_file = Path(f'{directory_paths["Results Path"]}/{batch}_S4.json')
        if out_file.is_file():
            pass
        else:
            results_dictionary = batch_4_layer_optimizer(
                batch_name=batch,
                file_paths=filepaths,
                parameters_path=Path(f'{root}/S4_parameters.json'),
                plot_figure=info['Plot Files'],
                out_path=Path(f'{directory_paths["Results Path"]}'),
                log_fom=info['Log FOM'])
            io.save_json_dicts(
                out_path=out_file,
                dictionary=results_dictionary)
