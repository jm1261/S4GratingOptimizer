import numpy as np
import matplotlib.pyplot as plt

from pathlib import Path


def fano_resonance(x, x0, gamma, q, amplitude):
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


def fanofitplot(wavelength,
                intensity,
                fano,
                out_path,
                peak_wavelength,
                peak_error,
                show=False):
    '''
    Plot wavelength, intensity, and fano peak fit on same axis.
    Args:
        wavelength: <array> wavelength array
        intensity: <array> intensity array
        fano: <array> fano fit intensity array
        out_path: <string> path to save
        peak_wavelength: <string> calculated peak wavelength
        peak_error: <string> calculated peak wavelength error
        show: <bool> if True, plot shows, always saves
    Returns:
        None
    '''
    fig, ax = plt.subplots(
        nrows=1,
        ncols=1,
        figsize=[10, 7])
    ax.plot(
        wavelength,
        intensity,
        'b',
        lw=2,
        label='Data')
    ax.plot(
        wavelength,
        fano,
        'r',
        lw=2,
        label='Fano Fit')
    ax.legend(
        frameon=True,
        loc=0,
        prop={'size': 14})
    ax.set_xlabel(
        'Wavelength [nm]',
        fontsize=14,
        fontweight='bold')
    ax.set_ylabel(
        'Intensity [au]',
        fontsize=14,
        fontweight='bold')
    ax.tick_params(
        axis='both',
        colors='black',
        labelsize=14)
    text_string = (
        f'Peak = ({round(peak_wavelength, 2)} +/- {round(peak_error, 2)})nm')
    props = dict(
        boxstyle='round',
        facecolor='wheat',
        alpha=0.5)
    ax.text(
        0.05,
        0.05,
        text_string,
        transform=ax.transAxes,
        verticalalignment='top',
        bbox=props)
    if show:
        plt.show()
    plt.savefig(out_path)
    fig.clf()
    plt.cla()
    plt.close(fig)


def fomplot(wavelengths,
            variables,
            figure_of_merit,
            fano_parameters,
            simulation_fano,
            experimental_fano,
            simulation_intensity,
            experimental_intensity,
            variable_names,
            out_path,
            show=False):
    '''
    Plot figure of merit log and simulation output comparison from S4 optimizer.
    Args:
        wavelengths: <array> [start, stop, step]
        variables: <array> nested array of simulation/optimizer variables for
                    iterations
        figure_of_merit: <array> figure of merit for iterations
        fano_parameters: <array> nested array of optimzied fano parameters for
                            iterations
        simulation_fano: <array> simulation best fano parameters
        experimental_fano: <array> experimental best fano parameters
        simulation_intensity: <array> simulation intensity for best parameters
        experimental_intensity: <array> experimental intensity for best
                                parameters
        variable_names: <array> variable name strings for variables array, must
                        be same length
        out_path: <string> path to save
        show: <bool> if True, graph show, always saves
    Returns:
        None
    '''
    fano_names = ['x0', 'gamma', 'q', 'amplitude']
    wavelength_range = np.arange(
        wavelengths[0],
        wavelengths[1] + (wavelengths[2] / 10),
        wavelengths[2])
    sim_fano = [
        fano_resonance(
            x=x,
            x0=simulation_fano[0],
            gamma=simulation_fano[1],
            q=simulation_fano[2],
            amplitude=simulation_fano[3])
        for x in wavelength_range]
    exp_fano = [
        fano_resonance(
            x=x,
            x0=experimental_fano[0],
            gamma=experimental_fano[1],
            q=experimental_fano[2],
            amplitude=experimental_fano[3])
        for x in wavelength_range]
    x_values = range(len(figure_of_merit))
    fig, (ax1, ax2, ax3) = plt.subplots(
        nrows=3,
        ncols=1,
        figsize=[20, 14])
    line0, = ax1.plot(
        x_values,
        figure_of_merit,
        'k',
        lw=2,
        label='Figure of Merit')
    ax1.set_xlabel(
        'Iterations [au]',
        fontsize=14,
        fontweight='bold')
    ax1.yaxis.label.set_color(line0.get_color())
    ax1.set_xlim(0, len(figure_of_merit))
    ax1.set_ylabel(
        'Figure of Merit [log]',
        fontsize=14,
        fontweight='bold')
    ax1.tick_params(
        axis='y',
        colors=line0.get_color(),
        size=12)
    ax1.set_yscale('log')
    line1, = ax2.plot(
        x_values,
        figure_of_merit,
        'k',
        lw=2,
        label='Figure of Merit')
    ax2.set_xlabel(
        'Iterations [au]',
        fontsize=14,
        fontweight='bold')
    ax2.yaxis.label.set_color(line1.get_color())
    ax2.tick_params(
        axis='y',
        colors=line1.get_color(),
        size=12)
    ax2.set_xlim(0, len(figure_of_merit))
    ax2.set_ylabel(
        'Figure of Merit [log]',
        fontsize=14,
        fontweight='bold')
    ax2.set_yscale('log')
    for index, variable in enumerate(variables):
        twinax = ax1.twinx()
        twinax.spines.right.set_position(('axes', 1 + (index / 10)))
        twinline, = twinax.plot(
            x_values,
            variable,
            f'C{index}',
            lw=2,
            label=f'{variable_names[index]}')
        twinax.set_ylabel(
            f'{variable_names[index]}',
            fontsize=14,
            fontweight='bold')
        twinax.yaxis.label.set_color(twinline.get_color())
        twinax.tick_params(
            axis='y',
            colors=twinline.get_color(),
            size=12)
    ax1.legend(
        loc=0,
        prop={'size': 14})
    for index, variable in enumerate(fano_parameters):
        twinax = ax2.twinx()
        twinax.spines.right.set_position(('axes', 1 + (index / 10)))
        twinline, = twinax.plot(
            x_values,
            variable,
            f'C{index}',
            lw=2,
            label=f'{fano_names[index]}')
        twinax.set_ylabel(
            f'{fano_names[index]}',
            fontsize=14,
            fontweight='bold')
        twinax.yaxis.label.set_color(twinline.get_color())
        twinax.tick_params(
            axis='y',
            colors=twinline.get_color(),
            size=12)
    ax2.legend(
        loc=0,
        prop={'size': 14})
    ax3.plot(
        wavelength_range,
        sim_fano,
        'r',
        lw=2,
        label='Simulated Fano')
    ax3.plot(
        wavelength_range,
        exp_fano,
        'b',
        lw=2,
        label='Experimental Fano')
    ax3.set_xlabel(
        'Wavelength [nm]',
        fontsize=14,
        fontweight='bold')
    ax3.set_ylabel(
        'Intensity [au]',
        fontsize=14,
        fontweight='bold')
    ax3.tick_params(
        axis='both',
        size=12)
    ax3.plot(
        wavelength_range,
        simulation_intensity,
        'g',
        lw=2,
        label='Simulation Intensity')
    ax3.plot(
        wavelength_range,
        experimental_intensity,
        'k',
        lw=2,
        label='Experimental Intensity')
    ax3.legend(
        loc=0,
        prop={'size': 14})
    fig.tight_layout()
    if show:
        plt.show()
    plt.savefig(out_path)
    plt.cla()
    fig.clf()
    plt.close(fig)


def S4_plot(wavelength,
            simulation_intensity,
            experimental_intensity,
            batch_name,
            sample_name,
            out_path):
    '''
    Plot S4 wavelength and intensity.
    Args:
        wavelength: <array> wavelength array in nm
        simulation_intensity: <array> simulation intensity array in au
        experimental_intensity: <array> experimental intensity array in au
        batch_name: <string> batch name identifier string
        sample_name: <string> sample name identifier string
        out_path: <string> path to save
    Returns:
        None
    '''
    fig, ax = plt.subplots(
        nrows=1,
        ncols=1,
        figsize=[10, 7])
    ax.plot(
        wavelength,
        simulation_intensity,
        'b',
        lw=2,
        label='simulation')
    ax.plot(
        wavelength,
        experimental_intensity,
        'r',
        lw=2,
        label='experimental')
    ax.legend(
        loc=0,
        prop={'size': 14})
    ax.set_xlabel(
        'Wavelength [nm]',
        fontsize=14,
        fontweight='bold')
    ax.set_ylabel(
        'Intensity [au]',
        fontsize=14,
        fontweight='bold')
    ax.tick_params(
        axis='both',
        size=12)
    ax.set_title(
        f'{sample_name}',
        fontsize=14,
        fontweight='bold')
    plt.savefig(Path(f'{out_path}/{batch_name}_{sample_name}_S4plot.png'))
    plt.cla()
    fig.clf()
    plt.close(fig)
