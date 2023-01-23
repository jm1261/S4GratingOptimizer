# S4 Grating Optimizer

This software is designed for analysis of resonant structures using rigorous coupled wave analysis (RCWA). Specifically the software uses Stanford Stratified Structure Solver (S<sup>4</sup>), a frequency domain code solving Maxwell's equations in layered periodic structures, to find any key parameter for a periodic structure. The software was designed with 1D resonant gratings in mind but can easily be adapted to include other resonant structures. Documentation for S<sup>4</sup> can be found [here](https://web.stanford.edu/group/fan/S4/index.html).

This software was created by Josh Male in December 2022 as part of the ongoing analysis of resonant structures. The code is under an MIT License and can be freely used.

## Table of Contents

* [General Information](#general-information)
* [Package Requirements](#package-requirements)
* [Launch](#launch)
* [Setup](#setup)
  * [Directory Paths](#directory-paths)
  * [File Names](#file-names)
    * [Film Thickness](#film-thickness)
    * [Grating Thickness](#grating-thickness)
    * [Resonance Position](#resonance-position)
    * [Grating Period](#grating-period)
  * [File Types](#file-types)
  * [Default info.json](#default-info)
* [General Work Flow](#general-work-flow)
  * [Organisation](#organisation)
  * [Data Handling](#data-handling)
  * [Parent Directory](#parent-directory)
  * [Batch Processing](#batch-processing)
  * [Find File paths](#find-file-paths)
* [Grating Optimizer](#grating-optimizer)
  * [Data Input](#data-input)
  * [Individual Parameter Optimization](#individual-parameter-optimization)
  * [Peak Analysis]
  * [Figure of Merit]
* [Acknowledgements]

## General Information

Automatic grating parameter optimizer built using Python 3, Lua, and S<sup>4</sup> for resonant structures. Suitable for all fano-resonance structures.

## Package Requirements

Language and package requirements can also be found in requirements.txt. The code was built using the following languages and module versions.

* Python 3.6
* tk 0.1.0
* numpy 1.21.4
* matplotlib 3.5.0
* scipy 1.7.2

## Launch

The code can be run from any terminal or editor. Main scripts are in the repository's main directory, while source code is stored safely in /src. The code relies on the use of info.json file for non-windows operating systems, this should be kept in the main repository directory.

## Setup

### Directory Paths

The code relies heavily on the info.json file for non-windows operating systems for data and results directory paths, this is to relieve pressure from any users from altering the script, and to ensure user alterations don't break the code. The info.json file should contain the following file paths:

* {
    "S4 Path": "/relative/path/to/data/files",
    "Results Path": "/relative/path/to/results/directory",
    "Plot Files": "True/False",
    "Log FOM": "True/False"
}

Where the relative paths are relative to the root directory (main directory of the repository). Windows devices make use of tkinter's interactive path selector and are able to select target files directly, but will require background and results paths to be present. Default paths are set to:

* {
    "S4 Path": "/S4",
    "Results Path": "/Results",
    "Plot Files": "False",
    "Log FOM": "True"
}

The code requires the data files directory to be called S4, i.e., files that intended for use by the S4 optimizer. No data can be processed if this directory is misnamed.

### File Names

Handling data and information is always a challenging aspect of coding. It's far easier to convey specific file information within the file name itself that to have a bank of dictionaries, text files, or other storage format files somewhere with all the relevant information required for data processing. This software is no different, and as such, there are specific pieces of information that are required from the file names that need to be included.

At least 1 measurement file is required for a grating parameter for the code to know what the identifier string for the process is called, though some modification could lead to a file-less system where the user inputs variables to optimize. The measurement files the code knows how to process are:

* Film thickness files
* Grating thickness files
* Resonance spectrum files
* Grating period files

These files are all .json files and they are each individual results files from the following links:

* [Film Thickness](https://github.com/jm1261/SurfaceProfileAnalysis)
* [Grating Thickness](https://github.com/jm1261/SurfaceProfileAnalysis)
* [Resonance Position](https://github.com/jm1261/PeakFinder)
* [Grating Period](https://github.com/jm1261/MicroscopyPeriodAnalysis)

Where the files require the following file name strings:

* 'Sample Identifier'_Film
* 'Sample Identifier'_Grating
* 'Sample Identifier'_Peak
* 'Sample Identifier'_Period

Respectively. The software knows how and what to do with these files and more could be added going forward. Providing at least 1 of these files is present in the /S4 folder, the code can run. There are checks in place to ensure that there are all the arguments required to run S4 and the code will run without all parameters with no output.

The same identifier string needs to be the same for all measured parameters for the same grating or gratings on the same chip. I.e., there can be multiple gratings on one sample identifier, and these should be indicated within the keys of the json files.

Assuming there are multiple gratings on each chip, let's look at what the individual files should look like. For this example, the chip name shall be A1, so all the files should be 'A1_Film.json', 'A1_Grating.json', etc. Keep in mind that any multiple grating identifier strings should be consistent over the separate measured parameter files, or the code will not know that the measurements are all for the same grating(s). The exception to this rule is the film thickness file that we shall cover first.

#### Film Thickness

The film thickenss file is by far the easiest to organise if you are not using the above linked software to process film thickness files. The film thickness file assumes that, because of resonant structures are on the same chip, then the film thickness must be the same for all structures. Therefore, the only keys required in the film thickness file for this software are:

* 'Average Result'
* 'Average Error'

In other words, the film thickness file only need contain a film thickness and an associated error key and value.

* {
    'Average Result': 150,
    'Average Error': 0.1
}

Where values must be given in nm. Film thickness can be given as a constant in the parameters file.

#### Grating Thickness

The grating thickness file, which is processed by the same software as the film thickness file, requires a list of all gratings on the chip under the key 'AFM Secondary String', and then a series of keys "(grating) Step Height" and "(grating) Step Height Error", where the (grating) is present in the original list of all gratings. For example:

* {
    'AFM Secondary String': [
        'P250',
        'P300',
        'P350'],
    'P250 Step Height': 150,
    'P300 Step Height': 156,
    'P350 Step Height': 156,
    'P250 Step Height Error': 1,
    'P300 Step Height Error': 2,
    'P350 Step Height Error': 3,
}

Where, again, the step heights and errors are in nm. Grating thicknesses can be given as a constant in the parameters file.

#### Resonance Position

The resonance position file is the single file that is required for the code to run, this file allows the code to cycle through individual gratings on chip and allows the code to find wavelength ranges, spectrum intensities, resonance peak target position, etc.

The files are so essential that we recommend using the following [software](https://github.com/jm1261/PeakFinder) to process experimentally obtained spectrum data and copy the output to the S4 directory here. The software is compatible with a wide range of spectrometers and spectrum capture software and should pose no problem for anyone wanting to use this software.

The README.md file there is easy to follow and guides the user through all attributes of using the software.

#### Grating Period

The grading period file, which is processed by the above linked software, but can also be manually entered, requires the following pieces of information. The software that processes measured grating periods uses SEM analysis, therefore requires a list of the gratings measured under 'SEM Secondary String', and then 'sample name Average' key, under which the grating key is listed to find the grating period. More information can be found on the above link, but grating periods can be given as a constant in the parameters file.

### File Types

All file types for measured parameters are stored in .json dictionary files. This is one of the preferred methods for storing and handling large amount of data in python and has been used for input data, output data, file paths, and input parameters.

### Default Info

As discussed in this section, all the user interface is done through the info.json file included in the main directory of the repository. This file can be opened using any text editor, and can be adjusted by the user with little-to-no consequence. The default info.json looks like this:

![info.json](./src/Images/info_dictionary.jpg)

## General Work Flow

This section contains general work flow details for code setup, data handling, batch processing, and a variety of other general features of the code.

### Organisation

The initial section of the code is concerned with finding directory and file paths required to process the data. This process is highly dependent on the operating system due to the relationship between tkinter and non-windows operating systems. The code makes use of pathlib's Path for file and directory paths, to maintain utility across operating systems, and sets the root directory as the main directory of the repository.

Details in the info.json setup file, including directory paths, are pulled in the function get_directory_paths, which builds the relative paths based on the set root path, which can be changed by adjusting the root variable. The directory paths is then treated as a dictionary containing the set directory paths by the user. File paths are returned as an array depending on operating system. If operating system is windows, the user can select which files they would like to process, if operating system is not windows, the user must have all desired files in the set directories as all files in those directories will be processed.

The first step of analysis is to make sure that all required variables are preset, if there are any missing files the code breaks the cycle and warns the user which variables are missing. Ensure that you have set the correct measured and set variables using S4_parameters.json file, and that the required measred variable files are present.

### Data Handling

Primary, secondary, and any other identifer strings are pulled from the file name string using sample_information function, which pulls all required parameters for the code from the file name and places them in a dictionary with clear identifier keys. The values of this dictionary are then used throughout the code to pull in required parameters. Therefore, this dictionary is referred to as sample_parameters throughout.

Measured parameters are pulled in using separate sample_information functions depending on the type of measured parameter in question. Secondary keys are only used here to distinguish between measured parameter types so that the code can process the data in the relevant way.

### Parent Directory

Parent directory is discussed above in directory paths. The code relies on relevant data being stored in the correctly name directory, allowing the code to distinguish between different data types. The parent directory names should be clear as discussed above. The parent directory is found using the function get_parent_directory, which returns a string of the parent directory name, which is then used as a parameter identifier key throughout.

### Batch Processing

For batch processing, indicated by the batch_ in either the script name or function name, the code then looks for like-file names, typically with the same primary identifier string, and groups them together for processing.

Batches are found using find_S4_batches function, which matches primary identifier strings, groups the file names, file paths, and measured parameter strings and returns a batch dictionary.

Using the batch keys and file paths stored within the batches dictionary, the code begins by pulling file names, file paths, and measured parameter strings into a batch results dictionary and appending each subsequent file parameter into an array under the appropriate key. The parent directory is used from here as a key identifier. The batch results utilises sample_information function in filepaths to pull in this information.

### Find File Paths

As discussed above, finding file paths is operating system dependent. On windows operating systems, the code uses tkinter's interactive file selection tool and allows the user to select any of the files in a directory they would like to process. In other operating systems, where tkinter is not so native, the code looks for all suitable files within the data directory and will process all of them, unless reslts have already been obtained and the results file exists.

## Grating Optimizer

This section contains work flow details to better understand how the grating optimizer is put together and what processes are undertaken in finding the optimum parameters for the periodic structure. It should be said, the code is not perfect and much more work could be put in to guarantee the success of the optimizer 100% of the time. Should the optimizer fail, we suggest a manual approach using S<sup>4</sup>. For the most part, the code is successful and helps to reduce processing time by narrowing in on the optimimum parameters.

### Data Input

Data is pulled in using the S4_args function in filepaths.py which reads through the S4_parameters.json file and sets the constants, parameters, and measured parameters as two dictionaries of iteration constants and iteration variables. I.e., which parameters are constant throughout optimization and which are to be optimized. The function also pulls in a list of argument names and checks against a known list of required arguments for the lua script to run. If any arguments are missing, the code will break out and return a statement to the user stating which arguments are missing. This is a fail safe as S<sup>4</sup> simulations will run but return 0s if there are any arguments missing.

### S<sup>4</sup> Simulation

S<sup>4</sup> can be controlled by lua or python scripts, where API documentation is detailed for both languages. Often it is easier to use python to feed arguments to the lua script. A large portion of the code is built around gathering all arguments for the simulation in python and feeding those arguments to the corresponding lua script.

### Individual Parameter Optimization

The optimizer begins by varying ....

# FINISH THIS LATER

Not sure I like how these are coming together and it might be better to explain the processes more generally..... i.e., this part does this and just have that as a standard response for all read me files. Rethink is probably necessary. Also good that you can have single # headers so that it looks like this. That might be better for separating sections later.