# bci solver - Overview

The purpose of the [bci solver](https://github.com/molmodcs/bci_solver) program is to predict the bci values [[1](#references),[2](#references),[3](#references),[4](#references)] present in a chemical structure from its partial atomic charges, as well as presenting this data in a clear fashion. It is divided into four modules, where the bci_solver_main module is the one where in fact the user will be able to interact with. It comes in two versions, a python script version which is usable in the command prompt line and a jupyter notebook one. 

In a nutshell, independently of the version preferred by the user, to utilize the program the user must provide the following data:

- Either a single .mol2 file corresponding to a chemical structure, or a folder containing .mol2 files which are organized in a suitable way which will be detailed in the [Generating and Organizing the Input Data section](#generating-and-organizing-the-input-data) in this README.

- The partial charges of the corresponding chemical structures represented by the .mol2 files. The format for these charge files must be either .xyz, an ORCA .log file format, or a sheet format (Pandas Dataframe,Excel,LibreOffice,CSV) following a standard outlined [Generating and Organizing the Input Data section](#generating-and-organizing-the-input-data).

The user can have a look at the [test-files](https://github.com/molmodcs/bci_solver/tree/main/test-files) folder to get a better understanding of the formats supported by the program.

As an output, the program will generate the following data:

- .par files containing for each .mol2 the bci values obtained through the calculations.
- Graphs and histograms depicting the variation between the bci values calculated across different structures or methods if computations are done across multiple chemical structures.

Explicit examples of output data are shown in the [Output](#output) section in this README and can also be acessed in the [Output]() folder available with the main project.

# Dependencies 

To properly run the script, it is first necessary to install the following dependencies;

Python == 3.12.0 - https://www.python.org/.

NumPy ==  1.26.0 - https://numpy.org/. 

~~~
pip install numpy
~~~

~~~
conda install numpy
~~~

SciPy == 1.12.0 - https://scipy.org/.

~~~
pip install scipy
~~~

~~~
conda install scipy
~~~

Matplotlib == 3.8.3 - https://matplotlib.org/.

~~~
pip install matplotlib
~~~

~~~
conda install matplotlib.
~~~

SciencePlots == 1.0.1 - https://pypi.org/project/SciencePlots/.

~~~
pip install SciencePlots.
~~~

Pandas == 2.2.1 - https://pandas.pydata.org/.

~~~
pip install pandas
~~~

odfpy == 1.4.1 - https://pypi.org/project/odfpy/.

~~~
pip install odfpy
~~~

~~~
conda install odfpy
~~~



# Modules

## bci_solver_main

In this section, we describe how to utilize the main module, namely bci_solver_main for calculating bci values of chemical structures.

### Generating and Organizing the Input Data

#### The .mol2 File Format 

These are the most important files and are required for all the calculations. A .mol2 file essentially provides the connectivity information about its corresponding chemical structure. They also provide the partial atomic charges of each atom consituting tbe structure, but not necessarily the computations need to be carried through using the charges available in the .mol2 file itself, they can also be provided by an external source, this will be further explained in the [Charge File Formats section](#charge-file-formats). 

<p align="center">
<img src="https://github.com/molmodcs/bci_solver/blob/main/Images/The%20.mol2%20File%20Format/mol2-file-example%20(cispla1-PBEQIDHsapoDZPdzp).png" alt="drawing",width = 250/>
<p align="center">
  
#### Charge File Formats

##### .xyz Format

<p align="center">
<img src="https://github.com/molmodcs/bci_solver/blob/main/Images/Charge%20File%20Formats/xyz%20Format/xyz-file-format-example%20(cispla1-PBEQIDHsapoDZPdzp-PBE-QIDH-Sapporo-DZP-2012-OPT).png" alt="drawing" width="250"/>
<p align="center">

##### Orca .log Format

<p align="center">
<img src="https://github.com/molmodcs/bci_solver/blob/main/Images/Charge%20File%20Formats/Orca%20.log%20Format/log-file-format-example%20(PLA01).png" alt="drawing" width="250"/>
<p align="center">

##### Sheet Format

<p align="center">
<img src="https://github.com/molmodcs/bci_solver/blob/main/Images/Charge%20File%20Formats/Sheet%20Format/sheet-file-format-example%20(chelpg-complexos-M06L-ZORA-TZVP).png" alt="drawing" width="250"/>
<p align="center">


### Single File Calculations

<p align="center">
<img src="https://github.com/molmodcs/bci_solver/blob/main/Images/Single%20File%20Calculations/mol2-file-example%20(cispla1-PBEQIDHsapoDZPdzp).png" alt="drawing" width="250"/>
<p align="center">

<p align="center">
<img src="https://github.com/molmodcs/bci_solver/blob/main/Images/Single%20File%20Calculations/xyz-file-format-example%20(cispla1-PBEQIDHsapoDZPdzp-PBE-QIDH-Sapporo-DZP-2012-OPT).png" alt="drawing" width="250"/>
<p align="center">

### Calculations with Multiple Files - Folder Calculations

#### Organizing the .mol2 Folder 

<p align="center">
<img src="https://github.com/molmodcs/bci_solver/blob/main/Images/Calculations%20with%20Multiple%20Files%20-%20Folder%20Calculations/Organizing%20the%20.mol2%20Folder/organizing-mol2-folder-example-solid-cisplatin.png" alt="drawing" width="250"/>
<p align="center">

#### Organizing the Charge Folder

<p align="center">
<img src="https://github.com/molmodcs/bci_solver/blob/main/Images/Calculations%20with%20Multiple%20Files%20-%20Folder%20Calculations/Organizing%20the%20Charge%20Folder/cisplatin%20-%20charge%20-%20folder%20-%20example.png" alt="drawing" width="250"/>
<p align="center">

<p align="center">
<img src="https://github.com/molmodcs/bci_solver/blob/main/Images/Calculations%20with%20Multiple%20Files%20-%20Folder%20Calculations/Organizing%20the%20Charge%20Folder/cisplatin%20-%20mol2%20-%20folder%20-example.png" alt="drawing" width="250"/>
<p align="center">

### User Input (bci_solver_main.py)

The user must provide first if computations will be done using a single .mol2 file or using a folder of .mol2 files by choosing between either -f/--file or -F/--folder. 

Then, the user must specify if the charges which will be used for the computations will be the ones available in the .mol2 files or if they pertain to an external source. This external source must either be a single file, in case the user has previously provived a -f/--file argument, which we recomend to have the same name as the .mol2 file, or in case the user has previously provided a -F/--folder argument, it must be a folder where each file has the same name as one of the files in the folder containing the chemical structures and has charges corresponding to the .mol2 file sharing its name. The available charge formats are either .xyz, an orca .log file format or a sheet format (Pandas Dataframe,Excel,LibreOffice) following the standard outlined in: [Generating and Organizing the Input Data section](#generating-and-organizing-the-input-data).

The user can also choose for the output folder. If one is not provided, the generated files and folders will be saved in a default output folder generated by th script. If the user provides an output folder which doesn't exist, the script will handle it by creating the folder with the specified path.

Finally, in case the user has provided a -F/--folder argument, he can optionally choose if the labels for the chemical structures appearing in the generated graphs and histograms will be the full corresponding file names by choosing --full or if the names will be shortened to first 6 characters of the file names (default behaviour). In any case, we recommend users to choose file names with a maximum length of 6 characters.

### List of Commands

* `-h`, `--help`: shows a help message and exits.
  
* `-o` OUTPUT, `--output` OUTPUT: specifies the output folder.
  
* `--full`: specifies that the labels appearing in the generated graphs and histograms will be the full .mol2 file names.

* `-f` FILE, `--file` FILE specifies that computations will be done using a single .mol2 file. Argument must be the path to a .mol2 file or the file name itself if it is available in the current directory.
  
* `-F` FOLDER, `--folder` FOLDER: specifies that computations will be done for multiple .mol2 files inside a folder. Argument must be the path to a folder containing .mol2 files or the folder name itself if it is available in the current directory.

* `-l` LOG, `--log` LOG: specifies that the atom charges for the respective .mol2 files are all available in a Orca .log format.
  
* `-s` SHEET, `--sheet` SHEET: specifies that the atom charges for the respective .mol2 files are all available in a sheet format.
  
* `-xyz` XYZ: specifies that the atom charges for the respective .mol2 files are all available in a .xyz format.

### Examples of Possible User Inputs (bci_solver_main.py)

* `-f '.mol2 file path'` 

<p align="center">
<img src="https://github.com/molmodcs/bci_solver/blob/main/Images/Examples%20of%20Possible%20User%20Inputs/user-input-file-example.png" alt="drawing" width="250"/>
<p align="center">

* `-f '.mol2 file path' -xyz 'path containing the charges which will be used for computations in .xyz format'`

<p align="center">
<img src="https://github.com/molmodcs/bci_solver/blob/main/Images/Examples%20of%20Possible%20User%20Inputs/user-input-file-xyz-example.png" alt="drawing" width="250"/>
<p align="center">

* `-F '.mol2 folder path' -s 'path containing the charges which will be used for computations in a sheet format'` 

<p align="center">
<img src="https://github.com/molmodcs/bci_solver/blob/main/Images/Examples%20of%20Possible%20User%20Inputs/user-input-folder-sheet-example.png" alt="drawing" width="250"/>
<p align="center">

* `-F '.mol2 folder path' -o 'output folder'`

<p align="center">
<img src="https://github.com/molmodcs/bci_solver/blob/main/Images/Examples%20of%20Possible%20User%20Inputs/user-input-folder-output-example.png" alt="drawing" width="250"/>
<p align="center">

* `-F 'mol22 folder path' --full` 

<p align="center">
<img src="https://github.com/molmodcs/bci_solver/blob/main/Images/Examples%20of%20Possible%20User%20Inputs/user-input-folder-full-example.png" alt="drawing" width="250"/>
<p align="center">

### User Input (bci_solver_main.ipynb)

### Output

<p align="center">
<img src="https://github.com/molmodcs/bci_solver/blob/main/Images/Output%20-%20Images/bci_values_N-H.png" alt="drawing" width="250"/>
<p align="center">

<p align="center">
<img src="https://github.com/molmodcs/bci_solver/blob/main/Images/Output%20-%20Images/bci_values_N-H_hist.png" alt="drawing" width="250"/>
<p align="center">

<p align="center">
<img src="https://github.com/molmodcs/bci_solver/blob/main/Images/Output%20-%20Images/bci_par_file.png" alt="drawing" width="250"/>
<p align="center">

<p align="center">
<img src="https://github.com/molmodcs/bci_solver/blob/main/Images/Output%20-%20Images/mol2_charge_conversion.png" alt="drawing" width="250"/>
<p align="center">

## bci_solver_mol2tools

## bci_solver_optimization

## bci_solver_visualization

# Technical Remarks

In the current version of the script, the user can not directly provide information about the atom types of the chemical structures. This is so because we still haven't figured out an efficient way the user can provide this information. So, in all the computations being done in this current version, there's no distinction between an atom and its corresponding atom type, i.e the same atom will always have the same atom type regardless of its role in the chemical structure. This might make some of the results slightly inaccurate. But future implementations of this program will resolve this issue.

# References

<b>[1]</b> Halgren, Thomas A. The Representation of van Der Waals (VdW) Interactions in Molecular Mechanics Force Fields: Potential Form, Combination Rules, and VdW Parameters. Sept. 1992, https://doi.org/10.1021/ja00046a032.

<b>[2]</b> Halgren, Thomas A. Merck Molecular Force Field. I. Basis, Form, Scope, Parameterization, and Performance of MMFF94. Apr. 1996, https://doi.org/10.1002/(sici)1096-987x(199604)17:5/6<490::aid-jcc1>3.0.co;2-p.

<b>[3]</b> Halgren, Thomas A.. “Merck molecular force field. II. MMFF94 van der Waals and electrostatic parameters for intermolecular interactions.” Journal of Computational Chemistry 17 (1996): n. pag.

<b>[4]</b> Halgren, Thomas A.. “Merck molecular force field. V. Extension of MMFF94 using experimental data, additional computational data, and empirical rules.” Journal of Computational Chemistry 17 (1996): n. pag.

