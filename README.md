# bci solver

# Github Page

https://github.com/molmodcs/bci_solver

# Dependencies 

To properly run the script, it is first necessary to install the following dependencies;

Python == 3.12.0 - https://www.python.org/.

NumPy ==  1.26.0 - https://numpy.org/. 

~~~
pip install numpy - conda install numpy
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

# User Input

The user must provide first if computations will be done using a single .mol2 file or using a folder of .mol2 files by choosing between either -f/--file or -F/--folder. 

Then, the user must specify if the charges which will be used for the computations will be the ones available in the .mol2 files or if they pertain to an external source. This external source must either be a single file, in case the user has previously provived a -f/--file argument, which we recomend to have the same name as the .mol2 file, or in case the user has previously provided a -F/--folder argument, it must be a folder where each file has the same name as one of the files in the folder containing the chemical structures and has charges corresponding to the .mol2 file sharing its name. The available charge formats are either .xyz, an orca .log file format or a sheet format (Pandas Dataframe,Excel,LibreOffice) following the standard outlined in: https://github.com/molmodcs/bci_solver.

The user can also choose for the output folder. If one is not provided, the generated files and folders will be saved in a default output folder generated by th script. If the user provides an output folder which doesn't exist, the script will handle it by creating the folder with the specified path.

Finally, in case the user has provided a -F/--folder argument, he can optionally choose if the labels for the chemical structures appearing in the generated graphs and histograms will be the full corresponding file names by choosing --full or if the names will be shortened to first 6 characters of the file names (default behaviour). In any case, we recommend users to choose file names with a maximum length of 6 characters.

# Examples of possible user inputs

* -f '.mol2 file path' 
* -f '.mol2 file path' -xyz 'path containing the charges which will be used for computations in .xyz format'
* -F '.mol2 folder path' -s 'path containing the charges which will be used for computations in a sheet format'
* -F '.mol2 folder path' -o 'output folder'
* -F 'mol22 folder path' --full

# Technical Remarks

In the current version of the script, the user can not directly provide information about the atom types of the chemical structures. This is so because we still haven't figured out an efficient way the user can provide this information. So, in all the computations being done in this current version, there's no distinction between an atom and its corresponding atom type, i.e the same atom will always have the same atom type regardless of its role in the chemical structure. This might make some of the results slightly inaccurate. But future implementations of this program will resolve this issue.
