import numpy as np
import pandas as pd
import os
from itertools import chain

#From the absolute file path of a mol2 file, returns for the file the adjacency matrix corresponding to the associated chemical structure.
def convert(mol2_file_path):    
    with open(mol2_file_path, 'r') as f:
        lines = f.readlines()
        word1 = '@<TRIPOS>MOLECULE'
        word2 = '@<TRIPOS>BOND'
        index_for_number_of_atoms = 0
        number_of_atoms = 0
        number_of_bonds = 0
        index_for_bonds = 0
        for line in lines:
            if line.find(word1) != -1:
                index_for_number_of_atoms = lines.index(line) + 2
                number_of_atoms = int(lines[index_for_number_of_atoms].split()[0])
                number_of_bonds = int(lines[index_for_number_of_atoms].split()[1])
                break
        initial_array = np.zeros((number_of_atoms,number_of_atoms))
        for line in lines:
            if line.find(word2) != -1 :
                index_for_bonds = lines.index(line)
                break
        lines_to_iterate = list(filter(lambda line :  (lines.index(line) > index_for_bonds) and (lines.index(line) - index_for_bonds <= number_of_bonds), lines))
        for line in lines_to_iterate:
            line = line.split()
            bond_pair = (int(line[1]) - 1, int(line[2]) - 1)
            initial_array[bond_pair] = 1
            initial_array[bond_pair[::-1]] = 1
    return initial_array   
    
#From the absolute file path of a mol2 file, returns a list of all the unique atoms that appear in the file.
def atoms(mol2_file_path):
    with open(mol2_file_path, "r") as file:
        lines = file.readlines()

    index_for_atoms = 0
    for index, line in enumerate(lines):
        if line.startswith("@<TRIPOS>ATOM"):
            index_for_atoms = index
            break

    index_for_bonds = 0
    for index, line in enumerate(lines):
        if line.startswith("@<TRIPOS>BOND"):
            index_for_bonds = index
            break

    lines_to_iterate = [line for line in lines if index_for_atoms is not None and index_for_bonds is not None
                        and index_for_atoms < lines.index(line) < index_for_bonds]

    all_atoms = [line.split()[1] for line in lines_to_iterate]

    return all_atoms

#From the absolute path of a mol2 file, returns a standardized list of atom types.
def atom_types_standard(mol2_file_path):
    all_atoms = atoms(mol2_file_path)
    unique_atoms = list(dict.fromkeys(all_atoms))
    number_of_atoms = len(all_atoms)
    atom_types = [0]*number_of_atoms
    for atom in unique_atoms:
        index_of_atom = unique_atoms.index(atom)
        indices = [i for i, x in enumerate(all_atoms) if x == atom]
        for i in indices:
            atom_types[i] = index_of_atom + 1
    return np.array(atom_types)

#From a mol2 file and a dictionary with keys the atoms appearing in the file and whose associated values are the atom type numbering of the atoms appearing in the mol2 file,
#returns a numpy array of atom types.
def atom_types_from_dict(mol2_file_path,atom_types_dictionary):
    all_atoms = atoms(mol2_file_path)
    atom_types = [atom_types_dictionary[atom] for atom in all_atoms]
    return np.array(atom_types)

def charges2mol2(mol2_file_path,charges_path,output_folder=None):
    
    if type(output_folder) == type(None):
        output_folder = os.path.abspath(r"Output") 
    
    with open(mol2_file_path,'r') as f:
        lines_f = f.readlines()
        word1 = '@<TRIPOS>ATOM'
        word2 = '@<TRIPOS>BOND'
        index_for_atoms = 0
        index_for_bonds = 0
        index_for_charges = 0
        index_for_total_charge = 0
        for line in lines_f:
            if line.find(word1) != -1:
                index_for_atoms = lines_f.index(line)
                break
        for line in lines_f:
            if line.find(word2) != -1:
                index_for_bonds = lines_f.index(line)
                break
        filtered_lines_f = list(filter(lambda line : (lines_f.index(line) > index_for_atoms)
                                       and (lines_f.index(line) < index_for_bonds),lines_f))

        indices_for_atoms = [lines_f.index(line) for line in filtered_lines_f]

    keys = [line.split()[-1] for line in filtered_lines_f]

    keys = [(x,i) for i,x in enumerate(keys)]
    
    with open(charges_path,'r') as g:
        lines_g = g.readlines()[2:]
        values = [line.split()[1] for line in lines_g]
    
    mydict = {k: v for k, v in zip(keys, values)}
    
    mol2_file_path_copy_name = os.path.splitext(os.path.basename(mol2_file_path))[0] + '_charges.mol2'

    mol2_file_path_copy = os.path.join(output_folder,mol2_file_path_copy_name)
    
    with open(mol2_file_path, "r") as rf:
        with open(mol2_file_path_copy, "w") as wf:
            counter1 = 0
            counter2 = 0
            for line in rf.readlines():
                if counter1 not in indices_for_atoms:
                    wf.writelines(line)
                    counter1=counter1+1
                else:
                    original_charge = keys[counter2][0]
                    new_charge = mydict.get(keys[counter2])
                    if float(original_charge) < 0 and float(new_charge) < 0:
                        original_charge = original_charge[1:] #considero apenas a parte positiva da carga antiga
                        new_charge = new_charge[1:] #considero apenas a parte da carga antiga, assim no arquivo novo a nova carga 
                                                    #estará escrita de maneira correta como negativa, pois a parte negativa vem da carga antiga
                    if float(original_charge) < 0 and float(new_charge) >=0:
                        new_charge  = " " + new_charge
                    if float(original_charge) >= 0 and float(new_charge) < 0:
                        original_charge = " " + original_charge
                    line = line.replace(original_charge,new_charge)
                    wf.writelines(line)
                    counter1=counter1+1
                    counter2=counter2+1
    return mol2_file_path_copy


        
#From a mol2 file and a charges log file, returns a new mol2 file whose charges are given by the ones in the charge .log file.
def charges_log2mol2(mol2_file_path, charges_log_path,output_folder=None):

    if type(output_folder) == type(None):
        output_folder = os.path.abspath(r"Output") 

    with open(mol2_file_path,'r') as f:
        lines_f = f.readlines()
        word1 = '@<TRIPOS>ATOM'
        word2 = '@<TRIPOS>BOND'
        index_for_atoms = 0
        index_for_bonds = 0
        index_for_charges = 0
        index_for_total_charge = 0
        for line in lines_f:
            if line.find(word1) != -1:
                index_for_atoms = lines_f.index(line)
                break
        for line in lines_f:
            if line.find(word2) != -1:
                index_for_bonds = lines_f.index(line)
                break
        filtered_lines_f = list(filter(lambda line : (lines_f.index(line) > index_for_atoms)
                                       and (lines_f.index(line) < index_for_bonds),lines_f))

        indices_for_atoms = [lines_f.index(line) for line in filtered_lines_f]

        keys = [line.split()[-1] for line in filtered_lines_f]

        keys = [(x,i) for i,x in enumerate(keys)]

        with open(charges_log_path, 'r') as g:
            lines_g = g.readlines()

            word3 = 'Charges'
            word4 = 'Total charge:'

            for line in lines_g:
                    if line.find(word3) != -1:
                        index_for_charges = lines_g.index(line)
                        break
            for line in lines_g:
                    if line.find(word4) != -1:
                        index_for_total_charge = lines_g.index(line)
                        break

            filtered_lines_g = list(filter(lambda line : (lines_g.index(line) > index_for_charges + 1) 
                                 and (lines_g.index(line) < index_for_total_charge - 1),lines_g))

            values = [line.split()[-1] for line in filtered_lines_g]

        mydict = {k: v for k, v in zip(keys, values)}

        mol2_file_path_copy_name = os.path.splitext(os.path.basename(mol2_file_path))[0] + '_charges_log.mol2'

        mol2_file_path_copy = os.path.join(output_folder,mol2_file_path_copy_name)

        with open(mol2_file_path, "r") as rf:
            with open(mol2_file_path_copy, "w") as wf:
                counter1 = 0
                counter2 = 0
                for line in rf.readlines():
                    if counter1 not in indices_for_atoms:
                        wf.writelines(line)
                        counter1=counter1+1
                    else:
                        original_charge = keys[counter2][0]
                        new_charge = mydict.get(keys[counter2])
                        if float(original_charge) < 0 and float(new_charge) < 0:
                            original_charge = original_charge[1:] #considero apenas a parte positiva da carga antiga
                            new_charge = new_charge[1:] #considero apenas a parte da carga antiga, assim no arquivo novo a nova carga 
                                                        #estará escrita de maneira correta como negativa, pois a parte negativa vem da carga antiga
                        if float(original_charge) < 0 and float(new_charge) >=0:
                            new_charge  = " " + new_charge
                        if float(original_charge) >= 0 and float(new_charge) < 0:
                            original_charge = " " + original_charge
                        line = line.replace(original_charge,new_charge)
                        wf.writelines(line)
                        counter1=counter1+1
                        counter2=counter2+1
    return mol2_file_path_copy

#From a mol2 file and a dataframe with charge values, returns a mol2 file with charges given by the charges in the dataframe.
def dataframe2mol2(mol2_file_path,charges_dataframe,sheet=None,output_folder=None):
    if type(output_folder) == type(None):
        output_folder = os.path.abspath(r"Output")

    if type(sheet) == type(None):
        charges_ods = pd.read_excel(charges_dataframe,engine="odf")
        charges_ods = charges_ods[["Unnamed: 1", "Unnamed: 2"]][5:-1].reset_index(drop=True)
        charges_ods = charges_ods.rename(columns={"Unnamed: 1" : "Átomo", "Unnamed: 2": "Carga"})

        charges = charges_ods["Carga"].to_numpy()
    else:
        charges_ods = pd.read_excel(charges_dataframe,sheet_name = sheet, engine="odf")
        charges_ods = charges_ods[["Unnamed: 1", "Unnamed: 2"]][5:-1].reset_index(drop=True)
        charges_ods = charges_ods.rename(columns={"Unnamed: 1" : "Átomo", "Unnamed: 2": "Carga"})

        charges = charges_ods["Carga"].to_numpy()

        maximum_digits = max([len(str(abs(charge))) for charge in charges])

        new_charges = []
        for charge in charges:      
            new_charge = str(charge)
            if float(new_charge) < 0: 
                while len(new_charge) < maximum_digits+1:
                    new_charge = new_charge + '0' 
            else:
                while len(new_charge) < maximum_digits:
                    new_charge = new_charge + '0' 
            new_charges.append(new_charge)
    
    index_for_atoms = 0
    index_for_bonds = 0
    
    with open(mol2_file_path,'r') as f:
        lines_f = f.readlines()
        word1 = '@<TRIPOS>ATOM'
        word2 = '@<TRIPOS>BOND'
      
        for line in lines_f:
            if line.find(word1) != -1:
                index_for_atoms = lines_f.index(line)
            if line.find(word2) != -1:
                index_for_bonds = lines_f.index(line)
                break
        filtered_lines_f = list(filter(lambda line : (lines_f.index(line) > index_for_atoms)
                                       and (lines_f.index(line) < index_for_bonds),lines_f))

        indices_for_atoms = [lines_f.index(line) for line in filtered_lines_f]

        keys = [line.split()[-1] for line in filtered_lines_f]
        
        keys = [(x,i) for i,x in enumerate(keys)]

        mydict = {k: v for k, v in zip(keys, new_charges)}

        mol2_file_path_copy_name = os.path.splitext(os.path.basename(mol2_file_path))[0] + '_dataframe.mol2'

        mol2_file_path_copy = os.path.join(output_folder,mol2_file_path_copy_name)
        
        with open(mol2_file_path, "r") as rf:
            with open(mol2_file_path_copy, "w") as wf:
                counter1 = 0
                counter2 = 0
                for line in rf.readlines():
                    if counter1 not in indices_for_atoms:
                        wf.writelines(line)
                        counter1=counter1+1
                    else:
                        original_charge = keys[counter2][0]
                        new_charge = mydict.get(keys[counter2])
                        if float(original_charge) < 0 and float(new_charge) < 0:
                            original_charge = original_charge[1:] #considero apenas a parte positiva da carga antiga
                            new_charge = new_charge[1:] #considero apenas a parte da carga antiga, assim no arquivo novo a nova carga 
                                                        #estará escrita de maneira correta como negativa, pois a parte negativa vem da carga antiga
                        if float(original_charge) < 0 and float(new_charge) >=0:
                            new_charge  = " " + new_charge
                        if float(original_charge) >= 0 and float(new_charge) < 0:
                            original_charge = " " + original_charge
                        line = line.replace(original_charge,new_charge)
                        wf.writelines(line)
                        counter1=counter1+1
                        counter2=counter2+1
    return mol2_file_path_copy
#From a mol2 file, returns a numpy array with all the charges of the atoms appearing in the file.
def get_charges_from_mol2(mol2_file_path):
    with open(mol2_file_path, 'r') as f:
        lines = f.readlines()
        word1 = '@<TRIPOS>ATOM'
        word2 = '@<TRIPOS>BOND'
        charges = []
        index_for_atoms = 0
        index_for_bonds = 0
        for line in lines:
            if line.find(word1) != -1:
                index_for_atoms = lines.index(line)
            if line.find(word2) != -1:
                index_for_bonds = lines.index(line)
                break
        lines_to_iterate = list(filter(lambda line : (lines.index(line) > index_for_atoms)
                                   and (lines.index(line) < index_for_bonds),lines))
        for line in lines_to_iterate:
            charge = line.split()[-1]
            charges.append(charge)
        charges = np.array(charges)
    return charges 

#iterates over a folder and returns a tuple where the first entry corresponds to a list of all the absolute paths
#of the files inside the folder, while the second entry gives the total number of files inside it.
def get_list_of_file_paths(folder_path):
    files_directory = os.path.abspath(folder_path)
    files_in_folder = []
    for root, directories, files in os.walk(files_directory):
        for filename in files:
            file_path = os.path.join(root, filename)
            file_path = os.path.abspath(file_path)
            files_in_folder.append(file_path)

    files_in_folder = sorted(files_in_folder)
    number_of_files = len(files_in_folder)
    
    return files_in_folder, number_of_files


def func_unknown_bci(adjacency_matrix,atom_type):
    type_configurations = np.transpose(np.meshgrid(atom_type,atom_type),(2,1,0))
    type_configurations = type_configurations[adjacency_matrix.astype(bool)]
    num_rows, num_columns =np.shape(type_configurations)
    possible_indices = list(filter(lambda i: type_configurations[i][0] < type_configurations[i][1], range(0,num_rows)))
    unknown_bci = np.unique(type_configurations[possible_indices], axis=0)
    return unknown_bci
        
def all_bcis(mol2_folder_path,list_of_atom_types):
    #platinum_complexes_atom_types = [mol2_tools.atom_types_from_dict(file,atom_types) for file in platinum_complexes_file_path]
    all_bcis = set()
    
    mol2_files_path, number_of_structures = get_list_of_file_paths(mol2_folder_path)
    
    mol2_adjacency_matrices = [convert(file) for file in mol2_files_path]

    for i in range(number_of_structures):
        ith_bcis = func_unknown_bci(mol2_adjacency_matrices[i],list_of_atom_types[i])
        ith_bcis = set(tuple(x) for x in ith_bcis)
        all_bcis = all_bcis.union(ith_bcis)
    all_bcis = sorted(list(all_bcis))
    return all_bcis

#Gets the ith-bci
def ith_bci(mol2_folder_path,list_of_atom_types,i):

    mol2_files_path, number_of_structures = get_list_of_file_paths(mol2_folder_path)

    mol2_adjacency_matrices = [convert(file) for file in mol2_files_path]
    
    ith_bci = [tuple(row) for row in func_unknown_bci(mol2_adjacency_matrices[i],list_of_atom_types[i])]
    return ith_bci

#From a folder containing mol2 files and a list of atom types, returns a list whose elements are lists of atoms which share
#a given bci
def common_bcis(mol2_folder_path,list_of_atom_types):
    
    all_bcis_mol2 = all_bcis(mol2_folder_path,list_of_atom_types)
    
    mol2_files_path, number_of_structures = get_list_of_file_paths(mol2_folder_path)

    common_bcis = []

    for bci in all_bcis_mol2:
        structures_with_common_bci = list(filter(lambda i : bci in ith_bci(mol2_folder_path,list_of_atom_types,i), range(number_of_structures)))
        common_bcis.append(structures_with_common_bci)
    return common_bcis

#add comment later
def all_atoms(mol2_folder_path):
    mol2_files_path, number_of_structures = get_list_of_file_paths(mol2_folder_path)
    all_atoms = []
    for mol2_file in mol2_files_path:
        for atom in atoms(mol2_file):
            if atom not in all_atoms:
                all_atoms.append(atom)
    return all_atoms

#add comment later
def atom_types_dict(mol2_folder_path,list_of_atom_types):
    mol2_file_path,number_of_structures = get_list_of_file_paths(mol2_folder_path)
    
    appended_list_of_atom_types = []
    for num in range(number_of_structures):
        appended_list_of_atom_types.append(list_of_atom_types[num])
    appended_list_of_atom_types = list(chain(*appended_list_of_atom_types))

    appended_list_of_atoms = []
    for num in range(number_of_structures):
        appended_list_of_atoms.append(atoms(mol2_file_path[num]))
    appended_list_of_atoms = list(chain(*appended_list_of_atoms))

    dict_of_atom_types_correspondences = {(key,i): value for i, (key, value) in enumerate(zip(appended_list_of_atoms, appended_list_of_atom_types))}

    return dict_of_atom_types_correspondences

#reverses the key-value pairs of a dictionary    
def reverse_dict(mydict):
    new_dict = {}

    for k, v in mydict.items():
        new_dict[v] = k
    return new_dict

def list_of_keys_from_value(mydict,value_to_find):
    keys_list = []
    for key, value in mydict.items():
        if value == value_to_find:
            keys_list.append(key)
    return keys_list

#add comment later
def bci_atom_interactions(mol2_folder_path,list_of_atom_types):
    
    all_bcis_mol2 = all_bcis(mol2_folder_path,list_of_atom_types)
    
    atom_to_label = reverse_dict(atom_types_dict(mol2_folder_path,list_of_atom_types))

    atom_interactions = [(atom_to_label[bci[0]], atom_to_label[bci[1]]) for bci in all_bcis_mol2]

    bci_atom_interactions = {k: v for k, v in zip(all_bcis_mol2, atom_interactions)}

    return bci_atom_interactions

#add comment later
def atom_types_standard_folder(mol2_folder_path):
    mol2_files_path, number_of_structures = get_list_of_file_paths(mol2_folder_path)
    all_atoms_folder = all_atoms(mol2_folder_path)
    number_of_atoms = len(all_atoms_folder)
    numbering = list(range(1,number_of_atoms+1))
    mydict = {k: v for k, v in zip(all_atoms_folder, numbering)}
    list_of_atom_types = []
    
    for mol2_file in mol2_files_path:
        atom_type_for_file = atom_types_from_dict(mol2_file,mydict)
        list_of_atom_types.append(atom_type_for_file)
    return list_of_atom_types

def charges2mol2_folder(mol2_folder_path,charges_folder_path,output_folder=None):
    
    if type(output_folder) == type(None):
        output_folder = os.path.abspath(r"Output") 

    mol2_files_path, number_of_structures = get_list_of_file_paths(mol2_folder_path)
    charges_files_path = get_list_of_file_paths(charges_folder_path)[0]
    
    mydict = {k: v for k, v in zip(mol2_files_path,charges_files_path)}
    
    new_mol2_folder_name = os.path.basename(mol2_folder_path) + "_charges"
    new_mol2_folder_path = os.path.join(output_folder,new_mol2_folder_name)
    os.mkdir(new_mol2_folder_path)

    for mol2_file_path in mydict:
        charges_mol2_file_path = mydict.get(mol2_file_path)
        charges_mol2_file_path_name = os.path.basename(charges_mol2_file_path)
        print(charges_mol2_file_path_name)
        print(mol2_file_path)
        mol2_file_path_copy_source = charges2mol2(mol2_file_path,charges_mol2_file_path,output_folder=new_mol2_folder_path)

    return new_mol2_folder_path

#From a folder of mol2 files and charges which are represented by .log, both ordered so that the ith file of the 
#folder containing mol2 files has charges corresponding to the ith file in the folder containing .log files,
#returns a folder of new mol2 files each of which corresponds to the original mol2 files, but whose charges are 
#replaced by the ones in the charge .log files.

def charges_log2mol2_folder(mol2_folder_path,charges_folder_path,output_folder=None):

    if type(output_folder) == type(None):
        output_folder = os.path.abspath(r"Output") 

    mol2_files_path, number_of_structures = get_list_of_file_paths(mol2_folder_path)
    charges_files_path = get_list_of_file_paths(charges_folder_path)[0]
    
    mydict = {k: v for k, v in zip(mol2_files_path,charges_files_path)}
    
    new_mol2_folder_name = os.path.basename(mol2_folder_path) + "_charges_log"
    new_mol2_folder_path = os.path.join(output_folder,new_mol2_folder_name)
    os.mkdir(new_mol2_folder_path)

    for mol2_file_path in mydict:
        charges_log_file_path = mydict.get(mol2_file_path)
        charges_log_file_path_name = os.path.basename(charges_log_file_path)
        print(charges_log_file_path_name)
        print(mol2_file_path)
        mol2_file_path_copy_source = charges_log2mol2(mol2_file_path,charges_log_file_path,output_folder=new_mol2_folder_path)

    return new_mol2_folder_path
    
#From a folder with mol2 files and a dataframe with charge values, returns a folder with new mol2 files
#each of which whose are given by the charges in the corresponding sheet of the dataframe.

def dataframe2mol2_folder(mol2_folder_path,charges_dataframe,output_folder=None):

    if type(output_folder) == type(None):
        output_folder = os.path.abspath(r"Output") 

    mol2_files_path, number_of_structures = get_list_of_file_paths(mol2_folder_path)
    
    charges_sheets = pd.ExcelFile(charges_dataframe)
    
    charges_sheets_names = charges_sheets.sheet_names
    
    new_mol2_folder_name = os.path.basename(mol2_folder_path) + "_dataframe"
    new_mol2_folder_path = os.path.join(output_folder,new_mol2_folder_name)
    os.mkdir(new_mol2_folder_path)
    
    mydict = {k: v for k,v in zip(mol2_files_path,charges_sheets_names)}
    
    for mol2_file_path in mydict:
        charges_sheet_name = mydict.get(mol2_file_path)
        print(charges_sheet_name)
        print(mol2_file_path)
        mol2_file_path_copy_source = dataframe2mol2(mol2_file_path,charges_dataframe,charges_sheet_name,output_folder=new_mol2_folder_path)

    return new_mol2_folder_path