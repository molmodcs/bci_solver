import numpy as np
import scipy.linalg
from bci_solver_mol2tools import *
import os
    
def bci_matrix(adjacency_matrix,atom_type):
    unknown_bci = func_unknown_bci(adjacency_matrix,atom_type)
    number_of_unknowns = len(unknown_bci)
    type_configurations = np.transpose(np.meshgrid(atom_type,atom_type),(2,1,0))
    number_of_atoms = len(atom_type)
    bci_matrix = np.array([])
    row_number = 0
    for row in adjacency_matrix:
        type_configuration_in_row = type_configurations[row_number]
        for bci in unknown_bci:
            columns_where_type_configuration_equals_bci = np.fromiter(filter(lambda i : (type_configurations[row_number][i] == bci).all() or (type_configurations[row_number][i] == bci[::-1]).all(), range(0,number_of_atoms)), dtype=int)
            if atom_type[row_number] == np.min(bci):
                column_bci_matrix_row = -np.sum(row[columns_where_type_configuration_equals_bci])
                bci_matrix = np.append(bci_matrix,column_bci_matrix_row)
            else:
                column_bci_matrix_row = np.sum(row[columns_where_type_configuration_equals_bci])
                bci_matrix = np.append(bci_matrix,column_bci_matrix_row)
        row_number+=1    
    return np.reshape(bci_matrix,(number_of_atoms,number_of_unknowns))
       
def bci_solver(adjacency_matrix,atom_type,dq):
    least_squares_solution, residues, rank, singular_values = scipy.linalg.lstsq(bci_matrix(adjacency_matrix,atom_type),dq)
    return least_squares_solution

def bci_solver_mol2(mol2_file_path,atom_type=None,q0=0,show_atom_symbol=True,show_atom_type=False,get_solutions_only = False,output_folder=None,
                   decimal_places=4):
    mol2_file_path = os.path.abspath(mol2_file_path)
    
    if type(output_folder) == type(None):
        if not os.path.exists(r"Output"):
            os.mkdir(r"Output")
        output_folder = os.path.abspath(r"Output") 
    output_folder = os.path.abspath(output_folder)
    
    if type(atom_type) == type(None):
        atom_type = atom_types_standard(mol2_file_path)
    
    adjacency_matrix = convert(mol2_file_path)
    charges = get_charges_from_mol2(mol2_file_path)
    charges_float = np.array([float(charge) for charge in charges])
    dq = charges_float - q0
    solutions = bci_solver(adjacency_matrix,atom_type,dq)
    
    if get_solutions_only == True:
        return solutions
    
    else:

        duplicates = find_duplicates(output_folder)
        
        bcis = func_unknown_bci(adjacency_matrix,atom_type)
        bcis = [tuple(bci) for bci in bcis]

        mydict = {k: v for k, v in zip(bcis, solutions)}

        atoms_in_mol2 = atoms(mol2_file_path)

        atom_type_correspondence = {(key,i): value for i, (key, value) in enumerate(zip(atoms_in_mol2, atom_type))}
    
        mol2_file_path_copy_name = os.path.splitext(os.path.basename(mol2_file_path))[0] + "_bci"

        if mol2_file_path_copy_name in duplicates:
            count = duplicates[mol2_file_path_copy_name]
            mol2_file_path_copy_name = create_new_name(mol2_file_path_copy_name,count) + ".par"
            mol2_file_path_copy = os.path.join(output_folder,mol2_file_path_copy_name)
        else:
            mol2_file_path_copy = os.path.join(output_folder,mol2_file_path_copy_name + ".par")

        with open(mol2_file_path_copy, "w") as wf:
                first_entry_atom_symbols = []
                second_entry_atom_symbols = []
                atom_and_type_first = []
                atom_and_type_second = []
                first_entry_atom_type = []
                second_entry_atom_type = []
                bci_values = []
                for bci in bcis:
                    first_entry_type = bci[0]
                    first_entry_atom_type.append(first_entry_type)
                    atoms_with_corresponding_type_first = list_of_keys_from_value(atom_type_correspondence,
                                                                                                 first_entry_type)

                    first_entry_atom_symbol = atoms_with_corresponding_type_first[0][0]
                    first_entry_atom_symbols.append(first_entry_atom_symbol)
                    atom_and_type_first_for_bci = f'{first_entry_atom_symbol}({bci[0]})' 
                    
                    atom_and_type_first.append(atom_and_type_first_for_bci)
                    
                    second_entry_type = bci[1]
                    second_entry_atom_type.append(second_entry_type)
                    atoms_with_corresponding_type_second = list_of_keys_from_value(atom_type_correspondence,
                                                                                                 second_entry_type)
                    second_entry_atom_symbol = atoms_with_corresponding_type_second[0][0]
                    second_entry_atom_symbols.append(second_entry_atom_symbol)
                    atom_and_type_second_for_bci = f'{second_entry_atom_symbol}({bci[1]})'
                    
                    atom_and_type_second.append(atom_and_type_second_for_bci)
                    
                    bci_value = round(mydict.get(bci),decimal_places)

                    bci_values.append(bci_value)
                    
                maximum_digits = max([len(str(abs(bci_value))) for bci_value in bci_values])
                
                new_bci_values = []
                for bci_value in bci_values:
                    
                    if bci_value >=0:
                                
                        new_bci_value = ' '
                        
                    else:
                        
                        new_bci_value = ''
                        
                    new_bci_value = new_bci_value + str(bci_value)
                        
                    while len(new_bci_value) < maximum_digits+1:
                        new_bci_value = new_bci_value + '0'
                        
                    new_bci_values.append(new_bci_value)
                    
            
                if (show_atom_symbol == True) and (show_atom_type == False):
                    wf.write('*  atoms       bci \n')
                    
                    col_width = 3
                    col_spacing = 2

                    rows = []
                    for i in range(len(bcis)):
                        ith_row = ['0',first_entry_atom_symbols[i],second_entry_atom_symbols[i],new_bci_values[i]]
                        rows.append(ith_row)

                    for row in rows:
                        for col in row:
                            if row.index(col) != 0:
                                wf.write(col.ljust(col_width))
                                wf.write(' ' * col_spacing)
                            else:
                                wf.write(col.ljust(1))
                                wf.write(' ' * 3)
                        wf.write('\n')
                
                elif (show_atom_symbol == True) and (show_atom_type == True):
                    wf.write('*       atoms         bci \n')
                    
                    col_width = 5
                    col_spacing = 3

                    rows = []
                    for i in range(len(bcis)):
                        ith_row = ['0',atom_and_type_first[i],atom_and_type_second[i],new_bci_values[i]]
                        rows.append(ith_row)

                    for row in rows:
                        for col in row:
                            if row.index(col) != 0:
                                wf.write(col.ljust(col_width))
                                wf.write(' ' * col_spacing)
                            else:
                                wf.write(col.ljust(1))
                                wf.write(' ' * 3)
                        wf.write('\n')
        
                elif (show_atom_symbol == False) and (show_atom_type == True):
                    wf.write('* types      bci \n')
                    
                    col_width = 2
                    col_spacing = 2
                    
                    rows = []
                    for i in range(len(bcis)):
                        ith_row = ['0',first_entry_atom_type[i],second_entry_atom_type[i],new_bci_values[i]]
                        rows.append(ith_row)

                    for row in rows:
                        for col in row:
                            if row.index(col) != 0:
                                wf.write(str(col).ljust(col_width))
                                wf.write(' ' * col_spacing)
                            else:
                                wf.write(str(col).ljust(1))
                                wf.write(' ' * 2)
                        wf.write('\n')

                else:
                    print("The output file must show the atom symbols or the atom types")
                    
        print(mol2_file_path_copy_name + ' Done!')
        return solutions, mol2_file_path_copy_name
                                        
def bci_solver_mol2_folder(mol2_folder_path,list_of_atom_types=None,q0=None,show_atom_symbol=True,show_atom_type=False,
                             get_solutions_only=False,output_folder=None):

    mol2_folder_path = os.path.abspath(mol2_folder_path)
    
    if type(output_folder) == type(None):
        if not os.path.exists(r"Output"):
            os.mkdir(r"Output")
        output_folder = os.path.abspath(r"Output") 
    output_folder = os.path.abspath(output_folder)
    
    mol2_files_path, number_of_structures = get_list_of_file_paths(mol2_folder_path)
    if list_of_atom_types == None:
        list_of_atom_types = atom_types_standard_folder(mol2_folder_path)
    if q0 == None:
        q0 = np.array([0]*number_of_structures)   

    if get_solutions_only == True:
        solutions_to_all_problems = []
        for mol2_file in mol2_files_path:
            index_of_mol2_file = mol2_files_path.index(mol2_file)
            bcis_of_mol2_file = bci_solver_mol2(mol2_file,list_of_atom_types[index_of_mol2_file],q0[index_of_mol2_file],
                                                show_atom_symbol,show_atom_type,get_solutions_only)
            solutions_to_all_problems.append(bcis_of_mol2_file)

        list_of_arrays_of_solutions = [solution for solution in solutions_to_all_problems]

        return list_of_arrays_of_solutions

    else:

        duplicates = find_duplicates(output_folder)
        new_mol2_folder_name = os.path.basename(mol2_folder_path) + "_bci_solutions"
        new_mol2_folder_path = os.path.join(output_folder,new_mol2_folder_name)

        if new_mol2_folder_name in duplicates:
            count = duplicates[new_mol2_folder_name]
            new_mol2_folder_name = create_new_name(new_mol2_folder_name,count)
            new_mol2_folder_path = os.path.join(output_folder,new_mol2_folder_name)
        os.mkdir(new_mol2_folder_path)

        solutions_to_all_problems = []
        for mol2_file in mol2_files_path:
            
            index_of_mol2_file = mol2_files_path.index(mol2_file)
            bcis_of_mol2_file = bci_solver_mol2(mol2_file,list_of_atom_types[index_of_mol2_file],q0[index_of_mol2_file],
                                                show_atom_symbol,show_atom_type,get_solutions_only,output_folder=new_mol2_folder_path)
            solutions_to_all_problems.append(bcis_of_mol2_file)

        list_of_arrays_of_solutions = [solution[0] for solution in solutions_to_all_problems]
       

        return list_of_arrays_of_solutions, new_mol2_folder_name
