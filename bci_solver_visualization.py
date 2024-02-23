import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import scienceplots
from bci_solver_mol2tools import *
from bci_solver_optimization import *
import os

plt.style.use(['science', 'notebook', 'grid'])

def bci_solver_visualizer(mol2_folder_path,list_of_atom_types=None,q0=None,list_of_names=None,dpi_value=100,output_folder=None,show_atom_symbol=True,show_atom_type=False,save_fig=True,full_name=False):

    mol2_folder_path = os.path.abspath(mol2_folder_path)

    if type(output_folder) == type(None):
        if not os.path.exists(r"Output"):
            os.mkdir(r"Output")
        output_folder = os.path.abspath(r"Output")
    output_folder = os.path.abspath(output_folder)
    
    if type(list_of_atom_types) == type(None):
        list_of_atom_types = atom_types_standard_folder(mol2_folder_path)

    mol2_files_path, number_of_structures = get_list_of_file_paths(mol2_folder_path)

    mol2_files_names = [os.path.basename(mol2_file) for mol2_file in mol2_files_path]

    if list_of_names == None:
        if full_name == False:
            list_of_names = [mol2_file_name[:6] for mol2_file_name in mol2_files_names]
        else:
            list_of_names = [mol2_file_name for mol2_file_name in mol2_files_names]
            
    list_of_arrays_of_solutions = bci_solver_mol2_folder(mol2_folder_path,list_of_atom_types,q0,get_solutions_only=True)

    mol2_files_adjacency_matrices = [convert(mol2_file) for mol2_file in mol2_files_path]
    
    all_bcis_mol2 = all_bcis(mol2_folder_path,list_of_atom_types)
    
    common_bcis_mol2 = common_bcis(mol2_folder_path,list_of_atom_types)
    
    indexes_of_mol2_with_common_bci = 0
    
    if save_fig == True:

        duplicates = find_duplicates(output_folder)
        fig_save_folder_name = os.path.basename(mol2_folder_path) + "_bci_graphs"
        fig_save_folder = os.path.join(output_folder,fig_save_folder_name)

        if fig_save_folder_name in duplicates:
            count = duplicates[fig_save_folder_name]
            fig_save_folder_name = create_new_name(fig_save_folder_name,count)
            fig_save_folder = os.path.join(output_folder,fig_save_folder_name)        
        os.mkdir(fig_save_folder)

    for bci in all_bcis_mol2:
        index_of_bci = all_bcis_mol2.index(bci)
        mol2_names_with_common_bci = list(filter(lambda name : list_of_names.index(name) in common_bcis_mol2[index_of_bci], list_of_names))
        num_of_structures_with_bci = len(common_bcis_mol2[all_bcis_mol2.index(bci)])
        indexes_of_mol2_with_common_bci = common_bcis_mol2[index_of_bci]

        unknown_bcis = list(map(lambda i : func_unknown_bci(mol2_files_adjacency_matrices[i],list_of_atom_types[i]), indexes_of_mol2_with_common_bci))
        
        index_of_bci = []
        
        for i in indexes_of_mol2_with_common_bci:

            get_index_of_bci = [tuple(x) for x in unknown_bcis[indexes_of_mol2_with_common_bci.index(i)]]
            get_index_of_bci = get_index_of_bci.index(bci)
            index_of_bci.append(get_index_of_bci)

        bci_values = [list_of_arrays_of_solutions[i] for i in indexes_of_mol2_with_common_bci]
        
        bci_values = [bci_values[j][index_of_bci[j]] for j in range(len(indexes_of_mol2_with_common_bci))]
        
        mean_bci = np.mean(bci_values)

        bci_atom_interactions_mol2 = bci_atom_interactions(mol2_folder_path,list_of_atom_types)
        
        fig, axes = plt.subplots(figsize=(10, 5),dpi=dpi_value)
        fig.tight_layout(pad=5.0)
        x = list(range(1,num_of_structures_with_bci+1))
        plt.setp(axes, xticks=x, xticklabels=mol2_names_with_common_bci)
        plt.setp(axes.get_xticklabels(), fontsize=12)
        axes.set_xticklabels(axes.get_xticklabels(), rotation=45, ha='right')
        axes.axhline(mean_bci,color='red', linewidth=1, linestyle='dashed', label="Average bci value")
        plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        
        if (show_atom_symbol == True) and (show_atom_type == False):

            axes.set_title(label= f"bci values for the interactions {bci_atom_interactions_mol2[bci][0][0]}-{bci_atom_interactions_mol2[bci][1][0]}")  

        if (show_atom_symbol == True) and (show_atom_type == True):

            axes.set_title(label= f"bci values for the interactions {bci_atom_interactions_mol2[bci][0][0]}({bci[0]})-{bci_atom_interactions_mol2[bci][1][0]}({bci[1]})")
        
        if (show_atom_symbol == False) and (show_atom_type == True):

            axes.set_title(label= f"bci values for the interactions {bci[0]}-{bci[1]}")

        axes.plot(x,bci_values,'o--')  

        if save_fig == True:
            if (show_atom_symbol == True) and (show_atom_type == False):
                fig_name = f'bci_values_{bci_atom_interactions_mol2[bci][0][0]}-{bci_atom_interactions_mol2[bci][1][0]}.png'

            if (show_atom_symbol == True) and (show_atom_type == True):
                fig_name = f'bci_values_{bci_atom_interactions_mol2[bci][0][0]}({bci[0]})-{bci_atom_interactions_mol2[bci][1][0]}({bci[1]}).png'

            if (show_atom_symbol == False) and (show_atom_type == True):
                fig_name = f'bci_values_{bci[0]}-{bci[1]}.png'

            fig_name_path = os.path.join(fig_save_folder,fig_name)

            plt.savefig(fig_name_path)
            

mmff94_bcis = os.path.abspath(r"Openbabel-Files/mmffchg.par")

if os.path.exists(mmff94_bcis):

    with open(mmff94_bcis, 'r') as f:
        lines = f.readlines()
        word1 = '*  types       bci     Source'.split()
        bcis = []
        last_line = len(lines) - 1 
        for line in lines:
            if line.split() == word1:
                index_for_bcis = lines.index(line)
                break
        lines_to_iterate = list(filter(lambda line : (lines.index(line) > index_for_bcis) and (lines.index(line) < last_line), lines))
        for line in lines_to_iterate:
            bci = line.split()[3]
            bcis.append(bci)
        bcis = [float(bci) for bci in bcis]

    
def bci_solver_hist_openbabel(dpi_value=100,output_folder=None,show_atom_symbol=True,show_atom_type=False,
                          save_fig=True):
        
    fig, ax = plt.subplots(dpi=dpi_value)
    ax.set_xlim(-1,0.8)
    n_bins = 80
    ax.hist(bcis,bins=n_bins,density=False,histtype='step',linewidth=1.5,label='Histogram depicting the bci values available in OpenBabel',color='purple')
    plt.legend(loc='center left',bbox_to_anchor=(1,0.5))
    ax.set_xlabel('bci values')
    ax.set_ylabel('Frequency')

def bci_solver_hist_visualizer(mol2_folder_path,list_of_atom_types=None,q0=None,list_of_names=None,dpi_value=100,output_folder=None,show_atom_symbol=True,show_atom_type=False,save_fig=True,full_name=False):

    mol2_folder_path = os.path.abspath(mol2_folder_path)
    
    if type(output_folder) == type(None):
        if not os.path.exists(r"Output"):
            os.mkdir(r"Output")
        output_folder = os.path.abspath(r"Output")
    output_folder = os.path.abspath(output_folder)
    
    if type(list_of_atom_types) == type(None):
        list_of_atom_types = atom_types_standard_folder(mol2_folder_path)
    
    mol2_files_path, number_of_structures = get_list_of_file_paths(mol2_folder_path)

    mol2_files_names = [os.path.basename(mol2_file) for mol2_file in mol2_files_path]

    if list_of_names == None:
        if full_name == False:
            list_of_names = [mol2_file_name[:6] for mol2_file_name in mol2_files_names]
        else:
            list_of_names = [mol2_file_name for mol2_file_name in mol2_files_names]
            
    list_of_arrays_of_solutions = bci_solver_mol2_folder(mol2_folder_path,list_of_atom_types,q0,get_solutions_only=True)

    mol2_files_adjacency_matrices = [convert(mol2_file) for mol2_file in mol2_files_path]
    
    all_bcis_mol2 = all_bcis(mol2_folder_path,list_of_atom_types)
    
    common_bcis_mol2 = common_bcis(mol2_folder_path,list_of_atom_types)
    
    indexes_of_mol2_with_common_bci = 0
    
    if save_fig == True:

        duplicates = find_duplicates(output_folder)
        
        fig_save_folder_name = os.path.basename(mol2_folder_path) + "_bci_histograms"
        fig_save_folder = os.path.join(output_folder,fig_save_folder_name)

        if fig_save_folder_name in duplicates:
            count = duplicates[fig_save_folder_name]
            fig_save_folder_name = create_new_name(fig_save_folder_name,count)
            fig_save_folder = os.path.join(output_folder,fig_save_folder_name)        
        
        os.mkdir(fig_save_folder)

    for bci in all_bcis_mol2:
        index_of_bci = all_bcis_mol2.index(bci)
        mol2_names_with_common_bci = list(filter(lambda name : list_of_names.index(name) in common_bcis_mol2[index_of_bci], list_of_names))
        num_of_structures_with_bci = len(common_bcis_mol2[all_bcis_mol2.index(bci)])
        indexes_of_mol2_with_common_bci = common_bcis_mol2[index_of_bci]

        unknown_bcis = list(map(lambda i : func_unknown_bci(mol2_files_adjacency_matrices[i],list_of_atom_types[i]), indexes_of_mol2_with_common_bci))
        
        
        index_of_bci = []
        
        for i in indexes_of_mol2_with_common_bci:

            get_index_of_bci = [tuple(x) for x in unknown_bcis[indexes_of_mol2_with_common_bci.index(i)]]
            get_index_of_bci = get_index_of_bci.index(bci)
            index_of_bci.append(get_index_of_bci)
              
        bci_values = [list_of_arrays_of_solutions[i] for i in indexes_of_mol2_with_common_bci]
        
        bci_values = [bci_values[j][index_of_bci[j]] for j in range(len(indexes_of_mol2_with_common_bci))]
        
        mean_bci = np.mean(bci_values)

        bci_atom_interactions_mol2 = bci_atom_interactions(mol2_folder_path,list_of_atom_types)
        
        n_bins1 = 80
        n_bins2 = len(indexes_of_mol2_with_common_bci) 

        fig, ax = plt.subplots()
        ax.set_xlim(-1,0.8)
        ax.set_xlabel('Valores de bci')
        ax.set_ylabel('Frequência')

        ax.hist(bcis,bins=n_bins1,density=True,histtype='step',linewidth=1.5,label='Histogram depicting the bci values available in OpenBabel',color='purple')  
        
        if (show_atom_symbol == True) and (show_atom_type == False):

            ax.hist(bci_values,bins=n_bins2,density=True,histtype='step',linewidth=1.5,label=f'Histogram of the bci values for the interactions {bci_atom_interactions_mol2[bci][0][0]}-{bci_atom_interactions_mol2[bci][1][0]} ',color='green')
            
        if (show_atom_symbol == True) and (show_atom_type == True):
            
            ax.hist(bci_values,bins=n_bins2,density=True,histtype='step',linewidth=1.5,label=f'Histogram of the bci values for the interactions {bci_atom_interactions_mol2[bci][0][0]}({bci[0]})-{bci_atom_interactions_mol2[bci][1][0]}({bci[1]}) ',color='green')
            
        
        if (show_atom_symbol == False) and (show_atom_type == True):
            
            ax.hist(bci_values,bins=n_bins2,density=True,histtype='step',linewidth=1.5,label=f'Histogram of the bci values for the interactions {bci[0]}-{bci[1]} ',color='green')

        if save_fig == True:
            if (show_atom_symbol == True) and (show_atom_type == False):
                fig_name = f'bci_values_{bci_atom_interactions_mol2[bci][0][0]}-{bci_atom_interactions_mol2[bci][1][0]}_hist.png'

            if (show_atom_symbol == True) and (show_atom_type == True):
                fig_name = f'bci_values_{bci_atom_interactions_mol2[bci][0][0]}({bci[0]})-{bci_atom_interactions_mol2[bci][1][0]}({bci[1]})_hist.png'

            if (show_atom_symbol == False) and (show_atom_type == True):
                fig_name = f'bci_values_{bci[0]}-{bci[1]}_hist.png'

            fig_name_path = os.path.join(fig_save_folder,fig_name)

            plt.savefig(fig_name_path)
        
        plt.legend(loc='center left',bbox_to_anchor=(1,0.5)) 
