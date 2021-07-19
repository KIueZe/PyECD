import glob
import re
import math
from operator import itemgetter
import os

HARTREE_TO_KCAL = 627.509391
TEMPERATURE = 298.0 
GAS_CONSTANT = 0.001986

#Below are the index values for the master data structure.
NAME = 0; CONF_NUM = 1; ENERGY = 2; KCAL_E = 3; REL_E = 4; BOLTZMANN_FACTOR = 5; MOL_X = 6; IMAG_FREQUENCIES = 7

#Below are the index values for the proton and carbon chemical shift data substructures.
ATOM_NUMBER = 0 ; ISOTROPIC_VALUE = 1; REF_SHIFT = 2; WEIGHTED_SHIFT = 3;

def main():
    lofc = read_gaussian_outputfiles()
    filename = get_filename_prefix(lofc[-1]) + "-allconf_energy_and_imag_freq.csv"
    lofc_freq = read_gaussian_freq_outfiles(lofc)
    lofc_energy = read_gaussian_energy_outfiles(lofc)
    # name, conf-num, {state:[rlength, nm, uv], ...}    
    rlength_and_uv = get_rlength_and_uv_value(lofc_energy)
    
    # name, conf-num, free-energy in Hartree
    lofe = get_list_of_free_energies(lofc_freq)
    # name, conf-num, free-energy in Hartree, free-energy in kcal, relative energy, boltzmann weights, mol fraction
    lofe = boltzmann_analysis(lofe)
    # name, conf-num, free-energy in Hartree, free-energy in kcal, relative energy, boltzmann weights, mol fraction, imag_frequencies
    lofe = count_imaginary_freq(lofc_freq, lofe)
    # name, conf-num, {state:[rlength * mol fraction, nm, uv * mol fraction], ...}      
    boltzmann_rlength_and_uv = boltzmann_analysis_rlength_and_uv(rlength_and_uv, lofe)

    with open("test.cd.bil", 'w') as f:
        for entry in boltzmann_rlength_and_uv:
            for k in entry[2]:
                print(entry[2][k][1],"\t",entry[2][k][0], file=f)

    write_validate_csv(lofe, FILENAME=filename)
    os.popen("mkdir " +'for_sum_spectra')
    for f in glob.glob('*.heat'):
        os.popen( 'move ' + f + ' for_sum_spectra')
    for f in glob.glob('*.bil'):
        os.popen( 'move ' + f + ' for_sum_spectra')

    
def boltzmann_analysis(lofe):
    # name, conf-num, free-energy in Hartree, free-energy in kcal
    lofe = kcal_convert(lofe)
    minE = find_minimum_E(lofe)
    # name, conf-num, free-energy in Hartree, free-energy in kcal, relative energy
    lofe = calc_rel_E(lofe, minE)
    # name, conf-num, free-energy in Hartree, free-energy in kcal, relative energy, boltzmann weights
    lofe = calc_boltzmann_weights(lofe)
    denom = calc_boltzmann_denomenator(lofe)
    # name, conf-num, free-energy in Hartree, free-energy in kcal, relative energy, boltzmann weights, mol fraction
    lofe = calc_mol_fraction(lofe,denom)

    return lofe

def write_validate_csv(lofe, FILENAME):
    validate_writer = open(FILENAME,'w')
    lofe = sorted(lofe, key=itemgetter(REL_E))

    print("Conformers have been sorted in order of increasing relative energy with lowest energy structures first.", file=validate_writer)
    print("Filename, Energy (a.u.), Energy (kcal/mol), Relative Energy (kcal/mol), Boltzmann Factor, Equilibrium Mole Fraction, Number of Imaginary Frequencies", file=validate_writer)
    for conformation in lofe:
        print(conformation[NAME],",", conformation[ENERGY],",", conformation[KCAL_E],",", conformation[REL_E],",",conformation[BOLTZMANN_FACTOR],",", conformation[MOL_X],",", conformation[IMAG_FREQUENCIES], file=validate_writer)
    print(" ", file=validate_writer)
    
def kcal_convert(lofe):
    # Index: NAME = 0; CONF_NUM = 1; ENERGY = 2
    for entry in lofe:
        entry.append(entry[ENERGY] * HARTREE_TO_KCAL)

    return lofe

def find_minimum_E(lofe):
    minE = 0
    for entry in lofe:  #This finds the minimum energy
        if entry[KCAL_E] < minE:
            minE = entry[KCAL_E]
    
    return minE

def calc_rel_E(lofe, minE):
    for entry in lofe: 
        entry.append(entry[KCAL_E] - minE)

    return lofe

def calc_boltzmann_weights(lofe):
    for entry in lofe:
        entry.append(math.exp( (-1 * entry[REL_E]) / (TEMPERATURE * GAS_CONSTANT))) 

    return lofe

def calc_boltzmann_denomenator(lofe):
    Boltzmann_denomenator = 0
    for entry in lofe:
        Boltzmann_denomenator = Boltzmann_denomenator + entry[BOLTZMANN_FACTOR]

    return Boltzmann_denomenator

def calc_mol_fraction(lofe,Boltzmann_denomenator):
    for entry in lofe:
        entry.append(entry[BOLTZMANN_FACTOR]/Boltzmann_denomenator)

    return lofe

def count_imaginary_freq(lofc_freq, lofe):
    OPT_FREQ_OUTFILE_CONF_NUM = 1; LINE_POS_OF_FREQUENCY_A = 2; LINE_POS_OF_FREQUENCY_B = 3; LINE_POS_OF_FREQUENCY_C = 4;

    for file in lofc_freq:
        IMAG_FREQUENCIES = 0
        for line in file[2]:
            if "Frequencies -- " in line:
                freq_linesplit = line.split()
                if float(freq_linesplit[LINE_POS_OF_FREQUENCY_A]) < 0:
                    IMAG_FREQUENCIES += 1
                if float(freq_linesplit[LINE_POS_OF_FREQUENCY_B]) < 0:
                    IMAG_FREQUENCIES += 1
                if float(freq_linesplit[LINE_POS_OF_FREQUENCY_C]) < 0:
                    IMAG_FREQUENCIES += 1

        for entry in lofe:
            if entry[CONF_NUM] == file[OPT_FREQ_OUTFILE_CONF_NUM]:
                entry.append(IMAG_FREQUENCIES)       
    return lofe

def boltzmann_analysis_rlength_and_uv(rlength_and_uv, lofe):
    # CONF_NUM = 1; MOL_X = 6
    for file in rlength_and_uv:
        mol_fraction = 1
        for entry in lofe:
            if entry[1] == file[1]:
                mol_fraction = float(entry[6])
        for k in file[2]:
            file[2][k] = [float(file[2][k][0])*mol_fraction, float(file[2][k][1]), float(file[2][k][2])*mol_fraction]
    return rlength_and_uv

# 读取Gibbs自由能，生成并写入.heat文件
def get_list_of_free_energies(lofc_freq):
    LINE_POS_OF_FREE_ENERGY = 7  
    list_of_free_energies = []
    for file in lofc_freq:
        for line in file[2]:
            if "Free Energies=" in line:
                free_linesplit = line.split()
                free_energy = float(free_linesplit[LINE_POS_OF_FREE_ENERGY])
                list_of_free_energies.append([file[0],file[1],free_energy])
                with open("conf-%s.heat" % file[1], 'w') as heatwritor:
                    print("HEAT OF FORMATION  =  " + str(free_energy * HARTREE_TO_KCAL) + " kcal/mol", file=heatwritor)
    return list_of_free_energies

def get_conf_number(filename): 
    split_filename = re.findall(r'\w+', filename)
    rev_filename = split_filename[::-1]
    conf_number = rev_filename[1]
    return conf_number

# aquire the filename prefix
# TODO: check the consistency of all conf prefix
def get_filename_prefix(filename): 
    prefix_pos = filename.find("opt_freq")
    prefix = filename[:prefix_pos - 1]
    return prefix

# 匹配opt_freq和energy文件的数量，不匹配则退出进程,未调用
def count_file_types(lofc_freq, lofc_energy):
    if len(lofc_freq) > len(lofc_energy):
        print('''
   There are more opt-freq output files than energy output files.
   Please check that all energy calculations successfully finished 
   and that all .out files are properly contained within the
   proper working directory. To ensure that all conformer data 
   is considered, the program will now exit without creating the  
   final .csv files.
''')
        quit()

    elif len(lofc_freq) < len(lofc_energy):
        print('''
   There are more energy output files than opt-freq output files.
   Please check that all opt-freq .out files are properly 
   contained within the working directory. To ensure that 
   all conformer data is considered, the program will now 
   exit without creating the final .csv files.
''')
        quit()

# 读取各个激发态excited_state的编号，以及对应的CD振子强度R(length)和紫外强度f, 生成并写入.cd.bil文件
# read excited_state, and corresponding CD strength 'R(length)' and UV factor strength 'f'
def get_rlength_and_uv_value(lofc_energy):
    rlength_state_num = 0; rlength_value = 4
    state_num = 2; nm_value = 6; f_value = 8
    ret_rlength_and_uv_value = []
    for file in lofc_energy:
        excited_state_dic = {} # {'1': [R(length), nm_value, f_value], '2': [...], ...}
        i, j = 0, 0

        # locate, acquire value and write into dictionary (key <= state number, value <= R(length))
        while i < len(file[2]):
            if "R(length)" in file[2][i]:
                i += 1
                while "1/2[<0|del|b>*<b|r|0> + (<0|r|b>*<b|del|0>)*] (Au)" not in file[2][i]:
                    linesplit = file[2][i].split()
                    v = [linesplit[rlength_value]]
                    excited_state_dic[linesplit[rlength_state_num]] = v
                    i += 1 
            else:
                i += 1
        
        # locate, acquire value and write into dictionary (key <= state number, value <= nm_value, f_value)
        while j < len(file[2]):
            if " Excited State  " in file[2][j]:
                    linesplit = file[2][j].split()
                    l = excited_state_dic[linesplit[state_num][:-1]]
                    l.append(linesplit[nm_value])
                    l.append(linesplit[f_value][2:])
                    excited_state_dic[linesplit[state_num][:-1]] = l
            j += 1

        # output CD strength 'R(length)' .cd.bil file
        with open("conf-%s.cd.bil" % file[1], 'w') as f:
            for k in excited_state_dic:
                print(excited_state_dic[k][1],"\t",excited_state_dic[k][0], file=f)

        # output UV factor strength 'f' .uv.bil file, f_value
        with open("conf-%s.uv.bil" % file[1], 'w') as f:
            for k in excited_state_dic:
                print(excited_state_dic[k][1],"\t",excited_state_dic[k][2], file=f) 
        ret_rlength_and_uv_value.append([file[0],file[1],excited_state_dic])
    return ret_rlength_and_uv_value

def read_gaussian_freq_outfiles(list_of_files):
    list_of_freq_outfiles = []
    for file in list_of_files:
        if file.find('freq-') !=-1:
            list_of_freq_outfiles.append([file,int(get_conf_number(file)),open(file,"r").readlines()])
    return list_of_freq_outfiles

def read_gaussian_energy_outfiles(list_of_files):
    list_of_energy_outfiles = []
    for file in list_of_files:
        if file.find('energy-') !=-1:
            list_of_energy_outfiles.append([file,int(get_conf_number(file)),open(file,"r").readlines()])
    return list_of_energy_outfiles

def read_gaussian_outputfiles():
    list_of_files = []
    for file in glob.glob('*.out') or glob.glob('*.log'):
        list_of_files.append(file)
    return list_of_files

if __name__ == "__main__":
    main()
    print("""
        The script successfully performed the task of creating the %s file that shows the conformer number, conformer filename,
        total electronic energy, free energy, and total number of imaginary frequencies for each conformer.
        """ % "###")