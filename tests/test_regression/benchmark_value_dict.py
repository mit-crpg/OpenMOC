import string
import os
import sys

current_directory = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,current_directory)
print current_directory

def read_file_to_list():

    val_file = open('benchmark_values.txt','r')
    file_line_list = []
    
    for line in val_file:

        if line[0] == '#':
            continue
        line = line.strip()
        split_line = line.split(',')
        file_line_list.append(split_line)

    val_file.close()
    return file_line_list

def create_dict(line_list):

    benchmark_value_dict = {}
    
    for line in line_list:
        test_type = line[0]
        val = line[1]
        benchmark_value_dict[test_type] = float(val)

    return benchmark_value_dict

def get_Keff(line_list):

    for line in line_list:
        if line[0] == 'Keff':
            Keff_value = float(line[1])
            return Keff_value

def update_Keff(Keff_value):

    val_file = open('benchmark_values.txt','rw')
    
    for line in val_file:
        if line[0] == 'Keff':
            line.replace(line[1], Keff_value)

    val_file.close()
    

    

line_list = read_file_to_list()
benchmark_value_dictionary = create_dict(line_list)
print 'getting Keff'
benchmark_val = get_Keff(line_list)
print 'Keff found was', benchmark_val
