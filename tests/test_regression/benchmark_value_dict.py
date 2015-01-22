import string
import os
import sys

current_directory = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,current_directory)
print current_directory

def read_file():

    val_file = open('test_regression/benchmark_values.txt','r')

    benchmark_val_dict = {}
    
    for line in val_file:

        if line[0] == '#':
            continue
        line = line.strip()
        split_line = line.split(',')
        benchmark_name = split_line[0]
        test_type = split_line[1]
        value = float(split_line[2])

        benchmark_val_dict[(benchmark_name, test_type)] = value
    print benchmark_val_dict
    return benchmark_val_dict

benchmark_value_dictionary = read_file()
