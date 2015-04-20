#!/usr/bin/python -O
#Jason Matthew Torres
'''
Collection of python functions I created that do routine tasks 
'''
#libraries 
import sys, os
import subprocess as sp
import re



#globals
cout = sys.stdout.write

#functions

def list_from_file(input_file, col_num, rows_to_skip = None):
    '''
    create a list of desired column info from a data file
    can skip specified rows if there is one or more header lines
    '''
    file = open(input_file, 'r')
    if rows_to_skip:
        for i in range(rows_to_skip):
            file.readline()
    getlist = [] 
    for line in file:
        list = line.strip().split() 
        getlist.append(list[col_num-1]) 
    file.close()
    return getlist
    
def reduce_file(input_file, out_name, pyfield_to_eval, num_row_exempt = None, keep_list = None, remove_list = None, delim = "\t"):  
    '''
    function to evaluate content in a file with either a keep_list or remove_list of 
    strings. num_row_exempt is the number of row not to evaluate (i.e. headers). 
    out_name is the name of the desired reduced output file
    ''' 
    file = open(input_file, 'r') 
    out_file = open(out_name, 'w') 
    if num_row_exempt:
        for num in range(num_row_exempt):
            line = file.readline().strip() 
            out_file.write(line + "\n") 
    for line in file:
        list = line.strip().split() 
        if keep_list:
            if list[pyfield_to_eval] in keep_list: 
                out_file.write(delim.join(list) + "\n") 
        elif remove_list:
            if list[pyfield_to_eval] not in remove_list: 
                out_file.write(delim.join(list) + "\n")
        else:
            print "No Evaluation List Specified"
            break 
    file.close()
    out_file.close()
    return out_name, ("%s file has been reduced to %s file\n" % (input_file, out_name))
    
def eval_list_for_pattern(list, pattern, remove = False):
    '''
    Scan list for pattern matches and either return a list with all matches or without 
    matches if remove is specified 
    '''
    search_pattern = re.compile(pattern) 
    capture_list = []  
    for element in list:
        if re.search(search_pattern, element):
            capture_list.append(list.index(element))
    if remove:
        return_list = [i for j, i in enumerate(list) if j not in capture_list]
    else: 
        return_list = [i for j, i in enumerate(list) if j in capture_list] 
    return return_list 
         
def list_to_file(list, outname, header = None): 
    '''
    Write elements of list into an output file 
    if header is specified then that will be the first line in the output file
    ''' 
    fout = open(outname, 'w')
    if header:
        fout.write(header+"\n") 
    for i in list:
        fout.write(i+"\n")
    fout.close()
    return outname 
    
def replace_file_entry(filename, outname, pyfield, string_to_replace, replacer, delim="\t"): 
    '''
    Replaces a specific unwanted recurrent entry in a file field with a desired string 
    '''
    fin = open(filename, 'r')
    if filename == outname: 
        fout = open(outname+"_temp", 'w')
    else:
        fout = open(outname, 'w') 
    for line in fin:
        list = line.strip().split()
        if list[pyfield] == string_to_replace:
            list[pyfield] = replacer
            fout.write(delim.join(list)+"\n") 
        else: 
            fout.write(delim.join(list)+"\n") 
    fin.close()
    fout.close()
    if filename == outname:
        os.rename(outname+"_temp", filename) 
    return outname 


