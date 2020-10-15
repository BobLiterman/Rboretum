#!/usr/bin/env python3
import subprocess
import os
import sys
import pickle

def rb_run_align_comp(package_dir,path_to_align,spp_info,align_name):
    cmd = ["python", package_dir+"/Alignment_Composition.py",path_to_align,spp_info,align_name]
    
    if os.name == 'nt':
        result = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,shell=True)
        stdout,stderr = result.communicate()
        pickle_path = stdout.decode("utf-8").rstrip()
    else:
        result=subprocess.check_output(cmd)
        pickle_path = result.decode("utf-8").rstrip()
    
    with open(pickle_path, 'rb') as pickle_file:
        content = pickle.load(pickle_file)
    os.remove(pickle_path)
    return(content)
    
def rb_run_species_comp(package_dir,path_to_align,spp_info,align_name):
    cmd = ["python", package_dir+"/Species_Composition.py",path_to_align,spp_info,align_name]
    
    if os.name == 'nt':
        result = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,shell=True)
        stdout,stderr = result.communicate()
        pickle_path = stdout.decode("utf-8").rstrip()
    else:
        result=subprocess.check_output(cmd)
        pickle_path = result.decode("utf-8").rstrip()
    
    with open(pickle_path, 'rb') as pickle_file:
        content = pickle.load(pickle_file)
    os.remove(pickle_path)
    return(content)

def rb_run_align_patterns(package_dir,path_to_align,use_gaps,spp_info,align_name):
    cmd = ["python", package_dir+"/Alignment_Patterns.py",path_to_align,spp_info,use_gaps,align_name]
    
    if os.name == 'nt':
        result = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,shell=True)
        stdout,stderr = result.communicate()
        pickle_path = stdout.decode("utf-8").rstrip()
    else:
        result=subprocess.check_output(cmd)
        pickle_path = result.decode("utf-8").rstrip()
    
    with open(pickle_path, 'rb') as pickle_file:
        content = pickle.load(pickle_file)
    os.remove(pickle_path)
    return(content)

def rb_run_align_splits(package_dir,path_to_align,use_gaps,spp_info,align_name):
    cmd = ["python", package_dir+"/Alignment_Splits.py",path_to_align,spp_info,use_gaps,align_name]
    
    if os.name == 'nt':
        result = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,shell=True)
        stdout,stderr = result.communicate()
        pickle_path = stdout.decode("utf-8").rstrip()
    else:
        result=subprocess.check_output(cmd)
        pickle_path = result.decode("utf-8").rstrip()
    
    with open(pickle_path, 'rb') as pickle_file:
        content = pickle.load(pickle_file)
    os.remove(pickle_path)
    return(content)