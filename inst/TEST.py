#!/usr/bin/env python3
import subprocess
import os
import sys
import pickle

def main(package_dir,path_to_align,spp_info,use_gaps,align_name):
    cmd = ["python", package_dir+"/Alignment_Splits.py",path_to_align,spp_info,use_gaps,align_name]
    
    if os.name == 'nt':
        result = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,shell=True)
        stdout,stderr = result.communicate()
        pickle_path = stdout.decode("utf-8").rstrip()
    else:
        result=subprocess.check_output(cmd)
        pickle_path = result.decode("utf-8").rstrip()
    print(pickle_path)

if __name__ == "__main__":
    main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5])