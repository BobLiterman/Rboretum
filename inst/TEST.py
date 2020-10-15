#!/usr/bin/env python3
import subprocess
import os
import sys
import pickle

def main(package_dir,path_to_align,spp_info,align_name):
    cmd = ["python", package_dir+"/Alignment_Composition.py",path_to_align,spp_info,align_name]
    result = subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,shell=True)
    print(result.communicate())

if __name__ == "__main__":
    main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])