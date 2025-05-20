#!/usr/bin/env python3
import sys,os

def get_root(filename):
    return filename.split('_R')[0] # like: 22_WNV_A_S1_L003_R1_001.fastq.gz => 22_WNV_A_S1_L003, 1_001.fastq.gz

def read_R1_R2(fh):
    lines = []
    last_root = None
    
    for line in fh:
        fields = line.strip().split()
        root = get_root(fields[0])
        if last_root is not None and root != last_root:
            R2 = lines.pop()
            R1 = lines.pop()
            yield last_root, R1, R2

        lines.append(fields[1])
        last_root = root

    
    # last pair
    R1 = lines.pop()
    R2 = lines.pop()
    yield root, R1,R2

def regularize_path(r_pth):
    basepath = "/nfs/home/rsbg/01_fastq"
    fullpath = os.path.join(basepath, r_pth)
    if not os.path.exists(fullpath):
        print(f"{basepath=}", file=sys.stderr)
        print(f"{r_pth=}", file=sys.stderr)
        print(f"{fullpath=}", file=sys.stderr)
        print(f"{os.path.exists(fullpath)=}", file=sys.stderr)
        raise FileNotFoundError

    return os.path.normpath(fullpath)


def main():

    for root,R1,R2 in read_R1_R2(sys.stdin):
        print(root, regularize_path(R1), regularize_path(R2), sep="\t")
        
        

if __name__ == "__main__": main()
