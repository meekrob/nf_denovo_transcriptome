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



def main():

    for root,R1,R2 in read_R1_R2(sys.stdin):
        print(f"{root}: {os.path.basename(R1)}, {os.path.basename(R2)}")
        

if __name__ == "__main__": main()
