"""
Ronald Randolph
CS594: Bioinformatics - HW2 Resubmission
This program computes the end-gap free alignment score of two sequences
in linear space using the following parameters: 

+2 for a match, -1 for a mismatch, and -2 for a gap. 
"""
import sys
import numpy as np

#function reads and stores sequences for alignment
def get_sequences():
    with open(sys.argv[1]) as f:
        data1 = f.read().splitlines()
        data1.pop(0)
	
    with open(sys.argv[2]) as f:
        data2 = f.read().splitlines()
        data2.pop(0)

    seq1 = "".join(data1)
    seq2 = "".join(data2)

    return seq1, seq2

#function creates 1-D list for sequence alignment
def initialize(n,m):

    #create substitution and traceback matrices
    T = np.zeros((m), dtype=int)
    return T

#function computes optimal score
def fill(T,s,t):
    f = [0,0,0]
    max_s = 0

    #Fill substitution matrix
    for i in range(1,len(s)+1):
        A = np.zeros((len(t)+1),dtype=int)
        for j in range(1,len(t)+1):
            f[0] = (T[j-1] + score(s[i-1],t[j-1]))
            f[1] = (A[j-1] - 2)
            f[2] = (T[j] - 2)
            
            A[j] = max(f)
            
        if(A[len(t)] >= max_s):
            max_s = A[len(t)]
        
        T = A
        #A = np.zeros((len(t)+1),dtype=int)

    for j in range(0,len(t)+1):
        if(T[j] >= max_s):
            max_s = T[j]

    #print score and return filled matrices
    print("\nOptimal Allignment Score: ",max_s)
    print()
    
    return

#function returns the score of matching two nucleotides
def score(a,b):
    if a == b:
        return 2
    else:
        return -1

def main():

    #argument check
    if(len(sys.argv) != 3):
        print("Usage: python3 egfalign_linear.py seq1.fasta seq2.fasta")
        sys.exit()

    s,t = get_sequences()

    n = len(s)+1
    m = len(t)+1

    T = initialize(n,m)
    fill(T,s,t)

if __name__ == "__main__":
    main()