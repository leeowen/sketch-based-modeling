"""
solve a Euler-Lagrange PDE equation
"""


import scipy as sp
import numpy as np
def getpointcount(line):
    # type: (object) -> object
    index=line.find(line,0,len(line))
    count=-1
    if index!=-1:
        count = ''.join(x for x in line if x.isdigit())
        print int(count)
    return count

def main():
    #parse argument from .geo file
    f = open("position.geo","r")
    next=f.readline()
    pointcount = 0
    while next !="":
        pointcount=getpointcount(next) if getpointcount(next)!=-1 else 0
        next = f.readline()
    print pointcount
    f.close()


if __name__ == "__main__":
    main()

import sys
print "This is the name of the script: ", sys.argv[0]
print "Number of arguments: ", len(sys.argv)
print "The arguments are: " , str(sys.argv)


