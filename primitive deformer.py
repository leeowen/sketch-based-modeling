"""
solve a Euler-Lagrange PDE equation
"""


import scipy as sp
import numpy as np
import re


def getpointcount(line):
    # get pointcount number from .geo file
    index=line.find("pointcount",0,len(line))
    count = -1
    if index != -1:
        count = int(re.search(r'\d+', line).group())
    return count


def main():
    print"this is main function:"
    # parse argument from .geo file
    try:
        f = open("position_constraint.geo","r")
    except IOError as reason:
        print("error :( and the reason is:" + str(reason))

    nextline = f.readline()
    ptnum = 0
    primnum =0

    while "" != nextline:
        if ptnum == 0:
            ptnum=getpointcount(nextline) if getpointcount(nextline) != -1 else 0
        if primnum == 0:
            primnum=getprimtivecount(nextline) if 
        nextline = f.readline()

    f.close()

    print "ptnum is:"+ str(ptnum)


if __name__ == "__main__":
    main()

'''
import sys
print "This is the name of the script: ", sys.argv[0]
print "Number of arguments: ", len(sys.argv)
print "The arguments are: " , str(sys.argv)
'''



