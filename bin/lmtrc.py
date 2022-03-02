#! /usr/bin/env python
#
#   lmtrc.py  rcfile key=val [key2=val2....] > rcfile2
#
#   lmtrc.py  rcfile [rcfile2 ...]  key=val [key2=val2....]
#    

import sys
import math
import numpy as np		
import os
import glob


if len(sys.argv) < 3:
    print("usage: %s rcfile [rcfile2...] key=val [key=val...]" % sys.argv[0])
    print("Existing keyword can change their value in the rcfile(s)")
    print("For a single rcfile, output is stdout")
    print("For multiple rcfile's, rcfile is changed")
    sys.exit(0)


nrc = 0    
for av in sys.argv[1:]:
    if os.path.exists(av):
        nrc = nrc + 1
    else:
        break

# print("Found %d RC files" % nrc)

rcfiles = sys.argv[1:nrc+1]


keyval = []
for argv in sys.argv[nrc+1:]:
    kv = argv.split('=')
    if len(kv) == 2:
        keyval.append([kv[0],kv[1]])


# print(keyval)

if len(keyval)==0:
    print("No keyword changes requested for the following rc files:")
    print(rcfiles)
    sys.exit(0)


for rcfile in rcfiles:
    lines = open(rcfile).readlines()
    m=0
    for i in range(len(lines)):
        line1 = lines[i].strip()
        n=0
        for kv in keyval:
            if line1.find("%s=" % kv[0])==0:
                n = n + 1
                m = m + 1
                line2 = "%s=%s" % (kv[0],kv[1])
                if nrc==1:
                    print(line2)
                else:
                    lines[i] = line2
        if n==0 and nrc==1:
            print(line1)
    if nrc > 1 and m>0:
        print("Changing %s with %d change(s)" % (rcfile,m))
        fp = open(rcfile,"w")
        for line in lines:
            fp.write("%s\n" % line.strip())
        fp.close()
        
        
        

