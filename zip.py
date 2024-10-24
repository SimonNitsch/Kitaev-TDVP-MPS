import numpy as np
import sys
import os
import shutil

arguments = sys.argv


def Convert(curdir, folders, destination, quantity):

    total = []    
    raw_quantity = quantity + "_raw.txt"
    new_quantity = quantity + ".txt"

    for f in folders:
        file = os.path.join(curdir,f,raw_quantity)
        data = np.loadtxt(file,delimiter=",")
        total.append(data)
    
    total = np.vstack(total)
    means = np.mean(total,axis=0)
    stds = np.sqrt(np.var(total,axis=0))

    new_file = os.path.join(curdir,destination,new_quantity)

    new_result = np.transpose(np.vstack([means,stds]))[1:,:]
    np.savetxt(new_file,new_result,delimiter=",")

    print("New Dataset saved as: %s" % new_file)
    




curdir = os.getcwd()
folders = arguments[1:-1]
destination = arguments[-1]

os.mkdir(destination)

Susceptx = os.path.isfile(os.path.join(curdir,arguments[1],"Chix_raw"))
Suscepty = os.path.isfile(os.path.join(curdir,arguments[1],"Chiy_raw"))
Susceptz = os.path.isfile(os.path.join(curdir,arguments[1],"Chiz_raw"))

xsrc = os.path.join(curdir,folders[0],"xdata.txt")
xdst = os.path.join(curdir,destination,"xdata.txt")
shutil.copy(xsrc,xdst)


Convert(curdir,folders,destination,"E")
Convert(curdir,folders,destination,"C")
Convert(curdir,folders,destination,"S")
Convert(curdir,folders,destination,"W")
Convert(curdir,folders,destination,"Mx")
Convert(curdir,folders,destination,"My")
Convert(curdir,folders,destination,"Mz")
Convert(curdir,folders,destination,"Mx2")
Convert(curdir,folders,destination,"My2")
Convert(curdir,folders,destination,"Mz2")

if Susceptx:
    Convert(curdir,folders,destination,"Chix")
if Suscepty:
    Convert(curdir,folders,destination,"Chiy")
if Susceptz:
    Convert(curdir,folders,destination,"Chiz")








