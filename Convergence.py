# Importing all necessary stuff
import math
import numpy as np
import matplotlib.pyplot as plt

from os import listdir
from os.path import isfile, join

mypath = "./output_files"
in_fname_list = [mypath+"/"+f for f in listdir(mypath) if isfile(join(mypath, f))]

data_dict = dict() #change later
data = []

x = 1
for in_fname in in_fname_list:
    new_dict = dict()
    Nr = in_fname.split("_")[x+1]
    Nt = in_fname.split("_")[x+3]
    alpha = in_fname.split("_")[x+5]
    Dimensions = in_fname.split("_")[x+7]

    new_dict["Nr"] = int(Nr)
    new_dict["Nt"] = int(Nt)
    new_dict["alpha"] = float(alpha)
    new_dict[Dimensions] = float(alpha)

    ff_arr = np.loadtxt(in_fname, usecols=(2, 3, 4, 5, 6, 7, 8, 9), skiprows=1)
    new_dict["field_re"] = ff_arr[0]
    new_dict["field_im"] = ff_arr[1]
    new_dict["impulse_re"] = ff_arr[2]
    new_dict["impulse_im"] = ff_arr[3]
    new_dict["next_field_re"] = ff_arr[4]
    new_dict["next_field_im"] = ff_arr[5]
    new_dict["next_impulse_re"] = ff_arr[6]
    new_dict["next_impulse_im"] = ff_arr[7]

    data.append(new_dict) #change later
    data_dict[int(Nr)] = new_dict

print("entries")
for entry in data_dict:
    print (entry)
print()
 
comp_field = "field_re"
dif_33_65       = np.std(data_dict[   33][comp_field] - data_dict[   65][comp_field])
dif_65_129      = np.std(data_dict[   65][comp_field] - data_dict[  129][comp_field])
dif_129_257     = np.std(data_dict[  129][comp_field] - data_dict[  257][comp_field])
dif_257_513     = np.std(data_dict[  257][comp_field] - data_dict[  513][comp_field])
dif_513_1025    = np.std(data_dict[  513][comp_field] - data_dict[ 1025][comp_field])
dif_1025_2049   = np.std(data_dict[ 1025][comp_field] - data_dict[ 2049][comp_field])
dif_2049_4097   = np.std(data_dict[ 2049][comp_field] - data_dict[ 4097][comp_field])
dif_4097_8193   = np.std(data_dict[ 4097][comp_field] - data_dict[ 8193][comp_field])
dif_8193_16385  = np.std(data_dict[ 8193][comp_field] - data_dict[16385][comp_field])
dif_16385_32769 = np.std(data_dict[16385][comp_field] - data_dict[32769][comp_field])

print("difs")
print (dif_33_65, dif_65_129, dif_129_257, dif_257_513)
print ()
ans = np.array([dif_33_65     /  dif_65_129,
                dif_65_129    /  dif_129_257,
                dif_129_257   /  dif_257_513,
                dif_257_513   /  dif_513_1025,
                dif_513_1025  /  dif_1025_2049,
                dif_1025_2049 /  dif_2049_4097,
                dif_2049_4097 /  dif_4097_8193,
                dif_4097_8193 /  dif_8193_16385,
                dif_8193_16385/  dif_16385_32769])

print ("relations")
print (ans)
print ()

print ("Convergence degree")
for a in ans:
    print (np.log(a)/np.log(2))
