# Importing all necessary stuff
import math
import numpy as np
import matplotlib.pyplot as plt

from os import listdir
from os.path import isfile, join

filename = "output.txt"

data = []
with open(filename) as file:
    for line in file:
        new_dict = dict()

        line = line.rstrip()
        task_id    = line.split(":")[1].split(",")[0]
        Action_re  = line.split(":")[2].split(",")[0]
        Action_im  = line.split(":")[3].split(",")[0]
        Energy_re  = line.split(":")[4].split(",")[0]
        Energy_im  = line.split(":")[5].split(",")[0]
        Energy_std = line.split(":")[6].split(",")[0]

        new_dict["Nr"]         =   int(task_id.split("_")[1])
        new_dict["Nt"]         =   int(task_id.split("_")[3])
        new_dict["alpha"]      = float(task_id.split("_")[5])
        new_dict["Dimensions"] = float(task_id.split("_")[7])

        new_dict["Action_re"]  = float(Action_re.split(",")[0])
        new_dict["Action_im"]  = float(Action_im.split(",")[0])
        new_dict["Energy_re"]  = float(Energy_re.split(",")[0])
        new_dict["Energy_im"]  = float(Energy_im.split(",")[0])

        new_dict["Energy_std"] = float(Energy_std)
        data.append(new_dict)

for dict in data:
    for entry in dict:
        print(entry, dict[entry])
    print()        
