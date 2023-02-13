#from turtle import color
import build.lib.bifurcation as bifurk
import matplotlib.pyplot as plt
from math import *
from datetime import datetime
import time
from pathlib import Path
import pickle

output_dir = '../../res/'
input_dir = '../../res/data/'

input_dir = Path(input_dir)
output_dir = Path(output_dir)

#------------------BIFURCATION----------------------------
lmax = list()
wq = list()

print("Bifurcation begin.. ")

lmax = list()
wq = list()
bifurk.bifurcation(lmax, wq)
plt.scatter(wq, lmax,  s = 1.5, color= 'darkred')
plt.xlabel('ω$_s$')
plt.ylabel('Z$_m$$_a$$_x$') 
plt.show()
