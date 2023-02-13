#from turtle import color
import build.lib.lyapunov as lyapunov
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

#------------------LYAPUNOV-LOCAL-MAX-------------------------
lmax = list()
w_s = list()
lya1 = list()
lya2 = list()
lya3 = list()
w_lya = list()

lyapunov.lyapunov_solution(lmax, w_s, lya1, lya2, lya3, w_lya)

# print("record to the file")
# print("begin")
# with open('loc_max', 'wb') as fp:
#     pickle.dump(lmax, fp)

# with open('w_time', 'wb') as fp:
#     pickle.dump(w_s, fp)

# with open('lyapun1', 'wb') as fp:
#     pickle.dump(lya1, fp)

# with open('lyapun2', 'wb') as fp:
#     pickle.dump(lya2, fp)

# with open('lyapun3', 'wb') as fp:
#     pickle.dump(lya3, fp)

# with open('l_time', 'wb') as fp:
#     pickle.dump(w_lya, fp)
# print("end")

fig = plt.figure(figsize=(14,9))
ax1 = fig.add_subplot(1,1,1)
ax1.scatter(w_s, lmax, s = 1.5, color= 'darkred')
## ax1.set_xlabel('r')
ax2 = ax1.twinx()
ax2.scatter(w_lya,lya1,s = 3, color= 'blue', label = 'l1')#s = 3, color= 'blue', 
ax2.scatter(w_lya,lya2,s = 3, color= 'orange', label = 'l2')#s = 3, color= 'orange',
ax2.scatter(w_lya,lya3,s = 3, color= 'green', label = 'l3')
ax2.legend()
##ax1.set_xlim(0.0045, 0.0079)
##ax1.set_ylim(-0.08, 5.5)
ax2.grid('on')
ax1.set_xlabel('ω$_s$')
ax1.set_ylabel('Z$_m$$_a$$_x$') 
plt.show()
