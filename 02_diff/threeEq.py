#from turtle import color
import threeSystem
import matplotlib.pyplot as plt
from math import *
from datetime import datetime
import time
from pathlib import Path

output_dir = '../../res/'
input_dir = '../../res/data/'

# plt.rcParams["figure.figsize"] = (14,9)
# plt.rcParams.update({'font.size': 40})
# plt.subplots_adjust(left=None, right= None, top=0.98, bottom = 0.15)

input_dir = Path(input_dir)
output_dir = Path(output_dir)

fig, ax = plt.subplots()

# for i in range (0, 100):
#     q0 = i / 100
#     rX = list()
#     rY = list()
#     rZ = list()
#     rT = list()
#     w_q = 0.004
#     threeSystem.retuRelQQ(rX, rY, rZ, rT, w_q, q0)
#     ax.plot(rT,rX, color = 'red')
#     ax.plot(rT,rY, color = 'orange')
#     ax.plot(rT,rZ, color = 'green')
#     fig.savefig('res/realis_Qs_' + str(q0)+'.png')
#     ax.clear()



rX = list()
rY = list()
rZ = list()
rT = list()
w_q = 0.004
q0 = 0.1
threeSystem.retuRelQQ(rX, rY, rZ, rT, w_q, q0)
#plt.title("Qs = 0.01, θz = 5.5, θp =6, ω = 0.0")
plt.plot(rT,rX, label = 'График fz(t)')
plt.plot( rT,rY, label = 'График fp(t)')
plt.plot( rT,rZ, label = 'График fr(t)')
#plt.plot(rX,rZ)
#plt.plot( rX,rY)
plt.xlabel('Ось time ')
plt.ylabel('Ось Z, P, R')
#plt.xlabel('Z')
#plt.ylabel('R')
#plt.xlabel('ω$_s$')
#plt.ylabel('Z$_m$$_a$$_x$') 
plt.legend()

plt.show()

###################################################################################
# xera = list()
# retun = list()

# print("Bifurcation begin.. ")
# threeSystem.bifurcation(xera, retun)

# sizeL = len(xera)
# sizeT = len(retun)

#fig, ax = plt.subplots()

# for i in range (0, 21):
#     q0 = i / 10 + 3
#     xera = list()
#     retun = list()
#     threeSystem.bifurcation(xera, retun, q0)
#     ax.scatter(retun, xera,  s = 1.5, color= 'darkred')
#     fig.savefig('bifurk' + str(q0)+'.png')
#     ax.clear()

# plt.rcParams["figure.figsize"] = (14,9)
# plt.rcParams.update({'font.size': 40})
# plt.subplots_adjust(left=None, right= None, top=0.98, bottom = 0.15)

# xera = list()
# retun = list()
# threeSystem.bifurcation(xera, retun)
# plt.scatter(retun, xera,  s = 1.5, color= 'darkred')
# plt.xlabel('ω$_s$')
# plt.ylabel('Z$_m$$_a$$_x$') 
# plt.show()
