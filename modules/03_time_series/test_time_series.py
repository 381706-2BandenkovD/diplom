#from turtle import color
import build.lib.time_series as time_series
import matplotlib.pyplot as plt
from math import *
from datetime import datetime
import time
from pathlib import Path
import pickle

# output_dir = '../../res/'
# input_dir = '../../res/data/'

# input_dir = Path(input_dir)
# output_dir = Path(output_dir)

#------------------LYAPUNOV-REALISATION-----------------------------
# rX = list()
# rY = list()
# rZ = list()
# rT = list()
# w_q = 0.0065
# time_series.RK_ly(rX, rY, rZ, rT, w_q)
# #plt.title("Qs = 0.01, θz = 5.5, θp =6, ω = 0.0")
# plt.plot(rT,rX, label = 'График fz(t)')
# plt.plot( rT,rY, label = 'График fp(t)')
# plt.plot( rT,rZ, label = 'График fr(t)')
# #plt.plot(rX,rZ)
# #plt.plot( rX,rY)
# plt.xlabel('Ось time ')
# plt.ylabel('Ось Z, P, R')
# #plt.xlabel('Z')
# #plt.ylabel('R')
# #plt.xlabel('ω$_s$')
# #plt.ylabel('Z$_m$$_a$$_x$') 
# plt.legend()
# plt.show()

#----------------------REALISATION-----------------------------------
plt.rcParams["figure.figsize"] = (14,9)
plt.rcParams.update({'font.size': 40})
plt.subplots_adjust(left=None, right= None, top=0.98, bottom = 0.15)

fig, ax = plt.subplots()
# for i in range (0, 100):
#     q0 = i / 100
#     rX = list()
#     rY = list()
#     rZ = list()
#     rT = list()
#     w_q = 0.004
#     time_series.retuRelease(rX, rY, rZ, rT, w_q)
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
time_series.retuRelease(rX, rY, rZ, rT, w_q)
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
