#from turtle import color
import lyapunov
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

print("record to the file")
print("begin")
with open('loc_max', 'wb') as fp:
    pickle.dump(lmax, fp)

with open('w_time', 'wb') as fp:
    pickle.dump(w_s, fp)

with open('lyapun1', 'wb') as fp:
    pickle.dump(lya1, fp)

with open('lyapun2', 'wb') as fp:
    pickle.dump(lya2, fp)

with open('lyapun3', 'wb') as fp:
    pickle.dump(lya3, fp)

with open('l_time', 'wb') as fp:
    pickle.dump(w_lya, fp)
print("end")

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

#------------------LYAPUNOV-REALISATION-----------------------------
# rX = list()
# rY = list()
# rZ = list()
# rT = list()
# w_q = 0.0065
# lyapunov.RK_ly(rX, rY, rZ, rT, w_q)
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
# plt.rcParams["figure.figsize"] = (14,9)
# plt.rcParams.update({'font.size': 40})
# plt.subplots_adjust(left=None, right= None, top=0.98, bottom = 0.15)

#fig, ax = plt.subplots()

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

# rX = list()
# rY = list()
# rZ = list()
# rT = list()
# w_q = 0.004
# q0 = 0.1
# threeSystem.retuRelQQ(rX, rY, rZ, rT, w_q, q0)
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

#------------------BIFURCATION----------------------------
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
