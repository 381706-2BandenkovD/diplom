import matplotlib.pyplot as plt
from numpy.linalg import norm
import numpy as np
import pickle

def f(r, t, w_q):
    Z = r[0]
    P = r[1]
    R = r[2]

    fx = (-(a_z + g_p * P) * Z + b_z * (Z0 - (Z0 - Z1) / (1 + np.exp(-(Q0 + alfa_0 * Z + Qs * np.sin(w_q * t) - q_z) / k_z))))
    fy = (-a_p * P + b_p * (P0 - (P0 - P1) / (1 + np.exp(-(Q0 + alfa_0 * Z + Qs * np.sin(w_q * t) - q_p) / k_p))))
    fz = (-a_r * R + b_r * (R0 - (R0 - R1) / (1 + np.exp(-(Q0 + alfa_0 * Z + Qs * np.sin(w_q * t) - q_r) / k_r))))

    return np.array([fx,fy,fz], float)

def jacobian(r, t, w_q):
    M = np.zeros([3,3])
    Z = r[0]
    P = r[1]
    R = r[2]
    
    fx_coff = (np.exp((-Q0 - Qs*np.sin(t*w_q) - Z*alfa_0 + q_z)/k_z) + 1)
    fy_coff = (np.exp((-Q0 - Qs*np.sin(t*w_q) - Z*alfa_0 + q_p)/k_p) + 1)
    fz_coff = (np.exp((-Q0 - Qs*np.sin(t*w_q) - Z*alfa_0 + q_r)/k_r) + 1)

    fx_x = (-P*g_p - a_z - alfa_0*b_z*(Z0 - Z1)*np.exp((-Q0 - Qs*np.sin(t*w_q) - Z*alfa_0 + q_z)/k_z)/(k_z * fx_coff * fx_coff))
    fx_y = -Z*g_p
    fx_z = 0
    fy_x =  (-alfa_0*b_p*(P0 - P1)*np.exp((-Q0 - Qs*np.sin(t*w_q) - Z*alfa_0 + q_p)/k_p)/(k_p * fy_coff * fy_coff))
    fy_y =  -a_p
    fy_z =  0
    fz_x =  (-alfa_0*b_r*(R0 - R1)*np.exp((-Q0 - Qs*np.sin(t*w_q) - Z*alfa_0 + q_r)/k_r)/(k_r * fz_coff * fz_coff))
    fz_y =  0
    fz_z =  -a_r
    
    M[0,:] = [fx_x, fx_y, fx_z]
    M[1,:] = [fy_x, fy_y, fy_z]
    M[2,:] = [fz_x, fz_y, fz_z]

    return M

def g(d, r, t, w_q):
    dx = d[0]
    dy = d[1]
    dz = d[2]

    M = jacobian(r, t, w_q)

    dfx = np.dot(M, dx)
    dfy = np.dot(M, dy)
    dfz = np.dot(M, dz)

    return np.array([dfx, dfy, dfz], float)

def RkLyapunovRealis(x, y, z, time, w_q):
    d = np.array([[1,0,0], [0,1,0], [0,0,1]], float)
    r = np.array([x0, y0, z0], float)

    l1, l2, l3 = 0, 0, 0

    x.append(r[0])
    y.append(r[1])
    z.append(r[2])
    time.append(0)

    t = 0
    i = 0
    while t < T:
        k1  = f(r, t, w_q)                 
        k11 = g(d, r, t, w_q)

        k2_  = f(r + 0.5 * k1, t, w_q)
        k22 = g(d + 0.5 * k11, r + 0.5 * k1, t, w_q)

        k3  = f(r + 0.5 * k2_, t, w_q)
        k33 = g(d + 0.5 * k22, r + 0.5 * k2_, t, w_q)

        k4  = f(r + k3, t, w_q)
        k44 = g(d + k33, r + k3, t, w_q)

        r  += h * (k1  + 2 * k2_  + 2 * k3  + k4)  / 6
        d  += h * (k11 + 2 * k22 + 2 * k33 + k44) / 6
        x.append(r[0])
        y.append(r[1])
        z.append(r[2])
        time.append(t)

        orth_1 = d[0]
        l1 += np.log(norm(orth_1))
        d[0] = orth_1 / norm(orth_1)

        orth_2 = d[1] - np.dot(d[1], d[0]) * d[0]
        l2 += np.log(norm(orth_2))
        d[1] = orth_2 / norm(orth_2)

        orth_3 = d[2] - (np.dot(d[2], d[1]) * d[1]) - (np.dot(d[2], d[0]) * d[0]) 
        l3 += np.log(norm(orth_3))
        d[2] = orth_3 / norm(orth_3)
        t += h

    lya1 = l1 / (T) / np.log(2)
    lya2 = l2 / (T) / np.log(2)
    lya3 = l3 / (T) / np.log(2)
    print("l1 = ", lya1,"l2 = ", lya2,"l3 = ", lya3)

def RkLyapunov(loc_max, w_time, lyapun1, lyapun2, lyapun3, w_q):
    # Initial conditions
    d = np.array([[1,0,0], [0,1,0], [0,0,1]], float)
    r = np.array([x0, y0, z0], float)

    l1, l2, l3 = 0, 0, 0
    x, y, z  = [], [], []

    x.append(r[0])
    y.append(r[1])
    z.append(r[2])

    t = 0
    i = 0
    while t < T:
        k1  = f(r, t, w_q)                 
        k11 = g(d, r, t, w_q)

        k2_  = f(r + 0.5 * k1, t, w_q)
        k22 = g(d + 0.5 * k11, r + 0.5 * k1, t, w_q)

        k3  = f(r + 0.5 * k2_, t, w_q)
        k33 = g(d + 0.5 * k22, r + 0.5 * k2_, t, w_q)

        k4  = f(r + k3, t, w_q)
        k44 = g(d + k33, r + k3, t, w_q)

        r  += h * (k1  + 2 * k2_  + 2 * k3  + k4)  / 6
        d  += h * (k11 + 2 * k22 + 2 * k33 + k44) / 6
        x.append(r[0])
        y.append(r[1])
        z.append(r[2])

        orth_1 = d[0]
        l1 += np.log(norm(orth_1))
        d[0] = orth_1 / norm(orth_1)

        orth_2 = d[1] - np.dot(d[1], d[0]) * d[0]
        l2 += np.log(norm(orth_2))
        d[1] = orth_2 / norm(orth_2)

        orth_3 = d[2] - (np.dot(d[2], d[1]) * d[1]) - (np.dot(d[2], d[0]) * d[0]) 
        l3 += np.log(norm(orth_3))
        d[2] = orth_3 / norm(orth_3)

        if ((t > t_porog) and ((z[i - 2] - z[i - 1]) * (z[i - 1] - z[i]) < 0) and ((z[i - 2] - z[i - 1]) < 0)):
            loc_max.append(z[i-1])
            w_time.append(w_q)
        t += h
        i += 1

    lya1 = l1 / (T) / np.log(2)
    lya2 = l2 / (T) / np.log(2)
    lya3 = l3 / (T) / np.log(2)
    lyapun1.append(lya1)
    lyapun2.append(lya2)
    lyapun3.append(lya3)

def localMaxLyap():
    loc_max = list()
    w_time  = list()
    lyapun1 = list()
    lyapun2 = list()
    lyapun3 = list()
    l_time = list()

    while i_py <= i_end:
        RkLyapunov(loc_max, w_time, lyapun1, lyapun2, lyapun3, i_py)
        l_time.append(i_py)
        i_py += 0.00001
        if i_py == 0.005:
            print("POLOVINA")


    print("record to the file")
    print("begin")
    with open('loc_max', 'wb') as fp:
        pickle.dump(loc_max, fp)

    with open('w_time', 'wb') as fp:
        pickle.dump(w_time, fp)

    with open('lyapun1', 'wb') as fp:
        pickle.dump(lyapun1, fp)

    with open('lyapun2', 'wb') as fp:
        pickle.dump(lyapun2, fp)

    with open('lyapun3', 'wb') as fp:
        pickle.dump(lyapun3, fp)

    with open('l_time', 'wb') as fp:
        pickle.dump(l_time, fp)
    print("end")

    fig = plt.figure(figsize=(14,9))
    ax1 = fig.add_subplot(1,1,1)
    ax1.scatter(w_time, loc_max, s = 1.5, color= 'darkred')
    ## ax1.set_xlabel('r')
    ax2 = ax1.twinx()
    ax2.scatter(l_time,lyapun1,s = 3, color= 'blue', label = 'l1')#s = 3, color= 'blue', 
    ax2.scatter(l_time,lyapun2,s = 3, color= 'orange',  label = 'l2')#s = 3, color= 'orange',
    ax2.scatter(l_time,lyapun3,s = 3, color= 'green',  label = 'l2')
    ##ax1.set_xlim(0.0045, 0.0079)
    ##ax1.set_ylim(-0.08, 5.5)
    ax2.grid('on')
    ax1.xlabel('ω$_s$')
    ax1.ylabel('R$_m$$_a$$_x$') 
    plt.show()

def Realisation():
    rX = list()
    rY = list()
    rZ = list()
    rT = list()
    w_q = 0.004
    RkLyapunovRealis(rX, rY, rZ, rT, w_q)
    print("x = ", len(rX), "y = ", len(rY), "z = ", len(rZ))
    title_name = str("Q$_s$ = " + str(Qs) + " Q$_0$ = " + str(Q0) + " θ$_r$ = " + str(q_r) + " k$_r$ = " + str(k_r) + " ω$_s$ = " + str(w_q))
    plt.title(title_name)
    plt.plot(rT,rX, label = 'График fz(t)')
    plt.plot( rT,rY, label = 'График fp(t)')
    plt.plot( rT,rZ, label = 'График fr(t)')
    #plt.plot(rX,rZ)
    #plt.plot( rX,rY)
    plt.xlabel('Ось time ')
    plt.ylabel('Ось Z, P, R')
    #plt.xlabel('Z')
    #plt.ylabel('R')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    a_p = 0.001
    a_z = 0.0001
    alfa_0 = 0.23
    b_p = 0.01
    b_z = 0.01
    g_p = 0.001
    q_p = 6.
    q_z = 5.5
    Z0 = 0.
    Z1 = 1.
    P0 = 0.
    P1 = 1.
    k_z = 0.15
    k_p = 0.05
    Q0 = 1.8
    Qs = 3.2
    a_r = 0.01
    b_r = 0.01
    R0 = 2.
    R1 = 1.
    q_r = 5.6
    k_r = 0.1

    t_porog = 100000
    T  = 200000
    h = 0.05
    x0 = 0.2
    y0 = 0.2
    z0 = 0.2

    i_py = 0.002
    i_end = 0.008
    
    Realisation()
    #localMaxLyap()