from sympy import diff, symbols, cos, sin, exp, simplify

a_p, a_z, a_r, alfa_0, b_p, b_z, b_r, g_p, q_p, q_z, q_r, Z0, Z1, P0, P1, R0, R1, k_z, k_p, k_r, Q0, Qs, t, t0, w_q, Z, P, R = symbols('a_p a_z a_r alfa_0 b_p b_z b_r g_p q_p q_z q_r Z0 Z1 P0 P1 R0 R1 k_z k_p k_r Q0 Qs t t0 w_q Z P R')
fx = (-(a_z + g_p * P) * Z + b_z * (Z0 - (Z0 - Z1) / (1 + exp(-(Q0 + alfa_0 * Z + Qs * sin(w_q * t) - q_z) / k_z))))
fy = (-a_p * P + b_p * (P0 - (P0 - P1) / (1 + exp(-(Q0 + alfa_0 * Z + Qs * sin(w_q * t) - q_p) / k_p))))
fz = (-a_r * R + b_r * (R0 - (R0 - R1) / (1 + exp(-(Q0 + alfa_0 * Z + Qs * sin(w_q * t) - q_r) / k_r))))

   

fx_difZ = diff(fx, Z, 1)
fx_difP = diff(fx, P, 1)
fx_difR = diff(fx, R, 1)

# simplify(fx_difZ) 
# simplify(fx_difP) 

print("fx_x =", fx_difZ)
print("fx_y =", fx_difP)
print("fx_z =", fx_difR)

fy_difZ = diff(fy, Z, 1)
fy_difP = diff(fy, P, 1)
fy_difR = diff(fy, R, 1)

# simplify(fy_difZ) 
# simplify(fy_difP) 

print("fy_x = ", fy_difZ)
print("fy_y = ", fy_difP)
print("fy_z = ", fy_difR)

# simplify(fy_difZ) 
# simplify(fy_difP) 

fz_difZ = diff(fz, Z, 1)
fz_difP = diff(fz, P, 1)
fz_difR = diff(fz, R, 1)

print("fz_x = ", fz_difZ)
print("fz_y = ", fz_difP)
print("fz_z = ", fz_difR)


