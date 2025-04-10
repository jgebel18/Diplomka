import numpy as np
import matplotlib.pyplot as plt
cimport numpy as np
cimport cython
from libc.math cimport sqrt, pi

cdef double gamma_0 = 1
cdef double V = 0.35
cdef double Delta = 0.1
cdef double tolerance = 0.02
cdef int max_iterations = 300
cdef int max_iterations_2 = 300
cdef int n_max = 6000
cdef int m_max = 300


cdef double g_0_inverse(int n, double beta):
    return ((2 * n + 1) * pi) / sqrt(V * beta) + Delta * sqrt(beta / V) * np.sign(2 * n + 1)

cdef double g(int n, double beta):
    cdef double g_0_inv = g_0_inverse(n, beta)
    return 2 / (g_0_inv + np.sign(g_0_inv) * sqrt(g_0_inv**2 + 3 * gamma_0))

cdef double new_g(int n, double beta, double B):
    return 2 / (B + np.sign(B) * sqrt(B**2 + 3 * gamma_0))

cdef double B(int n, double beta, int max_m, dict tab_1, dict tab_2):
    return g_0_inverse(n, beta) + 0.75 * soucet_gamma(max_m, tab_1, tab_2, n)



cdef dict napln_tabulku_B(func, int max_n, int max_m, dict tab_1, dict tab_2, double beta):
    cdef dict table = {}
    for n in range(-max_n, max_n + 1):
        table[n] = func(n, beta, max_m, tab_1, tab_2)
    return table

cdef dict iteration_process_G(int iterace, double beta, double tol, int max_n, int max_m, dict tab_2):
    cdef list pole_tabulek = []
    cdef dict B_pole = {}
    cdef dict g_n = {n: g(n, beta) for n in range(-max_n, max_n + 1)}
    cdef np.ndarray[int, ndim=1] n_values = np.arange(-max_n, max_n + 1, dtype=np.int32)
    cdef int iteration
    cdef dict new_g_n
    cdef list g_n_values, new_g_n_values, deviation_values
    cdef double deviation
    cdef double hodnota=1.0
    pole_tabulek.append(g_n)
   
    for iteration in range(1, iterace + 1):
        B_pole = napln_tabulku_B(B, max_n, max_m, g_n, tab_2, beta)
        new_g_n = {n: new_g(n, beta, B_pole[n]) for n in range(-max_n, max_n + 1)}
        deviation_values = [abs(new_g_n[n] - g_n[n]) / g_n[n] for n in range(-max_n, max_n + 1)]
        deviation = max(deviation_values)
       
        pole_tabulek.append(g_n)
        g_n = new_g_n
        if deviation < tol:
            #vykresly_do_jednoho_G(beta, pole_tabulek,n_values,hodnota ,tolerance,iteration)
            print(iteration)
            break
           
            return g_n
    #vykresly_do_jednoho_G(beta, pole_tabulek,n_values,hodnota ,tolerance, iteration)        
            
    return g_n

cdef double soucet_gamma(int max_m, dict tabulka_1, dict tabulka_2, int n):
    cdef int m
    cdef double suma = 0
    for m in range(-max_m, max_m + 1):
        if (n + m) in tabulka_1 and m != 0:
            suma += tabulka_1[n + m] * tabulka_2[m]
    return suma

cdef double s(int m, int max_n, double beta):
    cdef int n
    cdef double suma = 0
    for n in range(-max_n, max_n + 1):
        suma += g_0_inverse(n, beta) * g_0_inverse(n + m, beta)
    return suma

cdef double s_2(int m, int max_n, int max_m, double beta, dict tab_2, double tol):
    cdef int n
    cdef double suma = 0
    cdef dict pole_1 = iteration_process_G(max_iterations, beta, tol, max_n, max_m, tab_2)
    for n in range(-max_n, max_n + 1):
        if ((n + m) < max_n + 1 and (n + m) >= -max_n):
            suma += pole_1[n] * pole_1[n + m]
    return suma

cdef double calculate_gamma(int m, int max_n, double beta):
    return s(m, max_n, beta) / (1 - s(m, max_n, beta))

cdef double gamma_2(int m, int max_n, int max_m, double beta, dict tab_2, double tol):
    cdef double s2_val = s_2(m, max_n, max_m, beta, tab_2, tol)
    return s2_val / (1 - s2_val)

cdef dict iteration_process_gamma(int max_n, int max_m, double beta, int iterace, np.ndarray[int, ndim=1] pole_m, double tol):
    cdef dict gamma_m_tab = {m: calculate_gamma(m, max_n, beta) for m in range(-max_m, max_m + 1)}
    cdef int iteration
    cdef dict new_gamma_m
    cdef list gamma_m_tab_values, new_gamma_m_values, deviation_values
    cdef double deviation
    
    for iteration in range(1, iterace + 1):
        new_gamma_m = {m: gamma_2(m, max_n, max_m, beta, gamma_m_tab, tol) for m in range(-max_m, max_m + 1)}
        deviation_values = [abs(new_gamma_m[m] - gamma_m_tab[m]) / gamma_m_tab[m] for m in range(-max_m, max_m + 1)]
        deviation = max(deviation_values)
        gamma_m_tab = new_gamma_m
        if deviation < tol:
            break
            return gamma_m_tab
            
          
    return gamma_m_tab

def vykresly(double beta, dict tabulka, np.ndarray[int, ndim=1] pole_promennych, double hodnota, double tol):
    cdef list nahradni = [tabulka[n] / hodnota for n in pole_promennych]
    plt.figure(figsize=(10, 4))
    plt.errorbar(pole_promennych, nahradni, [tol * abs(nahradni[i]) for i in range(len(nahradni))], label=f'gamma(m)/gamma(0), β = {beta}', marker='o')
    plt.xlabel('$m$')
    plt.ylabel('$γ(m)/γ(0)$')
    plt.title('')
    plt.legend()
    plt.show()

def vykresly_do_jednoho_G(double beta, list tabulky, np.ndarray[int, ndim=1] pole_promennych, double hodnoty, double tol, int pocet_iteraci):
    plt.figure(figsize=(10, 4))
    cdef list nahradni
    for i in range(0,pocet_iteraci+1):
        nahradni = [tabulky[i][n] / hodnoty for n in pole_promennych]
        plt.errorbar(pole_promennych, nahradni, [tol * abs(nahradni[n]) for n in range(len(nahradni))], label=f'{i}-tá iterace, β = {beta}', marker='o')
    plt.xlabel('$n$')
    plt.ylabel(' $g(n)$')
    plt.title('')
    plt.legend()
    plt.show()

def main():
    cdef list beta_values = [50000]
    cdef double beta
    cdef np.ndarray[int, ndim=1] m_values = np.arange(-m_max, m_max + 1, dtype=np.int32)
    for beta in beta_values:
        gamma_history = iteration_process_gamma(n_max, m_max, beta, max_iterations_2, m_values, tolerance)
        vykresly(beta, gamma_history, m_values, gamma_0, tolerance)
        break

if __name__ == "__main__": 

    main()