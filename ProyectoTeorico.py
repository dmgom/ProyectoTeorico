# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 22:34:24 2024

@author: mauri
"""
#%%
# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.integrate import odeint

# # Definir la función del modelo ODE (puedes ajustar la ecuación)
# def modelo(y, t, params):
#     S0, P = y
#     k1, k2 = params  # Parametros a ajustar
#     dS0_dt = -k1 * S0 / (k2 + S0)
#     dP_dt = k1 * S0 / (k2 + S0) - P
#     return [dS0_dt, dP_dt]

# # Definir función para calcular chi^2 (error cuadrático entre datos observados y modelo)
# def chi2(model, data, sigma):
#     return np.sum(((data - model) / sigma) ** 2)

# # Implementar el algoritmo de Metropolis-Hastings
# def metropolis_hastings(data, time_points, initial_params, num_iterations, sigma, step_size=0.1):
#     # Inicializar parámetros y arrays para guardar resultados
#     params_current = np.array(initial_params)
#     chi2_current = np.inf
#     accepted_params = []
    
#     for i in range(num_iterations):
#         # Proponer nuevos parámetros
#         params_proposed = params_current + np.random.normal(0, step_size, size=len(params_current))
        
#         # Resolver las ODEs con Runge-Kutta (odeint)
#         solution_proposed = odeint(modelo, [data[0], 0], time_points, args=(params_proposed,))
#         P_proposed = solution_proposed[:, 1]
        
#         # Calcular chi2 para los parámetros propuestos
#         chi2_proposed = chi2(P_proposed, data, sigma)
        
#         # Aceptar o rechazar los nuevos parámetros
#         if chi2_proposed < chi2_current:
#             params_current = params_proposed
#             chi2_current = chi2_proposed
#             accepted_params.append(params_current)
#         else:
#             alpha = np.exp(-(chi2_proposed - chi2_current))
#             if np.random.rand() <= alpha:
#                 params_current = params_proposed
#                 chi2_current = chi2_proposed
#                 accepted_params.append(params_current)
                
#     return np.array(accepted_params)

# # Datos observados (ficticios) y puntos en el tiempo
# time_points = np.linspace(0, 5, 100)
# data_observada = np.sin(time_points)  # AQUI VAN DATOS REALES
# sigma = 0.1  # Incertidumbre (asumida constante)

# # Parámetros iniciales
# initial_params = [1.0, 1.0]  # Inicializa los parámetros

# # Número de iteraciones del algoritmo
# num_iterations = 1000

# # Ejecutar el algoritmo de Metropolis-Hastings
# resultados = metropolis_hastings(data_observada, time_points, initial_params, num_iterations, sigma)

# # Graficar los parámetros aceptados durante las iteraciones
# plt.plot(resultados[:, 0], label='Parámetro k1')
# plt.plot(resultados[:, 1], label='Parámetro k2')
# plt.legend()
# plt.xlabel('Iteración')
# plt.ylabel('Valor de parámetro')
# plt.title('Evolución de parámetros aceptados')
# plt.show()


# #%%
# import numpy as np
# import matplotlib.pyplot as plt
# from scipy.integrate import solve_ivp

# # Definir la función para P_bar
# def P_bar(P, K_I, K_II):
#     return (P*K_I + K_I*K_II*P**2) / (1 + P*K_I + K_I*K_II*P**2)

# # Definir la ODE para ddes/dt
# def ddes_dt(t, des, nu, gamma, P, K_I, K_II):
#     P_bar_val = P_bar(P, K_I, K_II)
#     return nu * P_bar_val - gamma * des

# # Parámetros iniciales
# nu = 1.0     # Tasa de transcripción
# gamma = 0.5  # Coeficiente de dilución
# K_I = 1.0    # Constante K_I
# K_II = 2.0   # Constante K_II
# P = 1.0      # Valor de P
# des_0 = 0.0  # Valor inicial de des (en t=0)
# t_span = (0, 10)  # Intervalo de tiempo para la solución
# t_eval = np.linspace(t_span[0], t_span[1], 200)  # Puntos de evaluación

# # Resolver la ODE
# sol = solve_ivp(ddes_dt, t_span, [des_0], args=(nu, gamma, P, K_I, K_II), t_eval=t_eval)

# # Graficar la solución
# plt.figure(figsize=(8,6), dpi=80)
# plt.plot(sol.t, sol.y[0], label=r'$des$', color='b')
# plt.xlabel("Tiempo (t)", fontsize=14)
# plt.ylabel(r"$des$", fontsize=14)
# plt.title(r"Evolución temporal de $des$", fontsize=16)
# plt.legend()
# plt.show()

#%%
##parametros tabla
alpha=37.9245
beta=49.5627
KI=12.5548
KII=22.8607
S0=80.7903

##equation 5.13
import numpy as np
import matplotlib.pyplot as plt

## función del modelo ODE 5.13
def func_prime(t, p, alpha, beta, S0, KI, KII):
    return (alpha * (S0 - p) / (KI + S0 - p)) - (beta * p / (KII + p))# + q*noise

# Resolver la ODE utilizando el método de Runge-Kutta de cuarto orden
# Masa molar de la proteína en g/mol
M = 71.834  # g/mol

# Volumen de la solución en litros (45 μL = 45e-6 L)
V = 45e-6  # L


def RK(alpha, beta, S0, KI, KII, p0, t0):
    h=0.02
    n_points = int((5+h)/h)  # Número de puntos de tiempo
    t = np.zeros(n_points)
    p = np.zeros(n_points)
    
    # Inicialización
    t[0] = t0
    p[0] = p0
    
    # Iteración sobre cada paso de tiempo
    for i in range(1, n_points):
        K1 = func_prime(t[i-1], p[i-1], alpha, beta, S0, KI, KII)
        t1 = t[i-1] + (h/2.0)
        p1 = p[i-1] + (h/2.0)*K1
        
        K2 = func_prime(t1, p1, alpha, beta, S0, KI, KII)
        t2 = t[i-1] + (h/2.0)
        p2 = p[i-1] + (h/2.0)*K2
        
        K3 = func_prime(t2, p2, alpha, beta, S0, KI, KII)
        t3 = t[i-1] + h
        p3 = p[i-1] + h*K3
        
        K4 = func_prime(t3, p3, alpha, beta, S0, KI, KII)
        av_K = (1.0/6.0) * (K1 + 2.0*K2 + 2.0*K3 + K4)
        
        # Actualizar las variables
        t[i] = t[i-1] + h
        p[i] = p[i-1] + h * av_K
    
    return p, t

# Parámetros iniciales
# alpha = 1.0
# beta = 0.5
# S0 = 10.0
# KI = 1.0
# KII = 2.0
p0 = 0.0
t0 = 0.0
#t_end = 1000.0

# Resolver la ODE
p_sol, t_sol = RK(alpha, beta, S0, KI, KII, p0, t0)

# Convertir los valores de P(t) de U.A./μg a moles utilizando la fórmula
P_moles = (p_sol * V) / M *200 # Convertir a moles
# Graficar la solución computacional
plt.figure(figsize=(14,12), dpi=80)
plt.plot(t_sol, P_moles, label="p(t)", color='b')
plt.xlabel("Time (minutes)", fontsize=14)
plt.ylabel("p (U.A./μg)", fontsize=14)
plt.title("ODE dp/dt con RK4", fontsize=16)
plt.legend()
plt.show()

#%%
#ESTIMACION CONSTANTES DE CI FORMULA
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# Datos proporcionados
data = [
    [1.42, 0],
    [28.3, 0.02],
    [39.62, 0.03],
    [63.21, 0.06],
    [88.68, 0.08],
    [100.94, 0.1],
    [122.64, 0.12],
    [155.66, 0.16],
    [185.85, 0.18],
    [210.38, 0.2],
    [240.57, 0.21],
    [266.98, 0.22],
    [304.72, 0.22],
    [337.74, 0.21],
    [383.02, 0.19],
    [430.19, 0.15],
    [477.36, 0.13],
    [507.55, 0.11],
    [546.23, 0.08],
    [575.94, 0.07]
]

# Convertir a arrays de NumPy
data_array = np.array(data)
P_en_moles = data_array[:, 0] * (10e-6)  # Primer valor de cada par (P), convertimos a moles
CI_en_p = data_array[:, 1]  # Segundo valor de cada par (CI)

# Función C_I_bar para calcular la concentración de I en función de P, KI y KII
def C_I_bar(P, KI, KII):
    """Modelo para calcular CI en función de P, KI y KII"""
    return (P * KI) / (1 + P * KI + KI * KII * (P ** 2))

# Ajuste de parámetros usando mínimos cuadrados (curve fitting)
# Inicializamos KI y KII con valores razonables para la optimización
initial_guess = [0.1, 0.05]  # Suponemos que KI y KII están cerca de estos valores

# Usamos curve_fit para ajustar el modelo C_I_bar a los datos observados
params, covariance = curve_fit(C_I_bar, P_en_moles, CI_en_p, p0=initial_guess)

# Los parámetros estimados son KI y KII
KI_est, KII_est = params
print(f"KI estimado: {KI_est:.4f}")
print(f"KII estimado: {KII_est:.4f}")

# Calcular los valores ajustados de CI usando los parámetros estimados
CI_fit = C_I_bar(P_en_moles, KI_est, KII_est)

# Graficar CI en función de P
plt.figure(figsize=(8, 6))
plt.plot(P_en_moles, CI_fit, label=f"Modelo ajustado: KI={KI_est:.4f}, KII={KII_est:.4f}", color='b', lw=2)
plt.scatter(P_en_moles, CI_en_p, color='r', label="Datos observados", zorder=5)
plt.xlabel('P (mol)', fontsize=12)
plt.ylabel('CI (Proporcion de I)', fontsize=12)
plt.title('Proporción de CI en función de P', fontsize=14)
plt.legend()
plt.grid(True)
plt.show()
##RESULTADOS
# KI=163.7011
# KII=935.0452

#%%
#modelo con ruido, calculando q como sqrt(<f(x)>-<g(x)>)
P_values=p_sol

f_values = alpha * (S0 - np.array(P_values)) / (KI + S0 - np.array(P_values))
g_values = beta * np.array(P_values) / (KII + np.array(P_values))

# Calcular los promedios
f_bar = np.mean(f_values)
g_bar = np.mean(g_values)

# Calcular q
q = np.sqrt(f_bar + g_bar)
print("Coeficiente q:", q)

def func_prime_noise(t, p, alpha, beta, S0, KI, KII, q):
    noise = np.random.normal(0, 0.01)
    return (alpha * (S0 - p) / (KI + S0 - p)) - (beta * p / (KII + p)) + q*noise

# Resolver la ODE utilizando el método de Runge-Kutta de cuarto orden
def RK_noise(alpha, beta, S0, KI, KII, p0, t0, q):
    h=0.02
    n_points = int((5+h)/h)  # Número de puntos de tiempo
    t = np.zeros(n_points)
    p = np.zeros(n_points)
    
    # Inicialización
    t[0] = t0
    p[0] = p0
    
    # Iteración sobre cada paso de tiempo
    for i in range(1, n_points):
        K1 = func_prime_noise(t[i-1], p[i-1], alpha, beta, S0, KI, KII, q)
        t1 = t[i-1] + (h/2.0)
        p1 = p[i-1] + (h/2.0)*K1
        
        K2 = func_prime_noise(t1, p1, alpha, beta, S0, KI, KII, q)
        t2 = t[i-1] + (h/2.0)
        p2 = p[i-1] + (h/2.0)*K2
        
        K3 = func_prime_noise(t2, p2, alpha, beta, S0, KI, KII, q)
        t3 = t[i-1] + h
        p3 = p[i-1] + h*K3
        
        K4 = func_prime_noise(t3, p3, alpha, beta, S0, KI, KII, q)
        av_K = (1.0/6.0) * (K1 + 2.0*K2 + 2.0*K3 + K4)
        
        # Actualizar las variables
        t[i] = t[i-1] + h
        p[i] = p[i-1] + h * av_K
    
    return p, t


p_sol, t_sol = RK_noise(alpha, beta, S0, KI, KII, p0, t0,q)

# Graficar la solución computacional
plt.figure(figsize=(8,6), dpi=80)
plt.plot(t_sol, p_sol, label="p(t)", color='b')
plt.xlabel("time (min)", fontsize=14)
plt.ylabel("p(t) con ruido", fontsize=14)
plt.title("ODE 5.13 con RK", fontsize=16)
plt.legend()
plt.show()
#%%
##3 ruido aplicacion Gillespie
t_end = 5.0

def reaction_rates(p, alpha, beta, S0, KI, KII):
    rate_f = alpha * (S0 - p) / (KI + S0 - p)  # Tasa de producción
    rate_g = beta * p / (KII + p)              # Tasa de degradación
    return rate_f, rate_g

# Función para realizar la simulación de Gillespie con ruido
def gillespie_with_noise(alpha, beta, S0, KI, KII, p0, t0, t_end, q):
    t = [t0]
    p = [p0]
    
    while t[-1] < t_end:
        # Calcular las tasas de reacción
        rate_f, rate_g = reaction_rates(p[-1], alpha, beta, S0, KI, KII)
        
        # Calcular la tasa total
        total_rate = rate_f + rate_g
        
        # Determinar el tiempo hasta el siguiente evento
        if total_rate > 0:
            dt = np.random.exponential(1 / total_rate)
        else:
            break  # Termina la simulación si no hay eventos

        # Determinar el tipo de evento (producción o degradación)
        if np.random.rand() < rate_f / total_rate:
            dp = 1  # Incremento en P (evento de producción)
        else:
            dp = -1  # Decremento en P (evento de degradación)
        
        # Agregar el ruido gaussiano
        noise = q * np.random.normal(0, 0.01)
        dp += noise

        # Actualizar el tiempo y el estado
        t.append(t[-1] + dt)
        p.append(p[-1] + dp)
    
    return np.array(t), np.array(p)

# Ejecutar la simulación de Gillespie con ruido
t_sol, p_sol = gillespie_with_noise(alpha, beta, S0, KI, KII, p0, t0, t_end, q)

# Graficar la solución
plt.figure(figsize=(8, 6), dpi=80)
plt.step(t_sol, p_sol, label="p(t) con ruido", color='b')
plt.xlabel("Tiempo (minutos))", fontsize=14)
plt.ylabel("p(t) U.A./ug", fontsize=14)
plt.title("Simulación con ruido de Gillespie", fontsize=16)
plt.legend()
plt.show()
#%%
#estocastico con ruido para n poblacion de bacterias
nmax = 5  #n bacterias
imax = 100  #eventos máximo

t = np.zeros((imax, nmax))
p = np.zeros((imax, nmax))  
p[0, :] = 0.0               # condicion

#Gillespie cada bacteria
for n in range(nmax):
    for i in range(1, imax):
        #tasas de reacción cada evento
        s1 = alpha * (S0 - p[i-1, n]) / (KI + S0 - p[i-1, n])  #tasa de producción de p
        s2 = beta * p[i-1, n] / (KII + p[i-1, n])              #tasa de degradación de p
        stotal = s1 + s2                                       #tasa total

        # Determinar el tiempo hasta el siguiente evento
        t[i, n] = t[i-1, n] + (1 / stotal) * np.log(1 / np.random.rand())

        # evento q aleatorio
        Q = np.random.rand()
        if Q < s1 / stotal:
            #evento de produccióna
            p[i, n] = p[i-1, n] + 1
        else:
            #evento de degradación
            p[i, n] = max(0, p[i-1, n] - 1)  # Evitar valores negativos de p

#Graficar la evolución de p para cada bacteria en la población
plt.figure(figsize=(10, 6))
for n in range(nmax):        
    plt.step(t[:, n], p[:, n])
plt.xlabel("Tiempo (minutos)")
plt.ylabel("p(t)")
plt.title("Simulación de Gillespie para una población de bacterias")
plt.show()

#%%
##equation 5.20
import numpy as np
import matplotlib.pyplot as plt

# Definir la ecuación 5.20 para C_I_bar
def C_I_bar(P, KIm, KIIm):
    return (P*KIm) / (1 + P*KIm + KIm*KIIm*(P**2))

P_vals = np.linspace(0,0.1,1000)

KIm=163.7011
KIIm=935.0452
# Resolver la ecuación para diferentes valores de P
C_I_vals = C_I_bar(P_vals, KIm, KIIm)
print(type(C_I_vals))
# Graficar la solución
plt.figure(figsize=(8,6), dpi=80)
plt.plot(P_vals, C_I_vals, label=r'$\bar{C}_I$', color='b')
plt.xlabel("P", fontsize=14)
plt.ylabel(r"$\bar{C}_I$", fontsize=14)
plt.title(r"Solución de $\bar{C}_I$ en función de P", fontsize=16)
plt.legend()
plt.show()

#%%
# #ancho de membrana en funcion del tiempo
# #GPt= np.array(())
# GPt = -(1.09219240*10**-5)*t_sol+5.79056764

# plt.figure(figsize=(8,6), dpi=80)
# plt.plot(t_sol, GPt, label=r'$\bar{GP}_I$', color='b')
# plt.xlabel("tiempo s", fontsize=14)
# plt.ylabel(r"$\bar{GP}$", fontsize=14)
# plt.title(r"Solución de $\bar{GP}$ en función de t", fontsize=16)
# plt.legend()
# plt.show()
#%%
# #Relacionando P(t) con GP(t), al final queda P(GP)
# #depsejar t de GPt
# def t_GP(GPt):
#     return (5.79056764 - GPt) / (1.09219240 * 10**-5)

# t_val=t_GP(GPt)
# p_GP = np.interp(t_val, t_sol, p_sol)

# # Graficar p en función de GP(t)
# plt.figure(figsize=(8,6), dpi=80)
# plt.plot(GPt, p_GP, label="p(GP)", color='b')
# plt.xlabel(r"$\bar{GP}$", fontsize=14)
# plt.ylabel("p(GP)", fontsize=14)
# plt.title("p en función de GP", fontsize=16)
# plt.legend()
# plt.show()

#%%
T = np.array([11.16, 12.44, 15.28, 19.89, 24.94, 30.23, 34.77, 39.77, 44.52])
GP_T = np.array([0.3, 0.29, 0.27, 0.25, 0.23, 0.22, 0.21, 0.19, 0.18])
from sklearn.metrics import mean_squared_error, mean_absolute_error

model = np.poly1d(np.polyfit(T, GP_T, 2))
ynew = model(T)
plt.figure(figsize=(12,10))
plt.plot(T, GP_T, 'o', T, ynew, '-' , )
plt.ylabel('GP_T MM_40°C Glyc')
plt.xlabel('T (°C)')
plt.title('Interpolación anisotropía en función de T')
mse_pol = mean_squared_error(GP_T, ynew)
mae_poly = mean_absolute_error(GP_T, ynew)
print(mse_pol,mae_poly, model)
#modelo lineal -0.003452 x + 0.3272, 0.0016172839506172836 0.035308641975308655
#modelo cuadrado 5.453e-05 x^2 - 0.006409 x + 0.36
plt.figure()
plt.plot(T,5.453e-05*T**2 - 0.006409*T+ 0.36 )

plt.show()
#%%
from tqdm import tqdm
#modelo completo constantes sacadas de tesis Juanita #valores constantes originales
alpha=37.9245
beta=49.5627
KI=12.5548
KII=22.8607
S0=80.7903

#T arbitrario distribucion uniforme
T=np.random.uniform(13,54,2000)
GP_T_cuadratica=5.453e-05*T**2 - 0.006409*T+ 0.36
#si alpha aprox. GP

nmax = 2000  # Número de bacterias en la población
imax = 100  # Número máximo de eventos por simulación

# Inicialización de matrices
t = np.zeros((imax, nmax))  # Tiempo para cada bacteria
p = np.zeros((imax, nmax))  # Valor de p para cada bacteria
p[0, :] = 0.0           # Condición inicial de p
C_In= np.zeros((imax, nmax))

alpha_aprox=GP_T_cuadratica*164.89
#alpha_aprox=np.array([alpha]*len(GP_T_cuadratica))

        # p_sol, t_sol = RK(alpha, beta, S0, KI, KII, p0, t0)
        # P_values=p_sol
        
        # f_values = alpha * (S0 - np.array(P_values)) / (KI + S0 - np.array(P_values))
        # g_values = beta * np.array(P_values) / (KII + np.array(P_values))
        
        # # Calcular los promedios
        # f_bar = np.mean(f_values)
        # g_bar = np.mean(g_values)
        
        # # Calcular q
        # q = np.sqrt(f_bar + g_bar)

#valores de prueba
# beta = 122.2165
# KI = -3161.0182
# KII = -19.3591
# S0 = 246128.9891

# Simulación de Gillespie para cada bacteria
for n in tqdm(range(nmax)):
    for i in range(1, imax):
        # p_sol, t_sol = RK(alpha_aprox[n], beta, S0, KI, KII, p0, t0)
        # P_values=p_sol

        # f_values = alpha_aprox[n] * (S0 - np.array(P_values)) / (KI + S0 - np.array(P_values))
        # g_values = beta * np.array(P_values) / (KII + np.array(P_values))

        # # Calcular los promedios
        # f_bar = np.mean(f_values)
        # g_bar = np.mean(g_values)

        # # Calcular q
        # q = np.sqrt(f_bar + g_bar)
        
        # noise = q * np.random.normal(0, 0.1)
        # Calcular las tasas de reacción para cada evento
        s1 = alpha_aprox[n] * (S0 - p[i-1, n]) / (KI + S0 - p[i-1, n])# +q*noise # Tasa de producción de p
        s2 = beta * p[i-1, n] / (KII + p[i-1, n])#  +q*noise        # Tasa de degradación de p
        stotal = s1 + s2                                       # Tasa total

        # Determinar el tiempo hasta el siguiente evento
        t[i, n] = t[i-1, n] + (1 / stotal) * np.log(1 / np.random.rand())
        
        
        # Determinar el tipo de evento usando un número aleatorio
        Q = np.random.rand()
        #ser cuidadoso, inicialmente Q>s2/total
        if Q > s2 / stotal:
            # Evento de producción: p aumenta
            p[i, n] = p[i-1, n] + 1
        else:
            # Evento de degradación: p disminuye
            p[i, n] = max(0, p[i-1, n] - 1)  # Evitar valores negativos de p
#%%
# Graficar la evolución de p para cada bacteria en la población
plt.figure(figsize=(10, 6))
for n in range(nmax): 
    #print(len(t[:,n]), len(p[:,n]))    
    plt.step(t[:, n], p[:, n], where="post", label=f"Bacteria {n+1}")
plt.xlabel("Tiempo")
plt.ylabel("p(GP)")
plt.title("Gillespie para una población de bacterias")
plt.show()

#plt.figure(figsize=(8,6), dpi=80)
#%%
# Como hacer el promedio
promedio_p = np.mean(p, axis=1)

# Graficar el promedio
plt.figure(figsize=(10, 6))
plt.step(t[:, 0], promedio_p, where="post", label="Promedio de todas las bacterias")
plt.xlabel("Tiempo")
plt.ylabel("p(GP)")
plt.title("Promedio de Gillespie para una población de bacterias")
plt.legend()
plt.show()

#%%
#Ec. Complejo I
KIm=163.7011
KIIm=935.0452
M = 71.834  # g/mol
# Volumen de la solución en litros (45 μL = 45e-6 L)
V = 45e-6  # L


def C_I_bar(P, KI, KII):
    return (P*KI) / (1 + P*KI + KI*KII*(P**2))
C_I_n=[]

P_moles=( (p * V) /M)*100 # Convertir a moles #factor de 100, cambiar a 200
for n in range (nmax):
    C_I_n.append(C_I_bar(np.mean(P_moles[-11:,n]), KIm, KIIm))
#C_I_n=np.array((C_I_n))
#%%
C_In=C_I_n
# plt.figure()
# for n in range (nmax):
#     #P_vals = (np.mean(p[-30:,n])*V) / M
#     P_vals = (p*V) / M*10
#     # Resolver la ecuación para diferentes valores de P
#     C_I_vals = C_I_bar(P_vals, KIIm, KIm)
#     #C_In.append(C_I_vals)
#     #print(type(C_I_vals))
#     # Graficar la solución
#     plt.plot(P_vals, C_I_vals, 'o')
# plt.xlabel("P", fontsize=14)
# plt.ylabel(r"$\bar{C}_I$", fontsize=14)
# plt.title(r"Solución de $\bar{C}_I$ en función de P", fontsize=16)

#%%

#Definición de bins
C_I_values=np.array(C_In).flatten()

#El numero de bins es 41, pues se 1°C es lo que se puede diferenciar
num_bins = 41
T_bins = np.linspace(min(T), max(T), num_bins)
C_I_bins = np.linspace(min(C_I_values), max(C_I_values), num_bins)

T_repeated = T#

#Verificar las dimensiones antes del histograma
print(f"Longitud de T_repeated: {len(T_repeated)}")
print(f"Longitud de C_I_values: {len(C_I_values)}")

if len(T_repeated) != len(C_I_values):
    raise ValueError("Las longitudes de T_repeated y C_I_values no coinciden.")

#Histograma bidimensional para la distribución conjunta P(T, C_I)
joint_hist, T_edges, C_I_edges = np.histogram2d(
    T_repeated, C_I_values, bins=[T_bins, C_I_bins], density=True
)
joint_prob = joint_hist / np.sum(joint_hist)  # Normalizar para obtener probabilidad conjunta


# Calcular la información mutua
nonzero_idxs = joint_prob > 0  # Evitar ceros en el logaritmo
T_marginal = np.sum(joint_prob, axis=1)  # Marginal en T
C_I_marginal = np.sum(joint_prob, axis=0)  # Marginal en C_I

# Calcular los centros de los bins
T_bin_centers = (T_bins[:-1] + T_bins[1:]) / 2
C_I_bin_centers = (C_I_bins[:-1] + C_I_bins[1:]) / 2

# Probabilidades marginales
T_prob = T_marginal
C_I_prob = C_I_marginal

mutual_information = np.sum(
    joint_prob[nonzero_idxs] * 
    np.log2(joint_prob[nonzero_idxs] / (T_marginal[:, None] * C_I_marginal)[nonzero_idxs])
)

print(f"Información mutua entre T y C_I: {mutual_information:.4f} bits")

#graficar las distribuciones
plt.figure(figsize=(12, 6))

# P(T)
plt.subplot(1, 3, 1)
plt.bar(T_bin_centers, T_prob, width=np.diff(T_bins), edgecolor="k", alpha=0.7)
plt.title("Distribución de $T$")
plt.xlabel("$T$")
plt.ylabel("$P(T)$")

# P(C_I)
plt.subplot(1, 3, 2)
plt.bar(C_I_bin_centers, C_I_prob, width=np.diff(C_I_bins), edgecolor="k", alpha=0.7)
plt.title("Distribución de $C_I$")
plt.xlabel("$C_I$")
plt.ylabel("$P(C_I)$")

# P(T, C_I)
plt.subplot(1, 3, 3)
plt.imshow(joint_prob.T, origin="lower", aspect="auto", extent=[T_bins[0], T_bins[-1], C_I_bins[0], C_I_bins[-1]], cmap="viridis")
plt.colorbar(label="$P(T, C_I)$")
plt.title("Distribución conjunta $P(T, C_I)$")
plt.xlabel("$T$")
plt.ylabel("$C_I$")

plt.tight_layout()
plt.show()
#%%
import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm

alpha = 37.9245
beta = 49.5627
KI = 12.5548
KII = 22.8607
S0 = 80.7903

nmax = 2000
imax = 100


cantidad_datos = 100
valores_informacion = []

for iteration in tqdm(range(cantidad_datos), desc="Calculando información mutua"):
    # T arbitrario distribución uniforme
    T = np.random.uniform(13, 54, nmax)
    GP_T_cuadratica = 5.453e-05 * T**2 - 0.006409 * T + 0.36
    alpha_aprox = GP_T_cuadratica * 164.89
    t = np.zeros((imax, nmax))  # Tiempo para cada bacteria
    p = np.zeros((imax, nmax))  # Valores de p
    p[0, :] = 0.0  # Condición inicial

    for n in range(nmax):
        for i in range(1, imax):
            s1 = alpha_aprox[n] * (S0 - p[i - 1, n]) / (KI + S0 - p[i - 1, n])
            s2 = beta * p[i - 1, n] / (KII + p[i - 1, n])
            stotal = s1 + s2

            t[i, n] = t[i - 1, n] + (1 / stotal) * np.log(1 / np.random.rand())

            Q = np.random.rand()
            if Q > s2 / stotal:
                p[i, n] = p[i - 1, n] + 1
            else:
                p[i, n] = max(0, p[i - 1, n] - 1)

    # Calcular C_I
    KIm = 163.7011
    KIIm = 935.0452
    M = 71.834  # g/mol
    V = 45e-6  # L

    def C_I_bar(P, KI, KII):
        return (P * KI) / (1 + P * KI + KI * KII * (P**2))

    C_In = []
    P_moles = ((p * V) / M) * 100  # Convertir a moles
    for n in range(nmax):
        C_In.append(C_I_bar(np.mean(P_moles[-11:, n]), KIm, KIIm))

    # Definición de bins
    C_I_values = np.array(C_In).flatten()
    num_bins = 41
    T_bins = np.linspace(min(T), max(T), num_bins)
    C_I_bins = np.linspace(min(C_I_values), max(C_I_values), num_bins)

    # Histograma bidimensional
    joint_hist, T_edges, C_I_edges = np.histogram2d(
        T, C_I_values, bins=[T_bins, C_I_bins], density=True
    )
    joint_prob = joint_hist / np.sum(joint_hist)

    # Información mutua
    nonzero_idxs = joint_prob > 0
    T_marginal = np.sum(joint_prob, axis=1)
    C_I_marginal = np.sum(joint_prob, axis=0)

    mutual_information = np.sum(
        joint_prob[nonzero_idxs]
        * np.log2(joint_prob[nonzero_idxs] / (T_marginal[:, None] * C_I_marginal)[nonzero_idxs])
    )
    valores_informacion.append(mutual_information)

# Calcular promedio y desviación estándar
mean_mutual_information = np.mean(valores_informacion)
std_mutual_information = np.std(valores_informacion)

print(f"Promedio de información mutua: {mean_mutual_information:.4f} bits")
print(f"Desviación estándar: {std_mutual_information:.4f} bits")

# #Valores datos reales p(t)
# import numpy as np
# from scipy.integrate import solve_ivp
# from scipy.optimize import curve_fit
# import matplotlib.pyplot as plt

# # Datos experimentales
# t_exp = np.array([0.88, 0.94, 0.97, 1.03, 1.07, 1.12, 1.16, 1.23, 1.3, 1.37, 1.4, 1.46, 1.53, 1.63, 1.72, 1.81, 1.91, 2.02, 2.15, 2.25, 2.35, 2.5, 2.59, 2.71, 2.83, 2.97, 3.1, 3.27, 3.43, 3.6, 3.79, 3.92, 4.08, 4.18, 4.35, 4.47, 4.63, 4.81, 4.91, 5.04])
# p_exp = np.array([4.25, 5.98, 7.27, 8.55, 9.69, 10.93, 11.81, 13.05, 14.04, 15.03, 16.02, 16.81, 17.6, 18.98, 19.72, 20.96, 21.75, 22.39, 23.43, 24.27, 24.82, 25.81, 26.3, 26.5, 27.49, 27.88, 28.33, 28.82, 29.27, 30.06, 30.35, 30.7, 31.0, 31.34, 31.64, 31.79, 32.18, 32.28, 32.63, 33.12])

# # Función ODE
# def func_prime(t, p, alpha, beta, S0, KI, KII):
#     return (alpha * (S0 - p) / (KI + S0 - p)) - (beta * p / (KII + p))

# # Solución de la ODE
# def solve_ode(t, beta, KI, KII, S0, alpha=0.25):
#     p0 = p_exp[0]  # Condición inicial (primer valor experimental)
#     sol = solve_ivp(func_prime, [t[0], t[-1]], [p0], t_eval=t, args=(alpha, beta, S0, KI, KII))
#     return sol.y[0]  # Retornamos p(t)

# # Función objetivo para ajustar
# def model_to_fit(t, beta, KI, KII, S0):
#     return solve_ode(t, beta, KI, KII, S0)

# # Ajustar los parámetros
# initial_guess = [1.0, 1.0, 1.0, 35.0]  # Valores iniciales de beta, KI, KII, S0
# params, covariance = curve_fit(model_to_fit, t_exp, p_exp, p0=initial_guess)

# # Parámetros óptimos
# beta_opt, KI_opt, KII_opt, S0_opt = params

# # Resultados
# print(f"Parámetros ajustados:")
# print(f"beta = {beta_opt:.4f}")
# print(f"KI = {KI_opt:.4f}")
# print(f"KII = {KII_opt:.4f}")
# print(f"S0 = {S0_opt:.4f}")

# # Graficar resultados
# p_simulated = model_to_fit(t_exp, beta_opt, KI_opt, KII_opt, S0_opt)
# plt.plot(t_exp, p_exp, 'o', label='Datos experimentales')
# plt.plot(t_exp, p_simulated, '-', label='Modelo ajustado')
# plt.xlabel('Tiempo (t)')
# plt.ylabel('p(t)')
# plt.legend()
# plt.show()