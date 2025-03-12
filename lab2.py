#Laborator 2

#############################################################
#Aplicatie in consola scrisa in limbajul de Programare Python
#############################################################

#Sa se proiecteze un compresor cu gradul de comprimare de 4.15, debit 20 kg/s, presiune normala si temperatura 288 K. Functionare in conditii statice la sol.

#1) Determinarea vitezei de rotatie In principiu, se estimeaza o viteza maxima la varful paletei de rotor de la prima treapta. Presupunand ca nu exista dispozitiv de admisie, deci intrarea este axiala. Pentru cazul cand nu exista trepte de ventilator, viteza maxima la varful paletei este in jur de 350 m/s iar viteza axiala poate varia intre 150 si 200 m/s pentru a preveni cresterea numarului Mv. Raportul razei la baza si razei la varf este cuprins intre 0.4 si 0.6. Presupunem ca nu avem pierderi la intrare.

# Importam biblioteca math
import math
import numpy
from scipy.optimize import fsolve

# Biblioteca math contine functii matematice, cum ar fi sin, cos, tan, sqrt, log, etc.

# Constante:
R = 287 #J/kgK
gamma = 1.4
Cp = 1004 #J/kgK
Cv = 717 #J/kgK

#Declaram si solicitam datele de intrare
#Se declara ca fiind float, pentru a putea realiza operatii cu numere reale
Pfr_1 = 1.013 #float(input("Introduceti presiunea de intrare in bar: "))
#^ nume variabila, = operator de atribuire, float() functie de conversie a datelor introduse in float
Tfr_1 = 288 #float(input("Introduceti temperatura de intrare in Kelvin: "))
M_0 = 0 #float(input("Introduceti nr Mach: "))
C_1 = 150 #float(input("Introduceti viteza de intrare in m/s: "))
C_ai = 150 #float(input("Introduceti viteza axiala de intrare in m/s: "))
C_u1 = 0 #float(input("Introduceti viteza circumferentiala de intrare in m/s (default 0): ") or 0)
#^ daca nu se introduce nimic, se va lua valoarea 0
m_dot = 20#float(input("Introduceti debitul de aer in kg/s: "))
#^ debitul de aer
pi_c = 4.15 #float(input("Introduceti raportul de compresie al treptei: "))
#^ raportul de compresie dorit

T_1 = Tfr_1 - C_1**2 / (2 * Cp)
#^ calculam temperatura statica la intrare
print("Temperatura statica la intrare este: ", T_1, "K")
#^ afisam in consola temperatura statica la intrare
P_1 = Pfr_1 * (T_1 / Tfr_1)**(gamma / (gamma - 1))
#^ calculam presiunea statica la intrare 
print("Presiunea statica la intrare este: ", P_1, "bar")
#^ afisam in consola presiunea statica la intrare
rho_1 = P_1 * 10**5 / (R * T_1)
#^ calculam densitatea la intrare
print("Densitatea la intrare este: ", rho_1, "kg/m^3")
#^ afisam in consola densitatea la intrare

constant = m_dot / (math.pi * rho_1 * C_ai)


# Function to solve for rv
def equation(rv, rv_ratio):
    rv_prime = rv / rv_ratio
    return rv**2 - constant / (1 - (rv/rv_prime)**2)

# Calculate rv for different rv/rv' ratios
rv_ratios = [0.4, 0.45, 0.5, 0.55, 0.6]
results = []

for ratio in rv_ratios:
    # Initial guess for rv
    rv_guess = 0.2
    
    # Solve for rv
    rv = fsolve(lambda x: equation(x, ratio), rv_guess)[0]
    
    # Calculate angular velocity
    U_v = C_ai / ratio  # peripheral velocity
    omega = (U_v / rv) / (2 * math.pi)  
    
    
    results.append((ratio, rv, omega))



# Print the table
print("rv/rv' | rv (m)  | ω (rot/s)")
print("-" * 30)
for ratio, rv, omega in results:
    print(f"{ratio:.2f}  | {rv:.4f} | {omega:.1f}")

# Alegerea corecta tine seama de caracteristicile turbinei care antreneaza compresorul. De pilda, diametrul exterior al turbineiar trebui sa fie comparabil cu cel al compresorului. Stiind ca acest diametru este de 0.239 m, rezulta ca a alegere corecta ar fi raportul razelor de 0.5 pentru care se obtine o viteza unghiulara de 246.3 rot/s de unde rezulta viteza de transport la varf de

rap_final = float(input("Introduceti raportul razelor:  "))
#^ introducem raportul razelor
r_v_final = float(input("Introduceti raza la varf in m: "))
#^ introducem raza la varf
w_final = float(input("Introduceti viteza in rot/s: "))
#^ introducem viteza unghiulara

U_v_final = 2 * math.pi * r_v_final * w_final
print("Viteza de transport la varf este: ", U_v_final, "m/s")
#^ afisam in consola viteza de transport la varf

#Se verifica numarul mach la varf
M_v = math.sqrt(C_1**2 + U_v_final**2) / math.sqrt(gamma * R * T_1)
print("Numarul Mach la varf este: ", M_v)
#^ afisam in consola numarul Mach la varf

r_b = r_v_final * rap_final
print("Raza la baza este: ", r_b, "m")
#^ afisam in consola raza la baza

#Calculam Aria
A_1 = math.pi *( r_v_final**2 - r_b**2)
print("Aria este: ", A_1, "m^2")

# Presiunea si temperatura la iesire din compresor sunt:
P_2_fr = pi_c * Pfr_1 
print("Presiunea la iesire este: ", P_2_fr, "bar")
T_2_fr = Tfr_1 * (P_2_fr / Pfr_1)**((gamma - 1) / gamma)
print("Temperatura la iesire este: ", T_2_fr, "K")

# S-a considerat un randament politropic de 0.9. Daca la iesire se considera viteza axiala pura tot de 150 m/s, atunci

T_2 = T_2_fr - C_ai**2 / (2 * Cp)
print("Temperatura statica la iesire este: ", T_2, "K")
P_2 = P_2_fr * (T_2 / T_2_fr)**(gamma / (gamma - 1))
print("Presiunea statica la iesire este: ", P_2, "bar")
rho_2 = P_2 * 10**5 / (R * T_2)
print("Densitatea la iesire este: ", rho_2, "kg/m^3")

# Se calculeaza aria la iesire
A_2 = m_dot / (rho_2 * C_ai)
print("Aria la iesire este: ", A_2, "m^2")

# Presupunand ca avem o configuratie cu raza medie constanta
r_m = (r_v_final + r_b) / 2
print("Raza medie este: ", r_m, "m")

h_2 = A_2 / (2 * math.pi * r_m)
print("Inaltimea este: ", h_2, "m")
r_2v = r_m + h_2/2
print("Raza la varf este: ", r_2v, "m")
r_2b = r_m - h_2/2
print("Raza la baza este: ", r_2b, "m")