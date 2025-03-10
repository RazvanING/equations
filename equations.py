#Required Solar Panel Area Calculation
def calculate_solar_panel_area(F_T, I_sp, g_0, zeta, eta_T, S_0, eta_p, eta_as):
    return (F_T * I_sp * g_0 * zeta) / (2 * eta_T * S_0 * eta_p * eta_as)
#Spacecraft Side Area Calculation
def calculate_side_area(A_Rs, A_i):
    return 4 * A_Rs * A_i
#Required Array Area Calculation
def calculate_array_area(A_p, A_Rs, A_i):
    return (A_p - 4 * A_Rs * A_i) / math.pi
#Molecular speed calculation
def molecular_speed_ratio(u_inf, R, T_inf):
    return u_inf / math.sqrt(2 * R * T_inf)
#Drag Coefficient Calculation (Flat Plate in Free-Molecular Flow)
def drag_coefficient(B, epsilon, alpha, S_inf, T_r, T_inf):
    term1 = (B * (1 - epsilon * math.cos(2 * alpha))) / (math.sqrt(math.pi) * S_inf)
    term2 = math.exp(-S_inf**2 * math.sin(alpha)**2)
    term3 = math.sin(alpha) / S_inf**2
    term4 = 1 + 2 * S_inf**2 + epsilon * (1 - 2 * S_inf**2 * math.cos(2 * alpha))
    term5 = erf(S_inf * math.sin(alpha)) + (1 - epsilon) * (math.sin(alpha)**2) / (S_inf * math.sqrt(math.pi))
    correction = math.sqrt(T_r / T_inf)
    return term1 * (term2 + term3 * term4 * term5) * correction
#Intake Compression Ratio Calculation
def intake_compression_ratio(beta_0, A_Ri, T_w, eta_c):
    return 0.87 * beta_0 * (0.244 + 0.33 * math.log(A_Ri)) * ((189 / (T_w + 33)) + 0.435) * (1 - 1.625 * eta_c)
#Force Balance (Drag-Compensation Equation of Motion)
def drag_compensation(F_T, m_dot_t, u_inf, F_D_i, F_D_s, F_D_a):
    return F_T - m_dot_t * u_inf - (F_D_i + F_D_s + F_D_a)
#Drag-Compensation in Terms of I_sp and eta_T
def drag_compensation_quadratic(g_0, C_D_a, rho_inf, u_inf, zeta, a_P, eta_as, eta_T, I_sp):
    a = (g_0**2 * C_D_a * rho_inf * u_inf**2 * zeta) / (4 * a_P * eta_as * eta_T)
    b = -g_0
    c = u_inf * (1 + ...)
    return a * I_sp**2 + b * I_sp + c
#Thruster controll law
def thruster_control_law(r, r_p, r_a, C_T):
    return (abs(r) - r_p) / (r_a - r_p) > C_T
