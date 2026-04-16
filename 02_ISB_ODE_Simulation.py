"""
ISB Protocol: Thermodynamic Coupled ODE Integration & Mega-Simulation
---------------------------------------------------------------------
Evaluates the continuous temporal trajectory of modeled astrocytic 
bioenergetics across parameterized ancestral cohorts (N=30,000).
"""

import numpy as np
import pandas as pd
from scipy.integrate import odeint
import warnings

warnings.filterwarnings('ignore')

# ---------------------------------------------------------
# 1. BIOPHYSICAL CONSTANTS & BOUNDARY PARAMETERS
# ---------------------------------------------------------
V_MAX_BASE = 5.0   # Modeled maximal clearance velocity (mM/s)
K_ATP = 0.5        # Estimated Michaelis constant for ATP binding (mM)
K_M = 1.0          # Estimated affinity constant for extracellular glutamate (mM)
C_BASAL = 0.2      # Estimated minimum baseline ATP consumption
K_LEAK = 0.05      # Passive diffusion coefficient for mathematical tractability

# ---------------------------------------------------------
# 2. ODE SYSTEM DEFINITION
# ---------------------------------------------------------
def bioenergetic_system(state, t, PaCO2, R_stress, P_basal_genetic):
    """
    Computes the time derivatives for modeled astrocytic ATP (A) 
    and extracellular glutamate (G).
    """
    A, G = state
    
    # Boundary constraint preventing absolute zero arithmetic errors
    A, G = max(1e-5, A), max(1e-5, G)
    
    # Non-linear limitation functions (Vasoconstriction & Bohr Effect)
    # Constants strictly adhere to derived exact fractional constraints
    k_v = 7.0 / 200.0
    beta = 1.0 / 10.0
    
    phi_vaso = np.exp(k_v * (PaCO2 - 40.0))
    phi_bohr = 1.0 / (1.0 + 10**(beta * (40.0 - PaCO2)))
    
    # Modeled ATP Production Vector
    P_ATP_current = P_basal_genetic * phi_vaso * phi_bohr
    
    # ATP-Dependent Michaelis-Menten Clearance Kinetics
    V_GLT1 = V_MAX_BASE * (A / (K_ATP + A)) * (G / (K_M + G))
    
    # System of Equations
    dA_dt = P_ATP_current - C_BASAL - (4.0 * V_GLT1)
    dG_dt = R_stress - V_GLT1 - (K_LEAK * G)
    
    return [dA_dt, dG_dt]

# ---------------------------------------------------------
# 3. PAN-ANCESTRY MEGA-SIMULATION INITIALIZATION
# ---------------------------------------------------------
if __name__ == "__main__":
    print("Initializing in silico Pan-Ancestry cohort parameters...")
    np.random.seed(42) # Ensuring computational reproducibility
    N_PER_COHORT = 10000
    
    # Environmental Parameterization: Borderline Stress (Covert Hyperventilation)
    paco2_env = np.random.normal(39.2, 0.6, N_PER_COHORT) 
    
    # Genomic Parameterization: Baseline Capacity & Afferent Load
    # Cohort 1: European (EUR) Parameterization
    p_basal_eur = np.random.normal(8.0, 0.2, N_PER_COHORT)
    r_stress_eur = np.random.normal(0.6, 0.1, N_PER_COHORT)
    
    # Cohort 2: East Asian (EAS) Parameterization (Elevated R_stress)
    p_basal_eas = np.random.normal(8.0, 0.2, N_PER_COHORT)
    r_stress_eas = np.random.normal(1.0, 0.15, N_PER_COHORT) 
    
    # Cohort 3: African (AFR) Parameterization
    p_basal_afr = np.random.normal(8.4, 0.2, N_PER_COHORT)
    r_stress_afr = np.random.normal(0.6, 0.1, N_PER_COHORT)
    
    # Aggregate Synthetic Population Arrays
    pop_paco2 = np.concatenate([paco2_env, paco2_env, paco2_env])
    pop_pbasal = np.concatenate([p_basal_eur, p_basal_eas, p_basal_afr])
    pop_rstress = np.concatenate([r_stress_eur, r_stress_eas, r_stress_afr])
    pop_labels = np.array(['EUR']*N_PER_COHORT + ['EAS']*N_PER_COHORT + ['AFR']*N_PER_COHORT)
    
    # Apply physiological boundaries to input distributions
    pop_paco2 = np.clip(pop_paco2, 20.0, 45.0)
    pop_rstress = np.clip(pop_rstress, 0.1, 5.0)
    
    # ---------------------------------------------------------
    # 4. NUMERICAL INTEGRATION EXECUTION
    # ---------------------------------------------------------
    print(f"Executing numerical ODE integration for N={N_PER_COHORT*3} subjects. Please wait...")
    total_subjects = N_PER_COHORT * 3
    final_atp = np.zeros(total_subjects)
    final_glu = np.zeros(total_subjects)
    
    time_span = np.linspace(0, 100, 100)
    initial_state = [3.0, 0.2] # Initial homeostasis: [ATP_0, Glutamate_0]
    
    for i in range(total_subjects):
        solution = odeint(
            bioenergetic_system, 
            initial_state, 
            time_span, 
            args=(pop_paco2[i], pop_rstress[i], pop_pbasal[i])
        )
        # Extract terminal phase-space coordinates
        final_atp[i] = solution[-1, 0]
        final_glu[i] = solution[-1, 1]

    # ---------------------------------------------------------
    # 5. DATA COMPILATION & BIFURCATION ANALYSIS
    # ---------------------------------------------------------
    print("Integration complete. Compiling dataset...")
    df_results = pd.DataFrame({
        'Ancestry_Cohort': pop_labels,
        'Parameterized_PaCO2': pop_paco2,
        'Parameterized_R_stress': pop_rstress,
        'Terminal_ATP': final_atp,
        'Terminal_Glutamate': final_glu
    })
    
    # Mathematical criteria for simulated Saddle-Node Bifurcation
    # Defined exclusively by the topological collapse of the bioenergetic state variable (ATP < 0.5 mM)
    df_results['Simulated_Bifurcation'] = (df_results['Terminal_ATP'] < 0.5)
    
    # Output Statistical Summary
    print("\n--- INCIDENCE OF SIMULATED BIFURCATION BY COHORT ---")
    bifurcation_summary = df_results.groupby('Ancestry_Cohort')['Simulated_Bifurcation'].mean() * 100
    print(bifurcation_summary.to_string(float_format="%.2f%%"))
    
    # Export to CSV
    export_filename = 'ISB_PanAncestry_Simulation_Results.csv'
    df_results.to_csv(export_filename, index=False)
    print(f"\nPhase-space coordinates exported successfully to '{export_filename}'.")


