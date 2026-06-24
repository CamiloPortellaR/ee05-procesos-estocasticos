import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm

from core.config import RNG
from core.diagnostics import johansen_report, vectores_cointegracion
from core.plotting import guardar_figura


def run() -> dict:
    print("\n", "=" * 80)
    print("BLOQUE 5: PROYECTO INTEGRADOR - SENSORES MÚLTIPLES Y CPCo")
    print("=" * 80, "\n")

    print("\n--- EJERCICIO 9: Cointegración en sensores redundantes (4 giróscopos) ---\n")

    T9 = 500
    omega = np.cumsum(RNG.normal(0, 0.01, T9))
    sesgo1 = np.cumsum(RNG.normal(0, 0.001, T9))
    sesgo2 = 0.1 + np.cumsum(RNG.normal(0, 0.001, T9))
    sesgo3 = -0.05 + np.cumsum(RNG.normal(0, 0.001, T9))
    sesgo4 = 0.15 + np.cumsum(RNG.normal(0, 0.001, T9))

    n_fallo = T9 - 349
    sesgo2[349:] += 0.0005 * np.arange(1, n_fallo + 1) ** 2

    G1 = omega + sesgo1 + RNG.normal(0, 0.02, T9)
    G2 = omega + sesgo2 + RNG.normal(0, 0.02, T9)
    G3 = omega + sesgo3 + RNG.normal(0, 0.02, T9)
    G4 = omega + sesgo4 + RNG.normal(0, 0.02, T9)
    tiempo9 = np.arange(1, T9 + 1)

    fig, ax = plt.subplots()
    ax.plot(tiempo9, G1, label="G1", linewidth=0.7)
    ax.plot(tiempo9, G2, label="G2", linewidth=0.7)
    ax.plot(tiempo9, G3, label="G3", linewidth=0.7)
    ax.plot(tiempo9, G4, label="G4", linewidth=0.7)
    ax.axvline(350, linestyle="--", color="darkred", linewidth=1)
    ax.set_title("Lecturas de 4 giróscopos redundantes")
    ax.set_xlabel("Tiempo (muestras)")
    ax.set_ylabel("Velocidad angular (rad/s)")
    ax.legend(loc="lower center")
    guardar_figura(fig, "ej09_giroscopos.png")

    datos_giro = np.column_stack([G1, G2, G3, G4])
    johansen9 = johansen_report(datos_giro, ["G1", "G2", "G3", "G4"])

    r_giro = 3
    V9 = vectores_cointegracion(johansen9, r_giro, normalizar_idx=0)
    print("\nVectores de cointegración (normalizados por G1):")
    print(V9)

    ECTs_giro = datos_giro @ V9
    fig, axes = plt.subplots(1, r_giro, figsize=(5 * r_giro, 4), squeeze=False)
    for i, ax in enumerate(axes[0]):
        ax.plot(tiempo9, ECTs_giro[:, i], linewidth=0.8)
        ax.axhline(0, linestyle="--", color="gray")
        ax.axvline(350, linestyle="--", color="darkred", linewidth=0.8)
        ax.set_title(f"ECT_G{i + 1}")
        ax.set_xlabel("Tiempo (muestras)")
    fig.suptitle("ECTs de los 4 giróscopos")
    plt.tight_layout()
    guardar_figura(fig, "ej09_ects_giroscopos.png")

    print("\n--- EJERCICIO 10: Control Predictivo basado en Cointegración (CPCo) ---\n")

    T10, dt10 = 300, 0.1
    h1 = np.zeros(T10)
    h2 = np.zeros(T10)
    q1 = np.zeros(T10)
    q2 = np.zeros(T10)
    h1[0], h2[0], q1[0], q2[0] = 1.0, 0.8, 0.5, 0.4

    A1, A2, k12 = 1.0, 1.2, 0.2
    fuga = np.zeros(T10)
    fuga[149:] = 0.02 * np.arange(1, T10 - 149 + 1)

    for t in range(1, T10):
        dh1 = (q1[t - 1] - k12 * (h1[t - 1] - h2[t - 1])) / A1 * dt10
        dh2 = (q2[t - 1] + k12 * (h1[t - 1] - h2[t - 1]) - fuga[t]) / A2 * dt10
        h1[t] = h1[t - 1] + dh1 + RNG.normal(0, 0.005)
        h2[t] = h2[t - 1] + dh2 + RNG.normal(0, 0.005)
        q1[t] = q1[t - 1] + RNG.normal(0, 0.01)
        q2[t] = q2[t - 1] + RNG.normal(0, 0.01)

    tiempo10 = np.arange(1, T10 + 1)
    fig, ax = plt.subplots()
    ax.plot(tiempo10, h1, label="h1", linewidth=1)
    ax.plot(tiempo10, h2, label="h2", linewidth=1)
    ax.axvline(150, linestyle="--", color="darkred", linewidth=0.8)
    ax.set_title("Sistema de 2 tanques acoplados")
    ax.set_xlabel("Tiempo (muestras)")
    ax.set_ylabel("Nivel (m)")
    ax.legend(loc="lower center")
    guardar_figura(fig, "ej10_tanques.png")

    modelo_ect = sm.OLS(h1[:100], sm.add_constant(h2[:100])).fit()
    beta_cpco = modelo_ect.params[1]
    print(f"\nRelación de cointegración estimada: h1 = {beta_cpco:.3f} * h2 + cte")

    ECT_cpco = h1 - beta_cpco * h2
    fig, ax = plt.subplots()
    ax.plot(tiempo10, ECT_cpco, color="darkblue", linewidth=1)
    ax.axhline(0, linestyle="--", color="gray")
    ax.axvline(150, linestyle="--", color="darkred", linewidth=0.8)
    ax.set_title(f"ECT para CPCo (β={beta_cpco:.3f})")
    ax.set_xlabel("Tiempo (muestras)")
    ax.set_ylabel("ECT = h1 - β*h2")
    guardar_figura(fig, "ej10_ect_cpco.png")

    h1_cpco = h1.copy()
    h2_cpco = h2.copy()
    h1_cpco[149:] -= ECT_cpco[149:] * 0.5
    h2_cpco[149:] += ECT_cpco[149:] * 0.3

    fig, ax = plt.subplots()
    ax.plot(tiempo10, h1, label="h1 (sin CPCo)", linewidth=0.8, color="red")
    ax.plot(tiempo10, h1_cpco, label="h1 (con CPCo)", linewidth=0.8, color="blue")
    ax.plot(tiempo10, h2, label="h2", linewidth=0.8, color="darkgreen")
    ax.axvline(150, linestyle="--", color="darkred", linewidth=0.8)
    ax.set_title("Efecto del CPCo en el sistema de tanques")
    ax.set_xlabel("Tiempo (muestras)")
    ax.set_ylabel("Nivel (m)")
    ax.legend(loc="lower center")
    guardar_figura(fig, "ej10_cpco_comparativa.png")

    return {"johansen9": johansen9, "r_giro": r_giro}
