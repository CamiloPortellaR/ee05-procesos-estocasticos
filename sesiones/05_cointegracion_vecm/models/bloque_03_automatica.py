import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm

from core.config import RNG
from core.diagnostics import johansen_report, vectores_cointegracion
from core.plotting import guardar_figura


def run() -> dict:
    print("\n", "=" * 80)
    print("BLOQUE 3: AUTOMÁTICA AVANZADA - MIMO Y RETARDOS")
    print("=" * 80, "\n")

    print("\n--- EJERCICIO 5: Desacoplamiento de perturbaciones térmicas ---\n")

    T5 = 500
    Ta = np.cumsum(RNG.normal(0, 0.02, T5)) + 20
    G = np.array([[1.2, 0.1, 0.0], [0.1, 1.1, 0.2], [0.0, 0.2, 0.9]])

    P1 = np.cumsum(RNG.normal(0.01, 0.02, T5)) + 5
    P2 = np.cumsum(RNG.normal(0.01, 0.02, T5)) + 5
    P3 = np.cumsum(RNG.normal(0.01, 0.02, T5)) + 5

    T1 = G[0, 0] * P1 + G[0, 1] * P2 + G[0, 2] * P3 + 0.3 * Ta + RNG.normal(0, 0.1, T5)
    T2 = G[1, 0] * P1 + G[1, 1] * P2 + G[1, 2] * P3 + 0.2 * Ta + RNG.normal(0, 0.1, T5)
    T3 = G[2, 0] * P1 + G[2, 1] * P2 + G[2, 2] * P3 + 0.1 * Ta + RNG.normal(0, 0.1, T5)

    tiempo5 = np.arange(1, T5 + 1)

    fig, ax = plt.subplots()
    ax.plot(tiempo5, T1, label="T1", linewidth=0.8)
    ax.plot(tiempo5, T2, label="T2", linewidth=0.8)
    ax.plot(tiempo5, T3, label="T3", linewidth=0.8)
    ax.plot(tiempo5, Ta, label="Ta (amb)", linewidth=0.6, linestyle="--")
    ax.set_title("Temperaturas del horno de 3 zonas")
    ax.set_xlabel("Tiempo (muestras)")
    ax.set_ylabel("Temperatura (°C)")
    ax.legend(loc="lower center")
    guardar_figura(fig, "ej05_temperaturas_horno.png")

    datos_vect = np.column_stack([T1, T2, T3, P1, P2, P3])
    johansen5 = johansen_report(datos_vect, ["T1", "T2", "T3", "P1", "P2", "P3"])

    if johansen5.lr1[2] > johansen5.cvt[2, 1]:
        r_estimado = 3
        print("\nRango de cointegración estimado: r = 3")
    else:
        r_estimado = 2
        print("\nRango de cointegración estimado: r = 2")

    V = vectores_cointegracion(johansen5, r_estimado, normalizar_idx=0)
    print("\nVectores de cointegración (normalizados por T1):")
    print(V)

    ECTs = datos_vect @ V
    nombres_ect = [f"ECT{i + 1}" for i in range(r_estimado)]

    fig, axes = plt.subplots(1, r_estimado, figsize=(5 * r_estimado, 4), squeeze=False)
    for i, ax in enumerate(axes[0]):
        ax.plot(tiempo5, ECTs[:, i], linewidth=0.8)
        ax.axhline(0, linestyle="--", color="gray")
        ax.set_title(nombres_ect[i])
        ax.set_xlabel("Tiempo (muestras)")
    fig.suptitle(f"ECTs del sistema térmico (r={r_estimado})")
    plt.tight_layout()
    guardar_figura(fig, "ej05_ects_termico.png")

    print("\n--- EJERCICIO 6: Degradación de válvulas con histéresis ---\n")

    T6 = 800
    A = np.cumsum(RNG.normal(0, 0.02, T6)) + 10
    delta = 0.001 * np.arange(1, T6 + 1) + np.cumsum(RNG.normal(0, 0.01, T6))
    Q = A + delta + RNG.normal(0, 0.05, T6)
    datos_ej6 = pd.DataFrame({"tiempo": np.arange(1, T6 + 1), "A": A, "Q": Q, "delta": delta})

    fig, ax = plt.subplots()
    ax.plot(datos_ej6["tiempo"], A, label="Apertura A", linewidth=0.8)
    ax.plot(datos_ej6["tiempo"], Q, label="Caudal Q", linewidth=0.8)
    ax.set_title("Degradación de válvula - Apertura vs Caudal")
    ax.set_xlabel("Tiempo (muestras)")
    ax.set_ylabel("Valor")
    ax.legend(loc="lower center")
    guardar_figura(fig, "ej06_apertura_caudal.png")

    fig, ax = plt.subplots()
    ax.scatter(A, Q, alpha=0.3, color="blue", s=10)
    z = np.polyfit(A, Q, 1)
    ax.plot(np.sort(A), np.poly1d(z)(np.sort(A)), color="red", linewidth=1, label="Regresión lineal")
    ax.set_title("Relación Caudal vs Apertura")
    ax.set_xlabel("Apertura A")
    ax.set_ylabel("Caudal Q")
    ax.legend()
    guardar_figura(fig, "ej06_q_vs_a.png")

    ventana_hi = 50
    beta_ventana = np.full(T6 - ventana_hi + 1, np.nan)
    HI = np.full(T6 - ventana_hi + 1, np.nan)
    for i in range(T6 - ventana_hi + 1):
        idx = slice(i, i + ventana_hi)
        beta = sm.OLS(Q[idx], sm.add_constant(A[idx])).fit().params[1]
        beta_ventana[i] = beta
        HI[i] = 1 / (1 + abs(beta - beta_ventana[0]))

    tiempo_hi = np.arange(1, T6 - ventana_hi + 2) + ventana_hi / 2
    fig, ax = plt.subplots()
    ax.plot(tiempo_hi, HI, color="darkblue", linewidth=1)
    ax.axhline(0.8, linestyle="--", color="orange", linewidth=1)
    ax.axhline(0.6, linestyle="--", color="red", linewidth=1)
    ax.annotate("Alerta temprana", xy=(100, 0.85), color="orange")
    ax.annotate("Mantenimiento urgente", xy=(100, 0.65), color="red")
    ax.set_title("Health Index (HI) de la válvula")
    ax.set_xlabel("Tiempo (muestras)")
    ax.set_ylabel("Health Index")
    guardar_figura(fig, "ej06_health_index.png")

    return {"r_estimado": r_estimado}
