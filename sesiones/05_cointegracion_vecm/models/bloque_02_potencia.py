import matplotlib.pyplot as plt
import numpy as np
import statsmodels.api as sm
from statsmodels.tsa.stattools import adfuller
from statsmodels.tsa.vector_ar.vecm import VECM

from core.config import RNG
from core.diagnostics import estimar_rango_johansen, johansen_report
from core.plotting import guardar_figura


def run() -> dict:
    print("\n", "=" * 80)
    print("BLOQUE 2: SISTEMAS DE POTENCIA - ESTABILIDAD DE FRECUENCIA")
    print("=" * 80, "\n")

    print("\n--- EJERCICIO 3: Frecuencia y Potencia de Intercambio ---\n")

    M, D, f0 = 10, 0.5, 50
    T3, dt = 1000, 0.1

    P_m = np.ones(T3)
    P_e = np.ones(T3)
    P_e[499:] = 0.85

    f = np.zeros(T3)
    f[0] = f0
    for t in range(1, T3):
        df = (P_m[t] - P_e[t] - D * (f[t - 1] - f0)) / M
        f[t] = f[t - 1] + df * dt + RNG.normal(0, 0.02)

    P_link = 0.8 * (f - f0) + np.cumsum(RNG.normal(0, 0.01, T3)) + 0.5
    f = f + np.cumsum(RNG.normal(0, 0.001, T3))
    P_link = P_link + np.cumsum(RNG.normal(0, 0.001, T3))

    tiempo3 = np.arange(1, T3 + 1)

    fig, ax1 = plt.subplots()
    ax1.plot(tiempo3, f, label="Frecuencia (Hz)", linewidth=0.8)
    ax1.axvline(500, linestyle="--", color="darkred")
    ax1.set_xlabel("Tiempo (muestras)")
    ax1.set_ylabel("Frecuencia (Hz)")
    ax2 = ax1.twinx()
    ax2.plot(tiempo3, P_link, label="Potencia", color="orange", linewidth=0.8)
    ax2.set_ylabel("Potencia (p.u.)")
    ax1.set_title("Frecuencia y Potencia de Intercambio - Sistema Interconectado")
    fig.legend(loc="lower center", ncol=2)
    guardar_figura(fig, "ej03_frecuencia_potencia.png")

    datos_johansen3 = np.column_stack([f, P_link])
    johansen3 = johansen_report(datos_johansen3, ["Frecuencia", "Potencia"])

    beta_hat = johansen3.evec[:, 0].copy()
    beta_hat /= beta_hat[0]
    print("\nVector de cointegración estimado (normalizado):")
    print(beta_hat)

    ECT_ej3 = f - beta_hat[1] * P_link
    evento = np.where(tiempo3 < 500, "Normal", "Deslastre")

    fig, ax = plt.subplots()
    for est, color in [("Normal", "blue"), ("Deslastre", "red")]:
        mask = evento == est
        ax.plot(tiempo3[mask], ECT_ej3[mask], color=color, label=est, linewidth=0.8)
    ax.axhline(0, linestyle="--", color="gray")
    ax.axvline(500, linestyle="--", color="darkred", linewidth=0.8)
    ax.set_title(f"ECT del sistema frecuencia-potencia (β₂={-beta_hat[1]:.3f})")
    ax.set_xlabel("Tiempo (muestras)")
    ax.set_ylabel("ECT = f - β*P_link")
    ax.legend(loc="lower center")
    guardar_figura(fig, "ej03_ect_frecuencia.png")

    r_ej3 = max(1, estimar_rango_johansen(johansen3))
    vecm3 = VECM(datos_johansen3, k_ar_diff=1, coint_rank=r_ej3, deterministic="ci")
    vecm3_res = vecm3.fit()
    alpha_hat = vecm3_res.alpha[:, 0]
    print("\nVelocidades de ajuste (α):")
    print(alpha_hat)

    fig, ax = plt.subplots()
    colores = ["darkgreen" if a > 0 else "darkred" for a in alpha_hat]
    ax.bar(["Frecuencia", "Potencia"], alpha_hat, color=colores, width=0.6)
    ax.axhline(0, linestyle="--", color="gray")
    ax.set_title("Velocidades de ajuste (α) del VECM")
    ax.set_ylabel("Coeficiente de ajuste α")
    guardar_figura(fig, "ej03_alpha_vecm.png")

    print("\n--- EJERCICIO 4: Detección de Isla vía Cointegración de Fases ---\n")

    T4, fs = 1000, 60
    theta_red = np.cumsum(RNG.normal(0.01, 0.001, T4)) + 2 * np.pi * 50 / fs * np.arange(1, T4 + 1)
    theta_inv = theta_red + RNG.normal(0, 0.02, T4)

    f_isla = 50.5
    idx_isla = np.arange(600, T4 + 1)
    theta_inv[599:] += 2 * np.pi * (f_isla - 50) / fs * np.arange(1, len(idx_isla) + 1)

    delta_theta = theta_inv - theta_red
    fase4 = np.where(np.arange(1, T4 + 1) < 600, "Normal", "Isla")

    fig, ax = plt.subplots()
    for est, color in [("Normal", "blue"), ("Isla", "red")]:
        mask = fase4 == est
        ax.plot(np.arange(1, T4 + 1)[mask], delta_theta[mask], color=color, label=est, linewidth=0.8)
    ax.axvline(600, linestyle="--", color="darkred", linewidth=0.8)
    ax.set_title("Diferencia de fases entre inversor y red")
    ax.set_xlabel("Tiempo (muestras)")
    ax.set_ylabel("Δθ = θ_inv - θ_red (rad)")
    ax.legend(loc="lower center")
    guardar_figura(fig, "ej04_diferencia_fases.png")

    ventana_isla = 120
    estadisticos_tau = np.full(T4 - ventana_isla + 1, np.nan)
    for i in range(T4 - ventana_isla + 1):
        idx = slice(i, i + ventana_isla)
        t_win = np.arange(1, ventana_isla + 1)
        modelo = sm.OLS(delta_theta[idx], sm.add_constant(t_win)).fit()
        residuales = modelo.resid
        estadisticos_tau[i] = adfuller(residuales, maxlag=2, regression="n", autolag=None)[0]

    tiempo_ventana = np.arange(1, T4 - ventana_isla + 2) + ventana_isla / 2
    evento4 = np.where(tiempo_ventana < 600, "Normal", "Isla")

    fig, ax = plt.subplots()
    for est, color in [("Normal", "blue"), ("Isla", "red")]:
        mask = evento4 == est
        ax.plot(tiempo_ventana[mask], estadisticos_tau[mask], color=color, label=est, linewidth=0.8)
    ax.axhline(-3.5, linestyle="--", color="red", linewidth=1)
    ax.axhline(-2.9, linestyle="--", color="orange", linewidth=1)
    ax.axvline(600, linestyle="--", color="darkred", linewidth=0.8)
    ax.annotate("Valor crítico 5%", xy=(200, -3.2), color="red")
    ax.annotate("Valor crítico 10%", xy=(200, -2.6), color="orange")
    ax.set_title(f"Estadístico τ de Engle-Granger (ventana={ventana_isla} muestras)")
    ax.set_xlabel("Tiempo (muestras)")
    ax.set_ylabel("Estadístico τ")
    ax.legend(loc="lower center")
    guardar_figura(fig, "ej04_engle_granger_ventanas.png")

    rocof = np.diff(delta_theta) * fs
    rocof_abs = np.abs(rocof)
    evento_rocof = np.where(np.arange(2, T4 + 1) < 600, "Normal", "Isla")

    fig, ax = plt.subplots()
    for est, color in [("Normal", "blue"), ("Isla", "red")]:
        t_r = np.arange(2, T4 + 1)
        mask = evento_rocof == est
        ax.plot(t_r[mask], rocof_abs[mask], color=color, label=est, linewidth=0.5)
    ax.axhline(0.2, linestyle="--", color="red")
    ax.axvline(600, linestyle="--", color="darkred", linewidth=0.8)
    ax.set_title("ROCOF (Rate of Change of Frequency) - Método tradicional")
    ax.set_xlabel("Tiempo (muestras)")
    ax.set_ylabel("|dΔθ/dt| (Hz)")
    ax.legend(loc="upper right")
    guardar_figura(fig, "ej04_rocof.png")

    return {
        "johansen3": johansen3,
        "beta_hat": beta_hat,
        "alpha_hat": alpha_hat,
    }
