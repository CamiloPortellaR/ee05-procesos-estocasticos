import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.diagnostic import breaks_cusumolsresid

from core.config import RNG
from core.diagnostics import adf_report, johansen_report
from core.plotting import guardar_figura


def run() -> None:
    print("\n", "=" * 80)
    print("BLOQUE 1: FUNDAMENTOS - DIAGNÓSTICO DE INSTRUMENTACIÓN")
    print("=" * 80, "\n")

    print("\n--- EJERCICIO 1: Señales sintéticas con deriva ---\n")

    T = 500
    mu = np.cumsum(RNG.normal(0, 0.1, T))
    epsilon1 = RNG.normal(0, 0.05, T)
    epsilon2 = RNG.normal(0, 0.05, T)
    tiempo = np.arange(1, T + 1)

    I1 = mu + epsilon1
    I2 = mu + 0.05 * tiempo + epsilon2

    fig, ax = plt.subplots()
    ax.plot(tiempo, I1, label="Sensor 1", linewidth=1)
    ax.plot(tiempo, I2, label="Sensor 2", linewidth=1)
    ax.set_title("Evolución de las corrientes de los sensores")
    ax.set_xlabel("Tiempo (muestras)")
    ax.set_ylabel("Corriente (A)")
    ax.legend(loc="lower center")
    guardar_figura(fig, "ej01_series_sensores.png")

    adf_report(I1, "I1")
    adf_report(I2, "I2")

    datos_johansen = np.column_stack([I1, I2])
    johansen_report(datos_johansen, ["I1", "I2"])

    ECT_ej1 = I1 - I2
    sigma_inicial = np.std(ECT_ej1[:200])

    fig, ax = plt.subplots()
    ax.plot(tiempo, ECT_ej1, color="darkred", linewidth=1)
    ax.axhline(0, linestyle="--", color="gray")
    ax.axhline(2 * sigma_inicial, linestyle=":", color="blue")
    ax.axhline(-2 * sigma_inicial, linestyle=":", color="blue")
    ax.set_title("Término de Corrección del Error (ECT) - Ejercicio 1")
    ax.set_xlabel("Tiempo (muestras)")
    ax.set_ylabel("ECT")
    ax.annotate("Deriva del sensor 2\ndetectada", xy=(450, np.max(ECT_ej1) * 0.8), color="darkred")
    guardar_figura(fig, "ej01_ect_sensores.png")

    print("\n--- EJERCICIO 2: Modelo VECM en lazo de control ---\n")

    T2 = 500
    K1, K2 = 2.0, 2.5
    punto_cambio = 300

    u = np.cumsum(RNG.normal(0, 0.05, T2)) + 50
    y = np.zeros(T2)
    for t in range(T2):
        k_actual = K1 if t < punto_cambio else K2
        y[t] = k_actual * u[t] + RNG.normal(0, 0.1)

    datos_ej2 = pd.DataFrame({"tiempo": np.arange(1, T2 + 1), "u": u, "y": y})

    fig, ax = plt.subplots()
    ax.plot(datos_ej2["tiempo"], u, label="Consigna u(t)", linewidth=0.8)
    ax.plot(datos_ej2["tiempo"], y, label="Salida y(t)", linewidth=0.8)
    ax.axvline(punto_cambio, linestyle="--", color="darkred", linewidth=1)
    ax.annotate("Cambio de ganancia", xy=(punto_cambio + 20, np.max(y)), color="darkred")
    ax.set_title("Lazo de control hidráulico - Consigna vs Salida")
    ax.set_xlabel("Tiempo (muestras)")
    ax.set_ylabel("Valor")
    ax.legend(loc="lower center")
    guardar_figura(fig, "ej02_consigna_salida.png")

    K_nominal = 2
    ECT_ej2 = y - K_nominal * u
    fase = np.where(np.arange(1, T2 + 1) < punto_cambio, "Normal", "Falla")

    fig, ax = plt.subplots()
    for estado, color in [("Normal", "blue"), ("Falla", "red")]:
        mask = fase == estado
        ax.plot(np.arange(1, T2 + 1)[mask], ECT_ej2[mask], color=color, label=estado, linewidth=1)
    ax.axhline(0, linestyle="--", color="gray")
    ax.axvline(punto_cambio, linestyle="--", color="darkred", linewidth=0.8)
    ax.set_title("ECT en lazo de control - Detección de cambio de ganancia")
    ax.set_xlabel("Tiempo (muestras)")
    ax.set_ylabel("ECT = y - 2u")
    ax.legend(loc="lower center")
    guardar_figura(fig, "ej02_ect_lazo_control.png")

    ventana = 20
    ect_series = pd.Series(ECT_ej2)
    media_movil = ect_series.rolling(ventana, min_periods=ventana).mean()
    sd_movil = ect_series.rolling(ventana, min_periods=ventana).std()
    umbral_superior = media_movil + 2 * sd_movil
    umbral_inferior = media_movil - 2 * sd_movil
    alarma = ((ECT_ej2 > umbral_superior) | (ECT_ej2 < umbral_inferior)) & umbral_superior.notna()

    fig, ax = plt.subplots()
    ax.plot(np.arange(1, T2 + 1), ECT_ej2, color="gray", linewidth=0.5, label="ECT")
    ax.plot(np.arange(1, T2 + 1), media_movil, color="blue", linewidth=1, label="Media móvil")
    ax.fill_between(
        np.arange(1, T2 + 1),
        umbral_inferior,
        umbral_superior,
        color="blue",
        alpha=0.15,
        label="±2σ",
    )
    ax.scatter(np.arange(1, T2 + 1)[alarma], ECT_ej2[alarma], color="red", s=20, zorder=5, label="Alarma")
    ax.axvline(punto_cambio, linestyle="--", color="darkred", linewidth=0.8)
    ax.set_title("Detección de cambio de ganancia vía ECT + Alarma")
    ax.set_xlabel("Tiempo (muestras)")
    ax.set_ylabel("ECT")
    ax.legend(loc="upper left")
    guardar_figura(fig, "ej02_alarma_ect.png")

    res_cusum = sm.OLS(ECT_ej2, sm.add_constant(np.ones(T2))).fit()
    sup_b, pval_cusum, crit_cusum = breaks_cusumolsresid(res_cusum.resid, ddof=1)
    cusum_stat = np.cumsum(res_cusum.resid)
    escala = res_cusum.resid.std(ddof=1) * np.sqrt(np.arange(1, T2 + 1))
    cusum_stat = cusum_stat / escala
    print(f"CUSUM sup-B: {sup_b:.4f} | p-valor: {pval_cusum:.4f}")
    fig, ax = plt.subplots()
    ax.plot(np.arange(len(cusum_stat)), cusum_stat, color="darkblue")
    ax.axhline(0, linestyle="--", color="gray")
    ax.set_title("Test CUSUM para el ECT")
    ax.set_xlabel("Observación")
    ax.set_ylabel("Estadístico CUSUM")
    guardar_figura(fig, "ej02_cusum_ect.png")
