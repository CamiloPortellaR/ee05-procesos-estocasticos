import matplotlib.pyplot as plt
import numpy as np
from statsmodels.tsa.vector_ar.vecm import VECM

from core.config import RNG
from core.diagnostics import ajuste_arima_aic, estimar_rango_johansen, johansen_report, vectores_cointegracion
from core.plotting import guardar_figura


def run() -> dict:
    print("\n", "=" * 80)
    print("BLOQUE 4: MERCADOS ELÉCTRICOS - ARBITRAJE Y CONGESTIÓN")
    print("=" * 80, "\n")

    print("\n--- EJERCICIO 7: Cointegración de precios nodales ---\n")

    T7 = 1000
    gas = np.cumsum(RNG.normal(0.01, 0.02, T7)) + 30
    PA = gas + RNG.normal(0, 0.5, T7) + 5
    PB = gas + RNG.normal(0, 0.5, T7) + 3
    PC = gas + RNG.normal(0, 0.5, T7) + 2

    tiempo7 = np.arange(1, T7 + 1)
    congestion = np.where((tiempo7 >= 400) & (tiempo7 <= 600), 2 * np.sin((tiempo7 - 400) / 20), 0)
    PA[399:600] += congestion[399:600]
    PB[399:600] -= congestion[399:600]

    fig, ax = plt.subplots()
    ax.plot(tiempo7, PA, label="Nodo A", linewidth=0.8)
    ax.plot(tiempo7, PB, label="Nodo B", linewidth=0.8)
    ax.plot(tiempo7, PC, label="Nodo C", linewidth=0.8)
    ax.plot(tiempo7, gas, label="Gas", linewidth=0.6, linestyle="--")
    ax.axvline(400, linestyle="--", color="darkred", alpha=0.7)
    ax.axvline(600, linestyle="--", color="darkred", alpha=0.7)
    ax.set_title("Precios spot en 3 nodos del mercado eléctrico")
    ax.set_xlabel("Tiempo (horas)")
    ax.set_ylabel("Precio ($/MWh)")
    ax.legend(loc="lower center")
    guardar_figura(fig, "ej07_precios_nodales.png")

    datos_precios = np.column_stack([PA, PB, PC])
    johansen7 = johansen_report(datos_precios, ["PA", "PB", "PC"])

    V7 = vectores_cointegracion(johansen7, 2, normalizar_idx=0)
    print("\nVectores de cointegración (normalizados):")
    print(V7)

    ECT_CA = PC - V7[2, 0] * PA
    media_ect = np.mean(ECT_CA[:300])
    sd_ect = np.std(ECT_CA[:300])
    umbral_sup = media_ect + 2 * sd_ect
    umbral_inf = media_ect - 2 * sd_ect

    senal = np.full(T7, "Neutral", dtype=object)
    senal[ECT_CA > umbral_sup] = "Vender C / Comprar A"
    senal[ECT_CA < umbral_inf] = "Comprar C / Vender A"

    fig, ax = plt.subplots()
    for s, color in [
        ("Neutral", "gray"),
        ("Vender C / Comprar A", "red"),
        ("Comprar C / Vender A", "green"),
    ]:
        mask = senal == s
        ax.plot(tiempo7[mask], ECT_CA[mask], color=color, label=s, linewidth=0.8, marker=".")
    ax.axhline(umbral_sup, linestyle="--", color="red", linewidth=1)
    ax.axhline(umbral_inf, linestyle="--", color="green", linewidth=1)
    ax.axvline(400, linestyle="--", color="darkred", alpha=0.5)
    ax.axvline(600, linestyle="--", color="darkred", alpha=0.5)
    ax.set_title("ECT para trading de pares (Nodo C vs Nodo A)")
    ax.set_xlabel("Tiempo (horas)")
    ax.set_ylabel("ECT = PC - β*PA")
    ax.legend(loc="upper right", fontsize=8)
    guardar_figura(fig, "ej07_ect_trading.png")

    print("\n--- EJERCICIO 8: Predicción de precios VECM vs ARIMA ---\n")

    train_idx = slice(0, 700)
    test_idx = slice(700, T7)
    train_data = datos_precios[train_idx]
    test_len = T7 - 700

    vecm_model = VECM(train_data, k_ar_diff=1, coint_rank=2, deterministic="ci")
    vecm_res = vecm_model.fit()
    print(vecm_res.summary())

    pred_vecm = vecm_res.predict(steps=test_len)
    pred_PA_vecm = pred_vecm[:, 0]

    arima_res, arima_ord = ajuste_arima_aic(PA[:700], d=1)
    pred_arima = arima_res.forecast(steps=test_len)

    rmse_vecm = np.sqrt(np.mean((PA[700:] - pred_PA_vecm) ** 2))
    rmse_arima = np.sqrt(np.mean((PA[700:] - pred_arima) ** 2))

    print("\n--- Comparación de predicciones ---")
    print(f"RMSE VECM: {rmse_vecm:.4f}")
    print(f"RMSE ARIMA: {rmse_arima:.4f}")
    print(f"Mejora VECM vs ARIMA: {(rmse_arima - rmse_vecm) / rmse_arima * 100:.2f}%")

    fig, ax = plt.subplots()
    t_test = np.arange(701, T7 + 1)
    ax.plot(t_test, PA[700:], label="Real", linewidth=1)
    ax.plot(t_test, pred_PA_vecm, label="VECM", linewidth=0.8, linestyle="--")
    ax.plot(t_test, pred_arima, label="ARIMA", linewidth=0.8, linestyle=":")
    ax.set_title(f"Predicción de precios — RMSE VECM: {rmse_vecm:.3f} | ARIMA: {rmse_arima:.3f}")
    ax.set_xlabel("Tiempo (horas)")
    ax.set_ylabel("Precio nodo A ($/MWh)")
    ax.legend(loc="lower center")
    guardar_figura(fig, "ej08_prediccion_vecm_arima.png")

    return {
        "johansen7": johansen7,
        "rmse_vecm": rmse_vecm,
        "rmse_arima": rmse_arima,
    }
