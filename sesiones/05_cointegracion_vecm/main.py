# %% Configuración de rutas del mini-proyecto
import sys
from pathlib import Path

import numpy as np
import pandas as pd

RUTA_SESION = Path(__file__).resolve().parent
if str(RUTA_SESION) not in sys.path:
    sys.path.insert(0, str(RUTA_SESION))

from core.config import RUTA_GRAFICOS
from core.diagnostics import estimar_rango_johansen
from models.bloque_01_fundamentos import run as run_bloque_01
from models.bloque_02_potencia import run as run_bloque_02
from models.bloque_03_automatica import run as run_bloque_03
from models.bloque_04_mercados import run as run_bloque_04
from models.bloque_05_integrador import run as run_bloque_05

# %% Bloque 1 — Fundamentos (Ejercicios 1-2)
run_bloque_01()

# %% Bloque 2 — Sistemas de potencia (Ejercicios 3-4)
ctx_potencia = run_bloque_02()

# %% Bloque 3 — Automática avanzada (Ejercicios 5-6)
ctx_automatica = run_bloque_03()

# %% Bloque 4 — Mercados eléctricos (Ejercicios 7-8)
ctx_mercados = run_bloque_04()

# %% Bloque 5 — Proyecto integrador (Ejercicios 9-10)
ctx_integrador = run_bloque_05()

# %% Resumen final de resultados
print("\n", "=" * 80)
print("RESUMEN DE RESULTADOS Y ESTADÍSTICAS CLAVE")
print("=" * 80, "\n")

johansen3 = ctx_potencia["johansen3"]
beta_hat = ctx_potencia["beta_hat"]
alpha_hat = ctx_potencia["alpha_hat"]
r_estimado = ctx_automatica["r_estimado"]
johansen7 = ctx_mercados["johansen7"]
rmse_vecm = ctx_mercados["rmse_vecm"]
rmse_arima = ctx_mercados["rmse_arima"]
r_giro = ctx_integrador["r_giro"]

r_ej7 = estimar_rango_johansen(johansen7)
resultados = pd.DataFrame(
    {
        "Ejercicio": [1, 3, 4, 5, 7, 9],
        "Descripcion": [
            "Sensores con deriva",
            "Frecuencia-Potencia",
            "Detección de isla",
            "Horno térmico",
            "Precios nodales",
            "Giróscopos",
        ],
        "Variables": [2, 2, 2, 6, 3, 4],
        "Rango_estimado": [1, estimar_rango_johansen(johansen3), np.nan, r_estimado, r_ej7, r_giro],
        "Interpretacion": [
            "ECT detecta deriva del sensor 2",
            f"β = {beta_hat[1]:.3f}",
            "τ cae por debajo de -3.5 en isla",
            f"r = {r_estimado} relaciones",
            "Congestión reduce rango de 2 a 1",
            "Fallo en G2 detectado por ECTs",
        ],
    }
)
print(resultados.to_string(index=False))

print("\n--- Velocidades de ajuste (Ejercicio 3) ---")
print(f"α_frecuencia: {alpha_hat[0]:.4f}")
print(f"α_potencia: {alpha_hat[1]:.4f}")
print(
    "Interpretación: La frecuencia se ajusta",
    "más rápido" if abs(alpha_hat[0]) > abs(alpha_hat[1]) else "más lento",
    "que la potencia al desequilibrio",
)

print("\n--- Mejora de predicción (Ejercicio 8) ---")
print(f"RMSE VECM: {rmse_vecm:.4f}")
print(f"RMSE ARIMA: {rmse_arima:.4f}")
print(f"Mejora: {(rmse_arima - rmse_vecm) / rmse_arima * 100:.2f}%")

print("\n", "=" * 80)
print("FIN DEL SCRIPT - Todas las gráficas generadas correctamente")
print(f"Gráficos guardados en: {RUTA_GRAFICOS}")
print("=" * 80, "\n")
