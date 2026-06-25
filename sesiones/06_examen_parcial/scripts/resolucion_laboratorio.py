# %% Configuración
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statsmodels.tsa.api import VAR

RUTA_SESION = Path(__file__).resolve().parent.parent
RUTA_DATA = RUTA_SESION / "data" / "base_datos_var.xlsx"
RUTA_OUTPUTS = RUTA_SESION / "outputs"
RUTA_OUTPUTS.mkdir(parents=True, exist_ok=True)

VARIABLES = ["frecuencia_hz", "potencia_activa_mw", "potencia_reactiva_mvar"]
ETIQUETAS = {
    "frecuencia_hz": "Frecuencia (Hz)",
    "potencia_activa_mw": "Potencia Activa (MW)",
    "potencia_reactiva_mvar": "Potencia Reactiva (MVAR)",
}
MAXLAGS = 12

# %% L1 — Importación, visualización y correlaciones
df = pd.read_excel(RUTA_DATA)
df["timestamp"] = pd.to_datetime(df["timestamp"])
df = df.set_index("timestamp")
df.index = pd.DatetimeIndex(df.index, freq="MS")

print("Primeras observaciones:")
print(df.head())
print(f"\nObservaciones: {len(df)} | Frecuencia índice: {df.index.freq}")

fig, ax = plt.subplots(figsize=(12, 6))
for col in VARIABLES:
    ax.plot(df.index, df[col], label=ETIQUETAS[col], linewidth=1.5)
ax.set_title("Evolución conjunta: Frecuencia, Potencia Activa y Reactiva")
ax.set_xlabel("Fecha")
ax.set_ylabel("Magnitud")
ax.legend(loc="best")
ax.grid(True, alpha=0.3)
plt.tight_layout()
ruta_series = RUTA_OUTPUTS / "l1_series_temporales.png"
fig.savefig(ruta_series, dpi=150)
plt.close()
print(f"\nGráfico L1 guardado: {ruta_series}")

corr = df[VARIABLES].corr()
print("\nMatriz de correlaciones:")
print(corr.round(4))

# %% L2 — VARselect, estimación y validación de estabilidad
datos_var = df[VARIABLES].dropna()
modelo = VAR(datos_var)
seleccion = modelo.select_order(maxlags=MAXLAGS)

print("\n--- VARselect ---")
print(seleccion.summary())
p_aic = seleccion.aic
p_bic = seleccion.bic
print(f"\nOrden AIC: p = {p_aic}")
print(f"Orden BIC: p = {p_bic}")

resultado_var = modelo.fit(p_aic)
print("\n--- Resumen VAR estimado ---")
print(resultado_var.summary())

raices = resultado_var.roots
modulos = np.abs(raices)
estable = np.all(modulos < 1)
print("\n--- Estabilidad ---")
print(f"Raíces del VAR companion: {raices}")
print(f"Módulos: {modulos.round(4)}")
print(f"Sistema estable (|λ| < 1): {estable}")

if __name__ == "__main__":
    print("\nLaboratorio L1-L2 completado.")
