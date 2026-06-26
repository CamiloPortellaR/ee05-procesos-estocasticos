# %% Configuración
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from statsmodels.tsa.api import VAR
from statsmodels.tsa.vector_ar.vecm import VECM, coint_johansen

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
NOMBRES_CORTOS = ["Frecuencia", "P. Activa", "P. Reactiva"]
MAXLAGS = 12
PERIODOS_IRF = 20
HORIZONTE_PRONOSTICO = 12


def estimar_rango_traza(johansen: coint_johansen, nivel: str = "5pct") -> int:
    idx = {"10pct": 0, "5pct": 1, "1pct": 2}[nivel]
    r = 0
    for i in range(len(johansen.lr1)):
        if johansen.lr1[i] > johansen.cvt[i, idx]:
            r = i + 1
    return r


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

# %% L3 — FIR ortogonalizadas (Cholesky) y matriz 3×3
irf = resultado_var.irf(PERIODOS_IRF)
fig = irf.plot(orth=True, figsize=(14, 10))
fig.suptitle(
    "FIR Ortogonalizadas — Matriz 3×3 (20 periodos, Cholesky)",
    fontsize=14,
    y=1.01,
)
plt.tight_layout()
ruta_fir = RUTA_OUTPUTS / "l3_fir_matriz_ortogonalizada.png"
fig.savefig(ruta_fir, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"\nGráfico L3 guardado: {ruta_fir}")

orth = irf.orth_irfs
idx_f, idx_p = 0, 1
respuesta_f_a_shock_p = orth[:, idx_f, idx_p]
print("\n--- Respuesta de frecuencia a shock en P. Activa (ortogonalizado) ---")
for h in [0, 1, 5, 10, PERIODOS_IRF]:
    print(f"h={h:2d}: {respuesta_f_a_shock_p[h]: .4f}")

# %% L4 — Johansen, VECM y comparación de pronósticos
k_ar_diff = max(p_aic - 1, 1)
johansen = coint_johansen(datos_var.values, det_order=0, k_ar_diff=k_ar_diff)
r_traza = estimar_rango_traza(johansen, nivel="5pct")

print("\n--- Test de Johansen (estadístico de traza) ---")
for i in range(len(VARIABLES)):
    print(
        f"H0: r<={i} | Trace={johansen.lr1[i]:.4f} "
        f"| 5%={johansen.cvt[i, 1]:.2f} "
        f"| Rechaza: {johansen.lr1[i] > johansen.cvt[i, 1]}"
    )
print(f"\nRango de cointegración por traza (5%): r = {r_traza}")

beta_johansen = johansen.evec[:, 0].real.copy()
beta_johansen /= beta_johansen[0]
print("\nVector de cointegración (1.er vector, normalizado en frecuencia):")
print(beta_johansen.round(4))

r_vecm = min(max(r_traza, 1), len(VARIABLES) - 1)
vecm_res = VECM(
    datos_var,
    k_ar_diff=k_ar_diff,
    coint_rank=r_vecm,
    deterministic="ci",
).fit()
print(f"\n--- VECM estimado (r={r_vecm}, k_ar_diff={k_ar_diff}) ---")
print(vecm_res.summary())
print("\nMatriz β (VECM):")
print(vecm_res.beta.round(4))

train = datos_var.iloc[:-HORIZONTE_PRONOSTICO]
test = datos_var.iloc[-HORIZONTE_PRONOSTICO:]
var_train = VAR(train).fit(p_aic)
pronostico_var = var_train.forecast(train.values[-p_aic:], steps=HORIZONTE_PRONOSTICO)
rmse_var = np.sqrt(np.mean((test.values - pronostico_var) ** 2, axis=0))

vecm_train = VECM(
    train,
    k_ar_diff=k_ar_diff,
    coint_rank=r_vecm,
    deterministic="ci",
).fit()
pronostico_vecm = vecm_train.predict(steps=HORIZONTE_PRONOSTICO)
rmse_vecm = np.sqrt(np.mean((test.values - pronostico_vecm) ** 2, axis=0))

print("\n--- Comparación de pronósticos a 12 pasos (últimas 12 obs. como test) ---")
for i, nombre in enumerate(NOMBRES_CORTOS):
    print(
        f"{nombre:12s} | RMSE VAR={rmse_var[i]:.4f} "
        f"| RMSE VECM={rmse_vecm[i]:.4f}"
    )
print(f"RMSE medio global VAR : {rmse_var.mean():.4f}")
print(f"RMSE medio global VECM: {rmse_vecm.mean():.4f}")

fig, axes = plt.subplots(1, 3, figsize=(15, 4), sharex=True)
fechas_test = test.index
for i, (ax, nombre) in enumerate(zip(axes, NOMBRES_CORTOS)):
    ax.plot(fechas_test, test.iloc[:, i], "k-o", label="Real", markersize=4)
    ax.plot(fechas_test, pronostico_var[:, i], "--", label="VAR", linewidth=1.2)
    ax.plot(fechas_test, pronostico_vecm[:, i], ":", label="VECM", linewidth=1.2)
    ax.set_title(nombre)
    ax.grid(True, alpha=0.3)
    if i == 0:
        ax.legend(fontsize=8)
axes[0].set_ylabel("Magnitud")
axes[1].set_xlabel("Fecha")
plt.suptitle("Pronóstico a 12 meses: VAR vs VECM", y=1.02)
plt.tight_layout()
ruta_pronostico = RUTA_OUTPUTS / "l4_pronostico_var_vs_vecm.png"
fig.savefig(ruta_pronostico, dpi=150, bbox_inches="tight")
plt.close(fig)
print(f"\nGráfico L4 guardado: {ruta_pronostico}")

if __name__ == "__main__":
    print("\nLaboratorio L1-L4 completado.")
