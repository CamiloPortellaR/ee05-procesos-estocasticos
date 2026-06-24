# %% Configuración de rutas del mini-proyecto
import sys
from pathlib import Path

RUTA_SESION = Path(__file__).resolve().parent
if str(RUTA_SESION) not in sys.path:
    sys.path.insert(0, str(RUTA_SESION))

from core.config import MAXLAGS, PERIODOS_FIR, RUTA_GRAFICOS
from core.data_loader import cargar_datos_potencia, obtener_series_var
from core.plotting import plot_series_temporales
from models.impulse_response import generar_fir_ortogonalizada
from models.var_estimator import estimar_var, seleccionar_orden

# %% Requerimiento 1 — Importación y visualización
df = cargar_datos_potencia()

print("Primeras observaciones:")
print(df.head())
print(f"\nFrecuencia del índice: {df.index.freq}")
print(f"Observaciones: {len(df)}")

ruta_series = plot_series_temporales(df)
print(f"Gráfico guardado: {ruta_series}")

# %% Requerimiento 2 — Selección de rezagos y estimación VAR
datos_var = obtener_series_var(df)
modelo, seleccion, p_aic = seleccionar_orden(datos_var, maxlags=MAXLAGS)

print("\n--- Selección de orden (VARselect equivalente) ---")
print(seleccion.summary())
print(f"\nRezago óptimo según AIC: p = {p_aic}")

resultado_var = estimar_var(modelo, p_aic)
print("\n--- Resumen del modelo VAR estimado ---")
print(resultado_var.summary())

# %% Requerimiento 3 — FIR ortogonalizadas
ruta_fir = generar_fir_ortogonalizada(resultado_var, periodos=PERIODOS_FIR)
print(f"\nFIR ortogonalizadas calculadas para {PERIODOS_FIR} periodos.")
print("Gráficos guardados en:", RUTA_GRAFICOS)
