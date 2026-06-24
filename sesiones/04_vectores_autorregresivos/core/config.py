from pathlib import Path

RUTA_SESION = Path(__file__).resolve().parent.parent
RUTA_DATA = RUTA_SESION / "data"
RUTA_EXCEL = RUTA_DATA / "base_datos_var.xlsx"
RUTA_GRAFICOS = RUTA_SESION / "graficos"

VARIABLES = ["frecuencia_hz", "potencia_activa_mw", "potencia_reactiva_mvar"]
ETIQUETAS = {
    "frecuencia_hz": "Frecuencia (Hz)",
    "potencia_activa_mw": "Potencia Activa (MW)",
    "potencia_reactiva_mvar": "Potencia Reactiva (MVAR)",
}

MAXLAGS = 12
PERIODOS_FIR = 20
