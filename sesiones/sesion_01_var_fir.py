# %% Imports
import pandas as pd
import matplotlib.pyplot as plt
from statsmodels.tsa.api import VAR

# %% Carga de datos
RUTA_EXCEL = "../data/base_datos_var.xlsx"

try:
    df_potencia = pd.read_excel(RUTA_EXCEL)
    print(df_potencia.head())
except FileNotFoundError:
    print(f"Error: No se encontró el archivo en {RUTA_EXCEL}")
except ValueError as exc:
    print(f"Error al parsear el Excel: {exc}")
except Exception as exc:
    print(f"Error inesperado al leer la base de datos: {exc}")
