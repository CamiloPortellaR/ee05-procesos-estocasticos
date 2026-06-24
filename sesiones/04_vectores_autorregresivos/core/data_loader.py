import pandas as pd

from .config import RUTA_EXCEL, VARIABLES


def cargar_datos_potencia() -> pd.DataFrame:
    """Carga y preprocesa la base de datos VAR local (data/base_datos_var.xlsx)."""
    try:
        df = pd.read_excel(RUTA_EXCEL)
    except FileNotFoundError as exc:
        raise FileNotFoundError(f"No se encontró el archivo en {RUTA_EXCEL}") from exc
    except ValueError as exc:
        raise ValueError(f"Error al parsear el Excel: {exc}") from exc

    df["timestamp"] = pd.to_datetime(df["timestamp"])
    df = df.set_index("timestamp")
    df.index = pd.DatetimeIndex(df.index, freq="MS")
    return df


def obtener_series_var(df: pd.DataFrame) -> pd.DataFrame:
    """Retorna las columnas endógenas del VAR sin valores faltantes."""
    return df[VARIABLES].dropna()
