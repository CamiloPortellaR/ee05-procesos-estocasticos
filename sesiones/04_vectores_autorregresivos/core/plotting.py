from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from .config import ETIQUETAS, RUTA_GRAFICOS, VARIABLES


def asegurar_directorio_graficos() -> Path:
    RUTA_GRAFICOS.mkdir(parents=True, exist_ok=True)
    return RUTA_GRAFICOS


def plot_series_temporales(df: pd.DataFrame, nombre_archivo: str = "series_temporales_superpuestas.png") -> Path:
    """Grafica la evolución conjunta de las variables del sistema de potencia."""
    ruta = asegurar_directorio_graficos() / nombre_archivo
    fig, ax = plt.subplots(figsize=(12, 6))
    for col in VARIABLES:
        ax.plot(df.index, df[col], label=ETIQUETAS[col], linewidth=1.5)
    ax.set_title("Evolución conjunta: Frecuencia, Potencia Activa y Reactiva")
    ax.set_xlabel("Fecha")
    ax.set_ylabel("Magnitud")
    ax.legend(loc="best")
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(ruta, dpi=150)
    plt.close(fig)
    return ruta
