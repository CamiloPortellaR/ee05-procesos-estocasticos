from pathlib import Path

import matplotlib.pyplot as plt

from core.plotting import asegurar_directorio_graficos


def generar_fir_ortogonalizada(
    resultado_var,
    periodos: int = 20,
    nombre_archivo: str = "fir_matriz_ortogonalizada.png",
) -> Path:
    """Calcula y guarda la matriz de FIR ortogonalizadas."""
    ruta = asegurar_directorio_graficos() / nombre_archivo
    irf = resultado_var.irf(periodos)
    fig = irf.plot(orth=True, figsize=(14, 10))
    fig.suptitle("FIR Ortogonalizadas — Matriz 3×3 (20 periodos)", fontsize=14, y=1.01)
    plt.tight_layout()
    fig.savefig(ruta, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return ruta
