import matplotlib.pyplot as plt

from .config import INTERACTIVO, RUTA_GRAFICOS


def guardar_figura(fig, nombre: str) -> None:
    """Guarda figura en graficos/; en modo interactivo también la muestra inline."""
    fig.savefig(RUTA_GRAFICOS / nombre, dpi=150, bbox_inches="tight")
    if INTERACTIVO:
        plt.show()
    else:
        plt.close(fig)
