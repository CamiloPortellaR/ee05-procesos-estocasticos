from statsmodels.tsa.api import VAR


def seleccionar_orden(datos_var, maxlags: int = 12):
    """Equivalente a VARselect: retorna el objeto de selección y el orden AIC."""
    modelo = VAR(datos_var)
    seleccion = modelo.select_order(maxlags=maxlags)
    return modelo, seleccion, seleccion.aic


def estimar_var(modelo, orden: int):
    """Ajusta el VAR con el orden de rezagos indicado."""
    return modelo.fit(orden)
