import numpy as np
import statsmodels.api as sm
from statsmodels.tsa.arima.model import ARIMA
from statsmodels.tsa.stattools import adfuller
from statsmodels.tsa.vector_ar.vecm import coint_johansen


def adf_report(serie, nombre: str, regression: str = "c", maxlag: int = 5) -> dict:
    """Equivalente a ur.df() de urca con resumen impreso."""
    stat, pval, lags, nobs, crit, _ = adfuller(serie, maxlag=maxlag, regression=regression, autolag="BIC")
    print(f"\n--- Test ADF para {nombre} ---")
    print(f"Estadístico: {stat:.4f} | p-valor: {pval:.4f} | Rezagos: {lags}")
    print(f"Valores críticos: 1%={crit['1%']:.2f}, 5%={crit['5%']:.2f}, 10%={crit['10%']:.2f}")
    return {"stat": stat, "pval": pval, "lags": lags, "crit": crit}


def johansen_report(datos: np.ndarray, nombres: list[str], k_ar_diff: int = 1) -> coint_johansen:
    """Equivalente a ca.jo(type='trace', ecdet='const', K=2) de urca."""
    resultado = coint_johansen(datos, det_order=0, k_ar_diff=k_ar_diff)
    print("\n--- Test de Johansen (Trace) ---")
    print(f"Variables: {', '.join(nombres)}")
    for i in range(len(nombres)):
        print(
            f"H0: r<={i} | Trace={resultado.lr1[i]:.4f} "
            f"| 5%={resultado.cvt[i, 1]:.2f} | 1%={resultado.cvt[i, 2]:.2f}"
        )
    return resultado


def estimar_rango_johansen(resultado: coint_johansen, nivel: str = "5pct") -> int:
    """Determina r comparando estadístico trace con valores críticos."""
    idx = {"10pct": 0, "5pct": 1, "1pct": 2}[nivel]
    r = 0
    for i in range(len(resultado.lr1)):
        if resultado.lr1[i] > resultado.cvt[i, idx]:
            r = i + 1
    return r


def vectores_cointegracion(resultado: coint_johansen, r: int, normalizar_idx: int = 0) -> np.ndarray:
    """Extrae y normaliza vectores β (equivalente a johansen@V)."""
    beta = resultado.evec[:, :r].copy()
    for i in range(r):
        beta[:, i] /= beta[normalizar_idx, i]
    return beta


def ajuste_arima_aic(serie: np.ndarray, d: int = 1, max_p: int = 3, max_q: int = 3) -> tuple[ARIMA, np.ndarray]:
    """Selección AIC sobre rejilla pequeña (equivalente a auto.arima simplificado)."""
    mejor_aic, mejor_res, mejor_ord = np.inf, None, (0, d, 0)
    for p in range(max_p + 1):
        for q in range(max_q + 1):
            if p == 0 and q == 0:
                continue
            try:
                res = ARIMA(serie, order=(p, d, q)).fit()
                if res.aic < mejor_aic:
                    mejor_aic, mejor_res, mejor_ord = res.aic, res, (p, d, q)
            except Exception:
                continue
    print(f"ARIMA seleccionado: {mejor_ord} (AIC={mejor_aic:.2f})")
    return mejor_res, np.array(mejor_ord)
