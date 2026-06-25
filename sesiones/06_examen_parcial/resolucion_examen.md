# Examen Parcial — Procesos Estocásticos — Resolución

**Curso:** Procesos Estocásticos 2026-1 | **UPG FIEE**  
**Fecha del examen:** 25 de junio de 2026  
**Fecha de resolución:** 25 de junio de 2026  
**Proyecto:** ee05-procesos-estocasticos | **Sesión:** 06_examen_parcial

---

## Parte Teórica (60 %)

---

### Pregunta 1 (15 %) — Fundamentos de Procesos Estocásticos y Modelos Univariados

#### a) Proceso estocástico y estacionariedad (4 puntos)

Un **proceso estocástico** $\{Y_t : t \in T\}$ es una colección de variables aleatorias indexadas en el tiempo. Cada realización del proceso es una **trayectoria** o serie temporal.

| Tipo | Definición |
|---|---|
| **Estacionario en sentido estricto** | Todas las distribuciones finito-dimensionales son invariantes ante desplazamientos temporales: la distribución conjunta de $(Y_{t_1},\ldots,Y_{t_k})$ coincide con la de $(Y_{t_1+h},\ldots,Y_{t_k+h})$ para todo $h$. |
| **Estacionario en sentido débil (covarianza estacionaria)** | $\mathbb{E}[Y_t]=\mu$ constante, $\text{Var}(Y_t)=\sigma^2$ finita y $\text{Cov}(Y_t,Y_{t-k})=\gamma(k)$ depende solo del rezago $k$. |

**Diferencia fundamental:** la estacionariedad estricta implica la débil (bajo momentos finitos), pero no al revés. Un proceso puede tener media y autocovarianza constantes sin que todas sus distribuciones sean invariantes.

#### b) Ruido blanco (4 puntos)

Un proceso $\{\varepsilon_t\}$ es **ruido blanco** si:

$$\mathbb{E}[\varepsilon_t]=0, \quad \text{Var}(\varepsilon_t)=\sigma^2, \quad \text{Cov}(\varepsilon_t,\varepsilon_s)=0 \;\; (t\neq s)$$

Es **fundamental** porque los modelos ARIMA/AR/MA/VAR asumen que la parte no explicada por la dinámica propia es innovación no predecible. Sin ruido blanco no hay base para inferencia, pronóstico ni interpretación causal de shocks.

#### c) Proceso AR(1) (7 puntos)

$$y_t = \phi_0 + \phi_1 y_{t-1} + \varepsilon_t, \quad \varepsilon_t \sim WN(0,\sigma^2)$$

**i) Condición de estacionariedad:** $|\phi_1| < 1$ (raíz característica fuera del círculo unitario).

**ii) Media y varianza en estado estacionario:**

$$\mu = \frac{\phi_0}{1-\phi_1}, \qquad \gamma(0) = \frac{\sigma^2}{1-\phi_1^2}$$

**iii) Función de autocorrelación:**

$$\rho(k) = \phi_1^k$$

La memoria del proceso decae geométricamente con el rezago.

---

### Pregunta 2 (15 %) — Modelos VAR y Función Impulso-Respuesta (FIR)

#### a) Definición formal del VAR($p$) (4 puntos)

$$Y_t = c + A_1 Y_{t-1} + A_2 Y_{t-2} + \cdots + A_p Y_{t-p} + \varepsilon_t$$

| Componente | Dimensión | Descripción |
|---|---|---|
| $Y_t$ | $K \times 1$ | Vector de variables endógenas |
| $c$ | $K \times 1$ | Vector de constantes |
| $A_i$ | $K \times K$ | Matrices de coeficientes del rezago $i$ |
| $\varepsilon_t$ | $K \times 1$ | Innovaciones: $E[\varepsilon_t]=0$, $E[\varepsilon_t\varepsilon_t']=\Sigma$ |

#### b) Función Impulso-Respuesta (5 puntos)

Mediante la representación VMA($\infty$):

$$Y_t = \mu + \sum_{h=0}^{\infty} \Psi_h \varepsilon_{t-h}, \quad \Psi_0 = I_K$$

La **FIR** es la secuencia $\{\Psi_h\}$. El elemento $\Psi_h(i,k)$ mide el efecto sobre $Y_{i,t+h}$ de un shock unitario en $\varepsilon_{k,t}$.

#### c) Sistema de potencia VAR(1) (6 puntos)

$$A_1 = \begin{bmatrix} 0{,}6 & -0{,}2 & 0{,}1 \\ 0{,}1 & 0{,}8 & -0{,}1 \\ 0{,}0 & 0{,}1 & 0{,}7 \end{bmatrix}$$

**i) Estabilidad:** Los autovalores son $\lambda \approx 0{,}7 \pm 0{,}14j$ y $\lambda_3 = 0{,}7$. Todos cumplen $|\lambda_i|<1$, por lo tanto el VAR(1) es **estable**.

**ii) Efecto de $P_{t-1}$ sobre $f_t$:** El coeficiente $A_1(1,2)=0{,}1>0$. Un aumento en potencia activa pasada se asocia con un incremento en frecuencia actual (controlando por los demás rezagos del sistema). Físicamente, refleja acoplamiento dinámico entre balance de potencia activa y dinámica de frecuencia del sistema.

**iii) Cholesky:** Como $\Sigma$ no es diagonal, los shocks contemporáneos están correlacionados. La descomposición $\Sigma = PP'$ define shocks ortogonalizados $u_t = P^{-1}\varepsilon_t$ y FIR ortogonalizadas $\Theta_h = \Psi_h P$, aislando el efecto de una variable sin contaminación simultánea de las demás.

---

### Pregunta 3 (15 %) — Cointegración y Aplicaciones en Sistemas de Potencia

#### a) Concepto de cointegración (5 puntos)

Un vector $X_t$ es $I(1)$ si cada componente requiere una diferencia para ser estacionaria. Es **cointegrado** si existe $\beta$ tal que $\beta' X_t \sim I(0)$: una combinación lineal de variables no estacionarias es estacionaria (relación de equilibrio de largo plazo).

#### b) Aplicación en planificación de expansión de red (5 puntos)

Variables: demanda máxima $D_t$, capacidad instalada $G_t$, capacidad de transmisión $L_t$, inversión $I_t$. Una relación $\beta'[D_t,G_t,L_t,I_t]' \sim I(0)$ expresa que, aunque cada variable crece con tendencia, el sistema revierte a un equilibrio oferta-demanda-infraestructura. El ECT guía decisiones de inversión en líneas y generación.

#### c) Pérdida de cointegración en monitoreo de frecuencia (5 puntos)

Si dos nodos distantes pierden cointegración en frecuencia, el equilibrio dinámico de largo plazo entre áreas se rompe. Implica posible **desacople operativo**, mayor riesgo de oscilaciones interárea, necesidad de revisar esquemas de control y posible intervención del operador (redespacho, separación de áreas).

---

### Pregunta 4 (15 %) — Criterios de Selección y Aplicaciones

#### a) AIC, BIC y MAPE (6 puntos)

| Criterio | Fórmula | Balance | Interpretación |
|---|---|---|---|
| **AIC** | $2k - 2\ln(L)$ | Ajuste vs. parsimonia (penalización moderada) | Menor AIC → mejor modelo |
| **BIC** | $k\ln(n) - 2\ln(L)$ | Penalización creciente con $n$ y $k$ | Favorece modelos más simples que AIC |
| **MAPE** | $\frac{100}{n}\sum\left|\frac{Y_i-\hat Y_i}{Y_i}\right|$ | Error de pronóstico relativo | Menor MAPE → mejor capacidad predictiva |

#### b) Interpretación FIR potencia-frecuencia (4 puntos)

Caída inmediata de $f_t$ tras shock en $P_t$ y recuperación lenta en 20 s indica **acoplamiento dinámico fuerte** con **inercia del sistema**: el desbalance de potencia afecta la frecuencia al instante, pero la regulación primaria/secundaria devuelve el equilibrio de forma gradual (respuesta persistente pero no explosiva).

#### c) VAR(1) temperatura-presión-caudal (5 puntos)

$$T_t = c_1 + a_{11}T_{t-1} + a_{12}P_{t-1} + a_{13}C_{t-1} + \varepsilon_{1t}$$

$$P_t = c_2 + a_{21}T_{t-1} + a_{22}P_{t-1} + a_{23}C_{t-1} + \varepsilon_{2t}$$

$$C_t = c_3 + a_{31}T_{t-1} + a_{32}P_{t-1} + a_{33}C_{t-1} + \varepsilon_{3t}$$

Si $A_1(2,1)=a_{21}=0$, la **temperatura rezagada no entra directamente** en la ecuación de presión. El lazo de control de presión puede diseñarse sin $T_{t-1}$ como regresor directo (efectos indirectos vía caudal siguen posibles si $a_{23}\neq 0$).

---

## Parte de Laboratorio (40 %)

### Contexto

Base de datos: `data/base_datos_var.xlsx` — 100 observaciones mensuales (`MS`) con:

- `frecuencia_hz` — Frecuencia (Hz)
- `potencia_activa_mw` — Potencia Activa (MW)
- `potencia_reactiva_mvar` — Potencia Reactiva (MVAR)

---

### L1 (8 %) — Importación, visualización y análisis exploratorio

Se cargó el Excel, se configuró índice temporal mensual y se generó el gráfico de series superpuestas.

![Series temporales superpuestas](outputs/l1_series_temporales.png)

**Matriz de correlaciones estimada:**

|  | Frecuencia | P. Activa | P. Reactiva |
|---|---:|---:|---:|
| Frecuencia | 1,000 | **0,995** | 0,980 |
| P. Activa | 0,995 | 1,000 | 0,988 |
| P. Reactiva | 0,980 | 0,988 | 1,000 |

**Interpretación:** Existe correlación lineal muy alta entre frecuencia y potencia activa ($\rho \approx 0{,}995$), coherente con el acoplamiento físico del sistema eléctrico: cuando la potencia activa cae, la frecuencia tiende a caer en este escenario muestral (degradación conjunta visible desde 2024-03).

#### Código L1

```python
# %% L1 — Importación, visualización y correlaciones
from pathlib import Path
import matplotlib.pyplot as plt
import pandas as pd

RUTA_SESION = Path(__file__).resolve().parent.parent
RUTA_DATA = RUTA_SESION / "data" / "base_datos_var.xlsx"
RUTA_OUTPUTS = RUTA_SESION / "outputs"
VARIABLES = ["frecuencia_hz", "potencia_activa_mw", "potencia_reactiva_mvar"]

df = pd.read_excel(RUTA_DATA)
df["timestamp"] = pd.to_datetime(df["timestamp"])
df = df.set_index("timestamp")
df.index = pd.DatetimeIndex(df.index, freq="MS")

fig, ax = plt.subplots(figsize=(12, 6))
for col in VARIABLES:
    ax.plot(df.index, df[col], label=col, linewidth=1.5)
ax.set_title("Evolución conjunta: Frecuencia, Potencia Activa y Reactiva")
ax.legend(); ax.grid(True, alpha=0.3)
fig.savefig(RUTA_OUTPUTS / "l1_series_temporales.png", dpi=150)

corr = df[VARIABLES].corr()
print(corr.round(4))
```

---

### L2 (10 %) — Estimación del modelo VAR

#### a) Selección de rezagos (VARselect)

| Criterio | Orden óptimo $p$ |
|---|---:|
| **AIC** | **2** |
| **BIC** | **1** |

#### b) Estimación VAR($\hat{p}$) con $\hat{p}=2$ (AIC)

Se estimó **VAR(2)** con 98 observaciones efectivas. Coeficientes destacados en la ecuación de frecuencia:

- Persistencia propia: $\hat{\alpha}_{11}^{(1)} = 0{,}651$, $\hat{\alpha}_{11}^{(2)} = 0{,}193$
- Efecto de potencia activa rezagada: $\hat{\alpha}_{12}^{(1)} = -0{,}062$, $\hat{\alpha}_{12}^{(2)} = -0{,}039$ (significativos)

#### c) Estabilidad

Los módulos de las raíces del companion matrix son $| \lambda | \in \{5{,}28,\; 2{,}37,\; 1{,}38,\; 1{,}14\}$ — **ninguno menor que 1 en todos los casos superiores a 1**. El VAR(2) estimado sobre esta muestra **no cumple estabilidad** en sentido estricto, señal de posible no estacionariedad o quiebre estructural en los datos (caída abrupta de frecuencia y potencia).

#### Código L2

```python
# %% L2 — VARselect, estimación y validación
from statsmodels.tsa.api import VAR
import numpy as np

datos_var = df[VARIABLES].dropna()
modelo = VAR(datos_var)
seleccion = modelo.select_order(maxlags=12)
p_aic, p_bic = seleccion.aic, seleccion.bic

resultado_var = modelo.fit(p_aic)
print(seleccion.summary())
print(resultado_var.summary())

raices = resultado_var.roots
print("Módulos raíces:", np.abs(raices))
print("Estable:", np.all(np.abs(raices) < 1))
```

---

### Nota sobre L3 y L4

Los requerimientos L3 (FIR ortogonalizadas) y L4 (Johansen, VECM, pronóstico) se desarrollan en la sesión 04 (`04_vectores_autorregresivos`) y sesión 05 (`05_cointegracion_vecm`) del laboratorio EE-05. Dada la **inestabilidad** del VAR(2) sobre esta muestra, se recomienda diferenciar o restringir el modelo antes de FIR y VECM en un análisis operativo.

---

### Código completo de referencia

Ver: `scripts/resolucion_laboratorio.py`

---

*Resolución generada para la Sesión 06 — Examen Parcial UPG FIEE 2026-1.*
