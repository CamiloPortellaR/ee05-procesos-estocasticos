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

**ii) Efecto de $P_{t-1}$ sobre $f_t$:** El coeficiente $A_1(1,2)=-0{,}2<0$. Un aumento en potencia activa pasada se asocia con una **disminución** en frecuencia actual (controlando por los demás rezagos del sistema). Físicamente, refleja el acoplamiento dinámico del balance potencia–frecuencia: mayor demanda de potencia activa tiende a reducir la frecuencia si no hay compensación inmediata de generación.

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
| **AIC** | $2k - 2\ln(L)$ | Verosimilitud vs. complejidad (penaliza cada parámetro con $2$) | Al comparar modelos estimados sobre **los mismos datos**: **menor AIC → mejor modelo** (mayor ajuste sin sobreparametrizar en exceso). No tiene sentido interpretar el valor absoluto; solo el **orden relativo** entre candidatos. |
| **BIC** | $k\ln(n) - 2\ln(L)$ | Verosimilitud vs. complejidad (penaliza cada parámetro con $\ln(n)$) | Al comparar modelos sobre **los mismos datos**: **menor BIC → mejor modelo**. La penalización crece con el tamaño muestral $n$, por lo que, con $n$ grande, BIC suele elegir **menos parámetros** que AIC ante el mismo conjunto de candidatos. |
| **MAPE** | $\frac{100}{n}\sum\left|\frac{Y_i-\hat Y_i}{Y_i}\right|$ | Error de pronóstico relativo (%) | **Menor MAPE → mejor capacidad predictiva** (menor error porcentual medio entre observado y pronosticado). |

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

![Series temporales superpuestas](../outputs/l1_series_temporales.png)

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

### L3 (10 %) — Función Impulso-Respuesta (FIR)

#### a) FIR ortogonalizadas (20 periodos, Cholesky)

Sobre el **VAR(2)** estimado en L2, se calcularon las FIR ortogonalizadas con descomposición de Cholesky de $\Sigma$ (orden de variables: frecuencia → potencia activa → potencia reactiva).

![FIR ortogonalizadas — matriz 3×3](../outputs/l3_fir_matriz_ortogonalizada.png)

#### b) Matriz de gráficos 3×3

La figura anterior muestra la respuesta de cada variable (filas) ante un shock unitario ortogonalizado en cada variable (columnas), para $h=0,1,\ldots,20$.

#### c) Interpretación: shock en potencia activa sobre frecuencia

Respuesta de $f_t$ ante shock en $P_t$ (FIR ortogonalizada, variable respuesta = frecuencia, shock = P. Activa):

| Horizonte $h$ | Respuesta |
|---:|---:|
| 0 | 0,000 |
| 1 | $-0{,}052$ |
| 5 | $-0{,}110$ |
| 10 | $-0{,}078$ |
| 20 | $-0{,}027$ |

**Conclusión:** El efecto es **negativo** (mayor potencia activa pasada/shock se traduce en menor frecuencia en el corto plazo). La respuesta es **transitoria**: se intensifica en los primeros periodos (máximo alrededor de $h=5$) y luego decae hacia cero, sin comportamiento explosivo persistente.

**Nota metodológica:** El VAR(2) de L2 es **inestable** ($|\lambda|>1$ en algunas raíces). Las FIR se reportan como requerimiento del examen, pero en un análisis operativo convendría reestimar sobre datos estacionarios o con restricciones.

#### Código L3

```python
# %% L3 — FIR ortogonalizadas
irf = resultado_var.irf(20)
fig = irf.plot(orth=True, figsize=(14, 10))
fig.savefig(RUTA_OUTPUTS / "l3_fir_matriz_ortogonalizada.png", dpi=150, bbox_inches="tight")

orth = irf.orth_irfs
respuesta_f_a_shock_p = orth[:, 0, 1]  # fila=frecuencia, shock=P. Activa
```

---

### L4 (12 %) — Cointegración y VECM

#### a) Test de Johansen (estadístico de traza)

Parámetro de rezagos en diferencias: $k_{\Delta}=1$ (equivalente a $p-1$ del VAR(2)).

| Hipótesis nula | Estadístico traza | Valor crítico 5 % | Decisión (5 %) |
|---|---:|---:|---|
| $r \leq 0$ | 795,45 | 29,80 | Rechazar |
| $r \leq 1$ | 51,13 | 15,49 | Rechazar |
| $r \leq 2$ | 5,90 | 3,84 | Rechazar |

**Rango de cointegración según traza (5 %):** $r = 3$.

Con $K=3$ variables, $r=3$ implica que el test detecta el máximo de relaciones posibles; en la práctica se interpreta la **relación principal** mediante el primer vector de cointegración.

#### b) VECM y vector de cointegración normalizado

Se estimó un **VECM** con $r=2$ (máximo no degenerado para $K=3$) y $k_{\Delta}=1$.

**Primer vector de cointegración (Johansen, normalizado en frecuencia):**

$$\beta' Y_t \approx f_t + 0{,}681\,P_t + 0{,}085\,Q_t$$

**Interpretación en sistema de potencia:** Existe una combinación lineal de frecuencia, potencia activa y potencia reactiva que es estacionaria (equilibrio de largo plazo). Las tres magnitudes eléctricas no divergen de forma independiente: un desvío del equilibrio $(f,P,Q)$ genera fuerzas de corrección vía los coeficientes de ajuste $\alpha$ del VECM (la frecuencia presenta velocidad de ajuste significativa: $\alpha_{f,ec1} \approx -0{,}15$).

#### c) Pronóstico a 12 pasos: VAR vs VECM

Se reservaron las **últimas 12 observaciones** como conjunto de prueba. Se reestimaron VAR(2) y VECM($r=2$) sobre el resto y se comparó el RMSE.

| Variable | RMSE VAR | RMSE VECM |
|---|---:|---:|
| Frecuencia | 0,151 | **0,131** |
| P. Activa | **0,675** | 0,755 |
| P. Reactiva | 1,017 | **1,002** |
| **Promedio global** | **0,614** | 0,630 |

![Pronóstico VAR vs VECM (12 meses)](../outputs/l4_pronostico_var_vs_vecm.png)

**Conclusión:** El **VAR sin restricciones** presenta menor RMSE medio global (0,614 vs 0,630). Sin embargo, el **VECM mejora el pronóstico de frecuencia** (0,131 vs 0,151), lo cual es coherente con que impone equilibrio de largo plazo entre variables acopladas. El VAR gana en potencia activa al no restringir la dinámica de largo plazo. En operación SCADA/analítica, la elección depende de la variable objetivo: si el foco es **frecuencia**, el VECM aporta valor; si se busca error global mínimo en las tres series, el VAR es preferible en esta muestra.

#### Código L4

```python
# %% L4 — Johansen, VECM y pronóstico
from statsmodels.tsa.vector_ar.vecm import VECM, coint_johansen

k_ar_diff = max(p_aic - 1, 1)
johansen = coint_johansen(datos_var.values, det_order=0, k_ar_diff=k_ar_diff)
r_traza = sum(johansen.lr1[i] > johansen.cvt[i, 1] for i in range(3))

beta = johansen.evec[:, 0].real
beta /= beta[0]

vecm_res = VECM(datos_var, k_ar_diff=k_ar_diff, coint_rank=2, deterministic="ci").fit()

train, test = datos_var.iloc[:-12], datos_var.iloc[-12:]
pron_var = VAR(train).fit(p_aic).forecast(train.values[-p_aic:], steps=12)
pron_vecm = VECM(train, k_ar_diff=k_ar_diff, coint_rank=2, deterministic="ci").fit().predict(steps=12)
rmse_var = np.sqrt(np.mean((test.values - pron_var) ** 2, axis=0))
rmse_vecm = np.sqrt(np.mean((test.values - pron_vecm) ** 2, axis=0))
```

---

### Código completo de referencia

Ver: `../scripts/resolucion_laboratorio.py`

---

*Resolución generada para la Sesión 06 — Examen Parcial UPG FIEE 2026-1.*
