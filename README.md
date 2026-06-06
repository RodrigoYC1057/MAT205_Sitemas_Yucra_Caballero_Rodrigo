# Método de Newton-Raphson — Cálculo de Momento Flector Máximo en Puente Peatonal

Proyecto de laboratorio de la materia **Métodos Numéricos** — Carrera de Ingeniería de Sistemas.  
Implementación del método numérico de Newton-Raphson para calcular el momento flector máximo en un puente peatonal ubicado en la Avenida Jaime Mendoza, ciudad de Sucre, Bolivia.

---

## Método Seleccionado

Se utilizó el **Método de Newton-Raphson** por su rápida convergencia hacia la solución, especialmente cuando se parte de una buena estimación inicial.

**Fórmula iterativa:**

```
x_(n+1) = x_n - f(x_n) / f'(x_n)
```

Donde:
- `x_n` — aproximación actual de la raíz
- `x_(n+1)` — nueva aproximación mejorada
- `f(x_n)` — valor de la función en `x_n`
- `f'(x_n)` — derivada de la función en `x_n`

El procedimiento se repite hasta que la diferencia entre `x_(n+1)` y `x_n` sea menor al criterio de tolerancia.

---

## Situación Problemática

El puente peatonal está diseñado para soportar una carga promedio de personas. Es necesario calcular el **momento flector máximo** que se generará bajo condiciones de carga inesperada, como la aglomeración masiva de personas durante eventos masivos (entrada de danzas folclóricas).

Por normativa, el peso admisible para puentes peatonales es de **350 kg/m²**, representado como carga rectangular distribuida.

---

## Planteamiento del Problema

¿Cómo aplicar el método de Newton-Raphson para el cálculo del momento flector máximo en el puente peatonal?

**Objeto de estudio:** Implementación del método numérico de Newton-Raphson para calcular el momento flector máximo en un puente peatonal de dos tramos:
- Tramo 1: **20 m**
- Tramo 2: **18 m**

---

## Objetivos

**General:** Calcular el momento flector máximo para estimar la carga máxima que soportará el puente peatonal.

**Específicos:**
- Establecer las ecuaciones de equilibrio que representan la distribución de momentos flectores.
- Desarrollar e implementar el algoritmo de Newton-Raphson para resolver las ecuaciones no lineales que describen el momento flector en el puente.

---

## Desarrollo — Aplicación del Método

### Datos del problema

- Carga distribuida: `w = 350 kg/m`
- Tramo 1: `L1 = 20 m`
- Tramo 2: `L2 = 18 m`
- Apoyos: A (extremo izquierdo), B (apoyo central), C (extremo derecho)

### Paso 1 — Ecuación de los 3 momentos

```
MA·LA + 2MB(L1 + L2) + MC·L2 + (A1·a1)/L1 + (A2·b1)/L2 = 6EI(h1/L1 + h3/L2)
```

Dado que `MA = 0`, `MC = 0`, `h1 = 0`, `h3 = 0`:

```
2MB(20 + 18) + 350(20)³/4 + 350(18)³/4 = 0
MB(76) + 1210300 = 0
MB = -15925 Kg·m
```

### Paso 2 — Cálculo de reacciones

**Tramo 1:**
- `V'A = 2703.75 Kg`
- `V'B = 4296.25 Kg`

**Tramo 2:**
- `V''B = 4034.722 Kg`
- `V'C = 2265.278 Kg`

**Reacciones totales:**
- `VA = 2703.75 Kg`
- `VB = 8330.972 Kg`
- `VC = 2265.278 Kg`

### Paso 3 — Ecuaciones de cortante y momento

**Tramo 1 (variable x):**
```
V = 2703.75 - 350x
M = 2703.75x - 350(x²/2)
```

**Tramo 2 (variable z):**
```
V = -2265.278 + 350z
M = 2265.278z - 350(z²/2)
```

### Paso 4 — Aplicando Newton-Raphson

**Tramo 1:**
```
f(x)  = 2703.75 - 350x
f'(x) = -350
X0 = 6

x1 = 6 - (2703.75 - 350·6) / (-350) = 7.725
```

**Tramo 2:**
```
f(x)  = -2265.278 + 350x
f'(x) = 350
X0 = -2

x1 = -2 - (-2265.278 + 350·(-2)) / 350 = 6.472
```

### Paso 5 — Prueba de control

```
2703.75 - 350(7.725) = 0  ✓
-2265.278 + 350(6.472) = 0  ✓
```

### Paso 6 — Momentos máximos

```
Mmax_1 = 2703.75(7.725) - 350(7.725²/2) = 10443.234 Kg·m
Mmax_2 = 2265.278(6.472) - 350(6.472²/2) = 7330.692 Kg·m
```

---

## Estructura del programa

```python
from sympy import symbols, log, diff, sin, cos, integrate, Rational
from sympy.parsing.sympy_parser import parse_expr

ERROR = 1E-6
x = symbols('x')

func = parse_expr(input('f(x): '))
x_0 = float(input('x_0: '))
dec = diff(func, x)

contador = 0
fx_1 = 100

while abs(fx_1) >= ERROR:
    fx_0 = func.subs(x, x_0).evalf()
    decx_0 = dec.subs(x, x_0).evalf()
    x_1 = (x_0) - (fx_0 / decx_0)
    fx_1 = func.subs(x, x_1).evalf()
    e = abs((x_1 - x_0) / x_1)
    print('iteración:', contador)
    print('x_0:', x_0)
    print('fx_0:', fx_0)
    print('x_1:', x_1)
    print('fx_1:', fx_1)
    print('*' * 15)
    if e > ERROR:
        x_0 = x_1
    contador += 1

print('Raíz: ', x_1)
print('f(Raíz):', fx_1)

# ─────────────────────────────────────────────────────────────────────────────
# CÁLCULO DEL MOMENTO MÁXIMO
# Se integra la función ingresada (cortante V) para obtener el momento M,
# luego se evalúa en la raíz encontrada (punto donde V = 0 → momento máximo).
# ─────────────────────────────────────────────────────────────────────────────

print()
print('=' * 40)
print('   CÁLCULO DEL MOMENTO MÁXIMO')
print('=' * 40)

# Integrar la función de cortante para obtener la función de momento
momento = integrate(func, x)
print(f'\nFunción de cortante V(x)  = {func}')
print(f'Función de momento   M(x) = {momento}')

# Evaluar el momento en la raíz (donde V = 0 → momento máximo)
momento_max = momento.subs(x, x_1).evalf()

print(f'\nRaíz encontrada (x donde V=0): x = {x_1}')
print(f'\nMomento Máximo M({float(x_1):.4f}) = {momento_max:.4f}')
print('=' * 40)
```

### Descripción del código

| Elemento | Descripción |
|----------|-------------|
| `sympy` | Librería para matemática simbólica (derivadas, logaritmos, etc.) |
| `parse_expr` | Convierte el texto ingresado por el usuario en expresión matemática |
| `ERROR = 1E-6` | Umbral de tolerancia para determinar convergencia |
| `diff(func, x)` | Calcula la derivada simbólica automáticamente |
| `func.subs(x, x_0).evalf()` | Evalúa la función en el punto actual |
| `x_1 = x_0 - fx_0/decx_0` | Aplica la fórmula iterativa de Newton-Raphson |
| `e = abs((x_1 - x_0) / x_1)` | Calcula el error relativo entre iteraciones |

---

## Dataset de prueba

| ID | Contexto / Aplicación | Función f(x) | x0 | Raíz Esperada |
|:--:|:---------------------:|:------------:|:--:|:-------------:|
| 01 | Ecuación lineal (caso del puente) | `2703.75 - 350*x` | 6 | 7.725 |
| 02 | Cálculo de tirante crítico (Hidráulica) | `x**3 - 5*x - 3` | 2.5 | 2.4908 |
| 03 | Esfuerzos en vigas (Logaritmos) | `x * log(x) - 10` | 5 | 5.7735 |
| 04 | Transiciones de peralte (Geometría Vial) | `x - 0.2*sin(x) - 0.5` | 1 | 0.6154 |
| 05 | Raíz cuadrada clásica | `x**2 - 25` | 10 | 5.0 |

---

## Diagrama del algoritmo

```
Inicio
  │
  ▼
Ingresar función f(x)
  │
  ▼
Ingresar x_0
  │
  ▼
Calcular f(x_n)
  │
  ▼
Calcular f'(x_n)
  │
  ▼
Aplicar Newton-Raphson: x_(n+1) = x_n - f(x_n)/f'(x_n)
  │
  ▼
Calcular error
  │
  ▼
¿Error < Tolerancia? ──No──→ Repetir
  │ Sí
  ▼
Mostrar raíz
```

---

## Conclusiones

- El uso del método de Newton-Raphson resultó efectivo para encontrar soluciones rápidas y precisas para el cálculo de momentos máximos en el puente peatonal.
- Este método iterativo permitió obtener resultados cercanos a la solución exacta con un número razonable de iteraciones.
- El éxito de su utilización depende en gran medida de una buena elección de la estimación inicial.
- El método se adapta bien a la resolución de ecuaciones no lineales derivadas de modelos estructurales.

## Recomendaciones

- Verificar los resultados con otros métodos numéricos o analíticos para asegurar que los momentos máximos calculados sean correctos.
- Elegir una estimación inicial cercana a la solución real basada en análisis preliminares.
- Implementar un criterio de convergencia claro y monitorear el número de iteraciones.

---

## Referencias

- Burden, R. L., & Faires, J. D. (2011). *Numerical Analysis* (9th ed.). Brooks/Cole.
- Chapra, S. C., & Canale, R. P. (2020). *Numerical Methods for Engineers* (8th ed.). McGraw-Hill Education.
- Hibbeler, R. C. (2017). *Structural Analysis* (10th ed.). Pearson Education.
- Beer, F. P., Johnston, E. R., DeWolf, J. T., & Mazurek, D. F. (2015). *Mechanics of Materials* (7th ed.). McGraw-Hill Education.

---

*UNIVERSIDAD MAYOR REAL Y PONTIFICIA DE SAN FRANCISCO XAVIER DE CHUQUISACA — Facultad Ciencias y Tecnología — 2026*
