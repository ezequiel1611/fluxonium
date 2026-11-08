# Fluxonium: mediciones de resonadores de readout

Este repositorio reúne el código, los notebooks y los datos de las mediciones del chip **fluxonium** del grupo, hechas en octubre y noviembre de 2026. Sirven para dos cosas:

- **Medir.** El código controla el VNA y la fuente de la bobina, y se puede reutilizar para mediciones parecidas.
- **Postprocesar.** Las bases de QCoDeS tienen los datos crudos y toda la metadata: potencia, atenuación, corriente y estado del VNA.

## Estructura

```
fluxonium/
├── codigo/       módulos de Python (control de equipos y ajustes)
├── notebooks/    notebooks de medición y análisis
├── datos/        bases de datos QCoDeS (SQLite, por Git LFS) y tablas CSV
├── informes/     informes y figuras
└── requirements.txt
```

## Cómo usarlo

1. **Clonar.** Las bases están en Git LFS. Hay que tener [git-lfs](https://git-lfs.com) instalado antes de clonar; si no, se bajan solo los punteros:
   ```bash
   git lfs install
   git clone https://github.com/ezequiel1611/fluxonium
   ```
2. **Instalar el entorno.** Es Python 3.10 con las versiones de `requirements.txt`, que son las del entorno conda `quantum` de la PC del laboratorio:
   ```bash
   pip install -r requirements.txt
   ```
3. **Abrir un notebook.** La primera celda de código de cada notebook agrega `codigo/` al path y se para en `datos/`. Por eso las rutas de las bases (por ejemplo `'fluxonium_2/punch_out.db'`) funcionan igual que en el laboratorio. Para postprocesar no hace falta tener los equipos: se saltean las celdas de conexión y se corren las de carga y análisis.

Ejemplo mínimo de carga, desde `datos/`:

```python
import sys; sys.path.insert(0, '../codigo')
from qcodes.dataset import initialise_or_create_database_at
from barrido_flujo import cargar_barrido_corriente, ajustar_barrido_corriente, saltos

initialise_or_create_database_at('fluxonium_2/punch_out.db')
d = cargar_barrido_corriente(1)        # run 1: barrido en corriente (ida)
a = ajustar_barrido_corriente(d)       # fr, Ql, Qi, Qc, κ por corriente
print(saltos(a['I_uA'], a['fr']))      # saltos periódicos de fr
```

## Montaje experimental

| | |
|---|---|
| VNA | Keysight PNA N5245A, `TCPIP::192.168.1.6::inst0::INSTR` (driver de QCoDeS) |
| Bobina (flujo) | Fuente de corriente Keithley, `GPIB0::15::INSTR` (comandos `I{valor}X` y `CURR?`) |
| Generador (two-tone) | Keysight N5173B, usado en `readout.ipynb` |
| Atenuación de entrada | 50 dB en la línea del criostato + 62 dB externos (dos atenuadores a la salida del VNA) = **−112 dB**. Algunas mediciones viejas usan −81 dB: ver `att_vna` en la metadata de cada run. |

Número medio de fotones en el resonador: ⟨n⟩ = 2 P_in Ql² / (Qc ħ ω0²), con `fit_resonancia.numero_fotones`. Para res_q4 con −112 dB de entrada:

| P VNA | ⟨n⟩ |
|---|---|
| −30 dBm | ≈ 0.02 |
| −10 dBm | ≈ 1.7 |
| 0 dBm | ≈ 17 |
| +10 dBm | ≈ 170 |

## Código (`codigo/`)

| Archivo | Equipo / uso | Contenido principal |
|---|---|---|
| `vna_single.py` | VNA N5245A | `configurar_vna`, `medir_single` (barrido simple guardado en QCoDeS con S21 complejo y metadata del VNA), `barrido_potencia` y `buscar_barrido` (barridos en potencia, un run por potencia), `leer_estado`, `traza_s21` |
| `barrido_flujo.py` | Fuente Keithley + VNA | `conectar_fuente` / `FuenteCorriente` (rampas seguras con límite `I_max` y `ir_a` para acercarse siempre desde el mismo lado), `barrido_corriente` (mapa 2D S21(I, f), ida y vuelta, gráfico en vivo), `cargar_barrido_corriente`, `ajustar_barrido_corriente`, `saltos` (detección de saltos de fr con flujo), `ajustar_runs` |
| `fit_resonancia.py` | Análisis | `ajustar_notch` (modelo de Probst con asimetría φ y retardo de cable: fr, Ql, Qi, Qc, κ con errores), `cargar_traza`, `numero_fotones` |
| `fotones_resonador.py` | Análisis | Número de fotones a partir de un run, con ajuste lineal o no lineal (Kerr, Swenson et al. 2013) para potencias altas |

## Notebooks (`notebooks/`)

| Notebook | Qué hace |
|---|---|
| `readout.ipynb` | Notebook base de las primeras mediciones: single-tone, barridos en flujo y two-tone. También incluye celdas de mediciones del transmon, cuyas bases **no** están en este repo. |
| `medir_resonadores.ipynb` | Medición de resonadores con el VNA en modo single sweep: configuración, medición, ajuste, barrido en potencia con tabla CSV y barrido en potencia para cada corriente. |
| `punch_out.ipynb` | Barrido en corriente de la bobina (ida y vuelta) para ver la respuesta del resonador al flujo, verificación en un punto y barrido en potencia a corriente fija para buscar punch-out. |

## Datos (`datos/`)

Todas las bases son SQLite de QCoDeS. Se abren con `initialise_or_create_database_at(ruta)` y `load_by_id(run_id)`. Los parámetros de cada run (potencia, atenuación, resonador, corriente, etc.) están en `ds.metadata`.

| Base | Fecha | Configuración | Contenido |
|---|---|---|---|
| `fluxonium_database/resonadores` | 12–13/10/2026 | resonadores λ/4 | 20 runs: trazas del VNA a distintas potencias |
| `fluxonium_database/resonancia_flux` | 12–13/10/2026 | resonadores λ/4 | 136 runs: barridos en flujo de la resonancia (S21 vs f e I) y two-tone vs corriente. Muestran saltos tipo diente de sierra con período ~60–65 µA. |
| `fluxonium_v2` | 06/11/2026 | λ/2 (corto a tierra cortado) | 26 runs: single-tone vs corriente cerca de 5.41 GHz y 5.99 GHz |
| `fluxonium_2/resonadores_vna.db` | 07/11/2026 | λ/2 | 374 runs: res_q1 (5.41 GHz), res_q2 (6.07 GHz) y res_q4 (5.99 GHz); barridos en potencia (`barrido`) y en potencia × corriente (`barrido_corriente`) |
| `fluxonium_2/punch_out.db` | 07–08/11/2026 | λ/2 | Barridos en corriente ida/vuelta de res_q4 (±200 µA, paso 2 µA). Runs 1–2 (07/11): −10 dBm con −112 dB de atenuación de entrada (⟨n⟩ ≈ 1.7), analizados en el informe. Runs 3–4 (08/11): −10 dBm con −81 dB (⟨n⟩ ≈ 2×10³). |
| `fluxonium_2/barrido_*.csv` | 07/11/2026 | λ/2 | Tablas de ajustes de los barridos en potencia (fr, Ql, Qi, Qc, ⟨n⟩ vs P) |

Las bases se subieron como copias consistentes, sin los archivos `-wal` y `-shm`.

## Informes (`informes/`)

- `informe_barrido_flujo_res_q4.pdf`: análisis del barrido en corriente de res_q4 (λ/2). Muestra saltos periódicos de fr de ~200 kHz con período 52.3 ± 0.8 µA y una histéresis de ~0.47 períodos entre ida y vuelta. Incluye la interpretación (un lazo superconductor con estados de fluxoide metaestables) y los pasos para determinar si hay un qubit acoplado.
- `fluxonium_fit.png`: mapa |S21| y fase vs corriente (0–70 µA) de res_q4 con la fr ajustada.

## Historia del chip

1. **Octubre 2026.** Resonadores λ/4. Barridos en flujo con saltos periódicos de ~0.7–1 MHz (`fluxonium_database`).
2. **Noviembre 2026.** Se cortó el corto a tierra y los resonadores pasaron a λ/2 (frecuencias más altas). res_q4: fr ≈ 5.9924 GHz, Ql ≈ 755, Qi ≈ 1020, Qc ≈ 2840, κ/2π ≈ 8 MHz.
