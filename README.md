# Fluxonium: mediciones de resonadores de readout

Este repositorio reúne el código, los notebooks y los datos de las mediciones del chip **fluxonium** del grupo, hechas en octubre y noviembre de 2026. Sirven para dos cosas:

- **Medir.** El código controla el VNA y la fuente de la bobina, y se puede reutilizar para mediciones parecidas.
- **Postprocesar.** Las bases de QCoDeS tienen los datos crudos y toda la metadata: potencia, atenuación, corriente y estado del VNA.

## Estructura

```
fluxonium/
├── codigo/       módulos de Python (control de equipos, dos tonos y ajustes)
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
| Generador (two-tone) | Keysight N5173B (`GPIB1::19::INSTR`), usado en `readout.ipynb` y `two_tone_flux2.ipynb`; drive por la antena con 20 dB de atenuación |
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
| `dos_tonos.py` | VNA + generador N5173B + fuente | `medir_dos_tonos` (espectroscopía de dos tonos vs corriente: busca el mínimo del resonador con el marcador del VNA y lee en CW mientras barre el drive; `buscar_cada=0` busca solo al principio y deja la f CW fija; un run por corriente; ventana en vivo con mapas I vs f drive), `medir_un_tono` (barrido de corriente leyendo S21 en CW sin drive, un solo run), `cargar_dos_tonos`, `cargar_un_tono`, `buscar_dos_tonos`, `tiempo_estimado` |
| `fotones_resonador.py` | Análisis | Número de fotones a partir de un run, con ajuste lineal o no lineal (Kerr, Swenson et al. 2013) para potencias altas |

## Notebooks (`notebooks/`)

| Notebook | Qué hace |
|---|---|
| `readout.ipynb` | Notebook base de las primeras mediciones: single-tone, barridos en flujo y two-tone. También incluye celdas de mediciones del transmon, cuyas bases **no** están en este repo. |
| `medir_resonadores.ipynb` | Medición de resonadores con el VNA en modo single sweep: configuración, medición, ajuste, barrido en potencia con tabla CSV y barrido en potencia para cada corriente. |
| `punch_out.ipynb` | Barrido en corriente de la bobina (ida y vuelta) para ver la respuesta del resonador al flujo, verificación en un punto y barrido en potencia a corriente fija para buscar punch-out. |
| `two_tone_flux2.ipynb` | Dos tonos (drive por la antena) vs corriente de la bobina con `dos_tonos.medir_dos_tonos`, y un tono con paso fino de corriente con `medir_un_tono`. Usa `datos/two_tone_flux2/two_tone_flux2.db`; los runs 651–723 están en `two_tone_flux2_651-723.db`. Incluye una celda para graficar Δ|S21| y Δfase de cualquier medición por su tag. |

## Datos (`datos/`)

Todas las bases son SQLite de QCoDeS. Se abren con `initialise_or_create_database_at(ruta)` y `load_by_id(run_id)`. Los parámetros de cada run (potencia, atenuación, resonador, corriente, etc.) están en `ds.metadata`.

| Base | Fecha | Configuración | Contenido |
|---|---|---|---|
| `fluxonium_database/resonadores` | 12–13/10/2026 | resonadores λ/4 | 20 runs: trazas del VNA a distintas potencias |
| `fluxonium_database/resonancia_flux` | 12–13/10/2026 | resonadores λ/4 | 136 runs: barridos en flujo de la resonancia (S21 vs f e I) y two-tone vs corriente. Muestran saltos tipo diente de sierra con período ~60–65 µA. |
| `fluxonium_v2` | 06/11/2026 | λ/2 (corto a tierra cortado) | 26 runs: single-tone vs corriente cerca de 5.41 GHz y 5.99 GHz |
| `fluxonium_2/resonadores_vna.db` | 07/11/2026 | λ/2 | 374 runs: res_q1 (5.41 GHz), res_q2 (6.07 GHz) y res_q4 (5.99 GHz); barridos en potencia (`barrido`) y en potencia × corriente (`barrido_corriente`) |
| `fluxonium_2/punch_out.db` | 07–10/11/2026 | λ/2 | Barridos en corriente ida/vuelta (`barrido_I`), VNA a −10 dBm salvo el run 30 (0 dBm). Runs 1–10: res_q4 (±200 µA y zooms; runs 1–2 con −112 dB de entrada, ⟨n⟩ ≈ 1.7, analizados en el informe; 3–4 con −81 dB; 5–10 con `att_vna` = −50). Runs 11–14: res_q2. Runs 15–33: res_q1 (el 32: −200 a −100 µA en 1001 pasos; el 33: ±500 µA). Runs 34–35: res_q4 (−200 a 0 µA). Se sacaron los runs 27, 29 y 31 (vacíos o cortados al empezar). |
| `two_tone_flux2/two_tone_flux2.db` | 10–12/11/2026 | λ/2 | Dos tonos y un tono, ver `informes/dos_tonos/README.md`. Runs 1–230: dos tonos res_q4 (un run por corriente, metadata `dos_tonos` = etiqueta del barrido). Runs 233–238: un tono (`un_tono`), paso de 5 nA de 0 a −2 µA (233 res_q1, 234–236 res_q4, 237–238 res_q2). Runs 239–619: dos tonos res_q4 con f CW fija, 1 MHz–4 GHz, 0 dBm, de −0.005 a −1.905 µA en pasos de 5 nA. Runs 620–650: el mismo barrido a +5 dBm (sigue en la base de abajo). Se sacaron los runs cortados al principio (3, 7, 10, 182, 231, 232) y el que estaba midiéndose; los runs 4, 49 y 619 quedaron con el 74–96 % del barrido de drive. |
| `two_tone_flux2/two_tone_flux2_651-723.db` | 12–13/11/2026 | λ/2 | Solo los runs nuevos, con los mismos `run_id` (se abre con `initialise_or_create_database_at('two_tone_flux2/two_tone_flux2_651-723.db')`). Runs 651–679: resto del barrido a +5 dBm (hasta −0.30 µA). Runs 681–713: dos tonos sensible, VNA 0 dBm, 10 pts CW, f CW en fr, drive +10 dBm, de −0.005 a −0.805 µA en pasos de 25 nA (solo ida). Runs 714–723: prueba de corriente contra tiempo, 2.75–2.88 y 2.95–3.05 GHz, 5 corrientes ×2 en orden alternado, 30 pts CW. Se sacó el run 680 (cortado al cambiar de configuración). |
| `fluxonium_2/barrido_*.csv` | 07/11/2026 | λ/2 | Tablas de ajustes de los barridos en potencia (fr, Ql, Qi, Qc, ⟨n⟩ vs P) |

Las bases se subieron como copias consistentes, sin los archivos `-wal` y `-shm`. En las bases nuevas se borraron los runs incompletos conservando el `run_id` del resto, así que puede haber huecos en la numeración.

## Informes (`informes/`)

- `informe_barrido_flujo_res_q4.pdf`: análisis del barrido en corriente de res_q4 (λ/2). Muestra saltos periódicos de fr de ~200 kHz con período 52.3 ± 0.8 µA y una histéresis de ~0.47 períodos entre ida y vuelta. Incluye la interpretación (un lazo superconductor con estados de fluxoide metaestables) y los pasos para determinar si hay un qubit acoplado.
- `fluxonium_fit.png`: mapa |S21| y fase vs corriente (0–70 µA) de res_q4 con la fr ajustada.
- `informe_barrido_flujo_res_q4_potencia.pdf`: el mismo barrido de res_q4 con −81 dB de entrada (⟨n⟩ ≈ 2×10³, `punch_out.db` runs 3–4). El diente de sierra no cambia con la potencia (saltos de 209 ± 25 kHz, mismo período e histéresis). Todas las ramas son la misma parábola desplazada en I, con período 50.2 µA y ventana de existencia de cada estado de fluxoide de ±67 µA.
- `informe_final/informe_fluxonium.tex`: informe en LaTeX de todas las mediciones (resonadores, ramas de fluxoide, espectroscopía de uno y dos tonos), las hipótesis físicas y por qué no se puede confirmar un qubit funcionando. Las figuras están en `informe_final/figuras/`; se compila con `pdflatex informe_fluxonium.tex` (dos pasadas, por el índice).
- `dos_tonos/`: búsqueda del qubit con dos tonos y un tono con paso fino de corriente. Resumen en `dos_tonos/README.md`, con las figuras de cada análisis. No se encontró señal del qubit. La modulación de fr con el flujo es menor que ~3–5 kHz. En el barrido sensible (runs 681–713) no hay nada que dependa del flujo por encima de ~3–4 kHz entre 1 MHz y 4 GHz; solo aparecen tres líneas fijas (0.293, 2.81 y 3.276 GHz).

## Historia del chip

1. **Octubre 2026.** Resonadores λ/4. Barridos en flujo con saltos periódicos de ~0.7–1 MHz (`fluxonium_database`).
2. **Noviembre 2026.** Se cortó el corto a tierra y los resonadores pasaron a λ/2 (frecuencias más altas). res_q4: fr ≈ 5.9924 GHz, Ql ≈ 755, Qi ≈ 1020, Qc ≈ 2840, κ/2π ≈ 8 MHz.
3. **10–13/11/2026.** Dos tonos con drive por la antena y un tono con paso fino de corriente (`two_tone_flux2.db`), y un barrido sensible de dos tonos con su prueba de control (`two_tone_flux2_651-723.db`). Sin señal del qubit.
