# Dos tonos y un tono con flujo fino: búsqueda del qubit (10–12/11/2026)

Datos: `datos/two_tone_flux2/two_tone_flux2.db`. Código: `codigo/dos_tonos.py`. Notebook: `notebooks/two_tone_flux2.ipynb`.

**Resultado:** no se encontró ninguna señal del qubit, ni en dos tonos ni en un tono.

## Período de flujo esperado

Con la misma bobina (4.5 kΩ), un lazo de ~1000 µm² tenía un período de 1.85 µA. Los lazos del fluxonium son de ~2700 µm², así que el período esperado es **≈ 0.685 µA**.

Los períodos de los saltos de fr vistos en los barridos de corriente (~24 y ~52 µA) corresponden a áreas de ~35–80 µm². No son el lazo del fluxonium.

## Dos tonos gruesos (runs 1–230, res_q4)

- **Configuración:** drive por la antena entre 0.1 y 13 GHz, entre −20 y +10 dBm a la salida del SG con 20 dB hasta la antena. Corrientes alrededor de −139 µA, −125 µA, −21 µA y 0 µA.
- **Lo que se ve:** caídas de |S21| sin cambio de fase cerca de 5.0, 5.41, 6.07 y 6.4 GHz, con el drive a +10 dBm (`vivo_run6_4-8GHz.png`). No dependen de la corriente. Se atribuyen a compresión del HEMT por el drive que llega a la línea de salida, no al qubit.
- **Runs 50–177** (paso de 0.1 µA): en los parámetros del resonador no hay ningún período entre 0.3 y 2 µA (`periodo_resonador.png`). Plegando los espectros de dos tonos con períodos de 0.40 a 1.20 µA tampoco aparece nada (`plegado.png`).
- **Limitación:** un paso de 0.1 µA es ~0.15 Φ0. Es grueso pero alcanza para seguir una línea que se mueva con el flujo.

## Un tono con paso fino (runs 233–238)

**Medición:** f CW fija sobre el flanco del resonador, barrido de 0 a −2 µA en pasos de 5 nA (~0.007 Φ0), VNA a −20 dBm con IFBW 1 Hz.

**Análisis:** cada punto se convierte a un corrimiento equivalente δfr con la derivada del ajuste notch de la traza de búsqueda. Después se hace FFT, periodograma de seno con nulo AR(1) y una prueba a período fijo. Por último se comparan runs repetidos del mismo resonador.

| Run | Resonador | Puntos CW | Resultado |
|---|---|---|---|
| 233 | res_q1 | 1 | pico marginal en P = 0.514 µA (p ≈ 0.08 global); no se repitió |
| 234 | res_q4 | 1 | nada |
| 235 | res_q4 | 5 | 4.9 ± 1.2 kHz en P = 0.549 µA; **no se reprodujo en el 236** |
| 236 | res_q4 | 10 | deriva lenta de ±20 kHz; en P = 0.549 µA da 0.2–3.4 kHz según la línea de base |
| 237 | res_q2 | 1 | ruido ~42 kHz por punto; bulto en 1.2–1.6 µA, del tamaño del barrido (deriva) |
| 238 | res_q2 | 3 | ruido ~22 kHz por punto; nada; las curvas lentas son **anticorrelacionadas** con el 237 (−0.76) |

**Límites:** cualquier modulación de fr tiene que ser menor que **~3–4 kHz en res_q4** (períodos de 0.4 a 0.8 µA) y **~4–5 kHz en res_q2** (períodos de 0.5 a 1.6 µA, promedio pesado de 237 y 238).

**Figuras:** `fft_23x.png`, `reproduce_q4.png` (234–236) y `reproduce_q2.png` (237 vs 238).

**Lección:** la deriva lenta de fr, de decenas de kHz en decenas de minutos, imita períodos de ~1 µA en un barrido lento. Conviene repetir barridos rápidos iguales y promediarlos (la deriva cambia en cada repetición) antes que sumar más puntos CW por corriente.

## Dos tonos con paso fino (runs 239–619, res_q4)

**Configuración:**
- Búsqueda del mínimo una sola vez (`buscar_cada=0`); la f CW queda fija en 5.990000 GHz, 2.1 MHz por debajo de fr = 5.99214 GHz (κ = 6.1 MHz).
- Drive de 1 MHz a 4 GHz en pasos de 5 MHz, a 0 dBm a la salida del SG con 20 dB hasta la antena.
- 1 punto CW a IFBW 10 Hz, n_ref = 5.
- 381 corrientes de −0.005 a −1.905 µA en pasos de 5 nA; duró 14.5 h.

**Análisis:**
- La referencia de cada run es la mediana compleja de su propio barrido de drive.
- Se proyecta a δfr y a la componente ortogonal.
- Se resta el fondo común (la mediana de todos los runs en cada f de drive) y se agrupa de a 10 corrientes (50 nA).

**Lo que da:**
- **Ruido:** 116 kHz por punto en δfr y 36 kHz por celda agrupada, lo esperado para ruido blanco.
- **Mapas** (`mapa_todo.png`): ninguna de las ~30 000 celdas pasa 4.5σ. Una línea que corriera el resonador más de ~150 kHz en una celda se vería; siguiendo una curva a lo largo de varias celdas, más de ~80 kHz.
- **Fondo común** (`fondo_ref_todo.png`, arriba): 7 kHz por f de drive, igual al ruido. No hay líneas fijas mayores que ~25 kHz.
- **Efecto global del drive** (prendido contra apagado): −4.3 ± 2.6 kHz, compatible con cero.
- **Referencia sin drive vs I** (`fondo_ref_todo.png`, medio): el resonador deriva de −350 a +390 kHz en las 14.5 h. El pico del periodograma depende de la tendencia que se resta (`ref_periodo.png`), y las dos mitades del barrido dan amplitudes y fases distintas. Es deriva, no flujo.

**Conclusión:** a 0 dBm (−20 dBm en la antena), entre 1 MHz y 4 GHz, no hay ninguna transición que corra el resonador más de ~100 kHz. Las explicaciones posibles:
- χ es más chico que eso.
- Al qubit le llega poco drive.
- La línea es más angosta que el paso de 5 MHz.
- La transición cae fuera de 0–4 GHz.

Los runs 620–679 repiten el barrido con el drive a +5 dBm, de −0.005 a −0.30 µA; se cortó ahí para pasar a la configuración más sensible.

## Dos tonos sensible (runs 681–713, res_q4)

Datos: `datos/two_tone_flux2/two_tone_flux2_651-723.db`.

**Configuración:**
- Lectura a 0 dBm en el VNA (⟨n⟩ ≈ 17, sin efecto en un qubit tan desacoplado), con 10 puntos CW a 10 Hz por f de drive.
- f CW en fr: 0.23 MHz debajo, con `desplazamiento_cw` = +2.1 MHz desde el mínimo de |S21|.
- Drive a +10 dBm a la salida del SG, de 1 MHz a 4 GHz en pasos de 5 MHz.
- 33 corrientes de −0.005 a −0.805 µA en pasos de 25 nA, solo ida (8.5 h; la vuelta no se midió). Cubre más de un período si es de 0.5–0.7 µA.

**Ruido:** ~2 kHz por punto en δfr, unas 60 veces menos que el barrido 239–619. Promediando las 33 corrientes baja a 0.35 kHz por punto.

**Tres líneas fijas** (`ida_681-713.png`, `lineas_fijas_run681.png`). Las tres tienen fr↓ y pérdida↑, el patrón de absorción o calentamiento:

| Línea | Forma | δfr | Pérdida |
|---|---|---|---|
| 0.293 GHz | 2 puntos de la grilla (≤ 10 MHz) | −4 kHz | +9 kHz |
| 2.81 GHz | rectangular, 2.796–2.826 GHz, bordes abruptos (¿artefacto del generador?) | −5 kHz | +9 kHz |
| 3.276 GHz | 1 punto, con hombros asimétricos ±50 MHz | −40 kHz | +46 kHz; ya vista en los runs 1–230 |

Ninguna cambia con la corriente: posición (0.293: +0.1 ± 0.8 MHz/µA, `linea_0293.png`), intensidad y forma son constantes.

**Búsqueda sin supuestos** (`busqueda_sin_supuestos_681-713.png`). Con nulo de datos desordenados en todos los casos:
- La posición fina de las líneas fijas, usando el cociente entre puntos vecinos.
- Cada f de drive: tendencia, curvatura y seno de período 0.05–1.6 µA.
- Rectas de hasta ±6 GHz/µA.
- Trayectorias libres por programación dinámica, con saltos de 5–400 MHz entre corrientes.

Nada supera lo esperable para el ruido. La mejor trayectoria libre no se reproduce entre corrientes pares e impares.

**Prueba de corriente contra tiempo** (runs 714–723, `prueba_I_vs_tiempo_714-723.png`). Quedaban dos indicios débiles: los bordes de 2.81 GHz en δfr y una curvatura en 2.996 GHz = f_CW/2. Se midieron de nuevo en 5 corrientes, dos veces cada una y en orden alternado:
- **Bordes de 2.81 GHz:** +0.8 ± 1.7 kHz/µA, contra −3.7 kHz/µA en la ida, pero cambian con el tiempo. Era deriva.
- **2.996 GHz:** la forma de la ida aparece con amplitud −0.2 ± 0.6. No se repite.

**fr vs I:** una joroba de +79 kHz cerca de −0.5 µA que no se repite en el barrido 239–619. Con solo la ida no se separa de la deriva.

**Conclusión:** entre 0 y −0.8 µA, con drive de +10 dBm entre 1 MHz y 4 GHz, no hay nada que dependa del flujo por encima de ~3–4 kHz de corrimiento del resonador. Ninguna línea fija se corre más de ~1 MHz.
