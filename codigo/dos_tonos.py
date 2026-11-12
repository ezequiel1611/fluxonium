r"""
Espectroscopía de dos tonos con flujo de bobina (lo usa two_tone_flux2.ipynb).

Instrumentos:
    pna    : Keysight PNA (driver QCoDeS KeysightN5245A) -> lectura del resonador
    sg     : Keysight N5173B (driver QCoDeS)             -> drive por la antena del qubit
    fuente : barrido_flujo.FuenteCorriente (Keithley)    -> corriente de la bobina (flujo)

Para cada corriente de la bobina:
    1) rampa de la fuente hasta la corriente y espera
    2) drive APAGADO; el VNA barre en frecuencia alrededor del resonador y el marcador busca
       el mínimo de |S21| (FUNC MIN con tracking, sobre la traza S21 en formato MLOG).
       Con CALC:MEAS:MARK:SET CWFReq el VNA copia la frecuencia del marcador a su frecuencia CW.
    3) VNA en CW en esa frecuencia; referencia con el drive apagado
    4) drive encendido: para cada potencia y cada frecuencia del SG, N puntos CW -> promedio complejo

El mínimo se busca solo en la primera corriente (buscar_cada=0, por defecto): después la f CW
queda fija y en cada corriente solo se toma la referencia sin drive y se barre el drive. Con
buscar_cada=n se vuelve a buscar cada n corrientes (1 = en todas). Nunca se busca al cambiar la
frecuencia o la potencia del drive.

Cada corriente es un run propio de QCoDeS (su propio run_id), con:
    p_drive [dBm], f_drive [Hz]          -> setpoints
    s21_re, s21_im, s21_mag [dB], s21_phase [deg], s21_err (error del promedio complejo)
    f_busqueda [Hz] y busqueda_re/_im/_mag  -> la traza con la que se buscó el mínimo
    metadata: corriente, f del mínimo, f CW, referencia sin drive, config del VNA, SG, etc.

Un tono (medir_un_tono): sin drive, solo barrido de la corriente de la bobina leyendo S21 en CW.
Todo el barrido es UN run de QCoDeS con setpoint I_bobina. La f CW puede quedar fija (se busca el
mínimo una vez al principio, o se da f_cw) o seguir al mínimo cada `buscar_cada` corrientes.
Se carga con cargar_un_tono(run_id).
"""

import json
import time

import numpy as np

from vna_single import _info_trazas, _leer_binario, disparar_single, eje_x, leer_estado


# =============================================================================
# 1) VNA: búsqueda del mínimo con el marcador y lectura en CW
# =============================================================================
def traza_s21_mlog(pna):
    """Número de medición de la traza S21 en formato MLOG (la que usa el marcador)."""
    for t in _info_trazas(pna):
        if t["sparam"] == "S21" and t["formato"] == "MLOG":
            return t["num"]
    raise ValueError("no hay una traza S21 en formato MLOG en el VNA: el marcador busca el mínimo "
                     "de |S21| en dB, agregá esa traza en el panel (o cambiá el formato de una S21)")


def configurar_busqueda(pna, centro, span, puntos, ifbw, potencia=None):
    """VNA en barrido LIN alrededor del resonador (promedios apagados)."""
    pna.sweep_type("LIN")
    # directo por SCPI: el validador del driver no deja spans < 10 MHz
    pna.write(f"SENS:FREQ:CENT {float(centro):.3f}")
    pna.write(f"SENS:FREQ:SPAN {float(span):.3f}")
    pna.points(int(puntos))
    pna.if_bandwidth(float(ifbw))
    if potencia is not None:
        pna.power(float(potencia))
    pna.averages_enabled(False)
    pna.ask("*OPC?")


def buscar_minimo(pna, num, centro, span, puntos, ifbw, potencia=None, marcador=1,
                  desplazamiento=0.0):
    """Un barrido LIN y el marcador busca el mínimo de |S21|; la f CW del VNA queda en el mínimo.

    num: número de medición de la traza S21/MLOG (traza_s21_mlog).
    desplazamiento [Hz]: la f CW queda en f_min + desplazamiento (0 = en el mínimo).
    Deja el VNA en HOLD, todavía en LIN. Devuelve dict(f_min, s21_min_db, f_cw_objetivo, f, z, estado).
    """
    # la configuración se envía 2 veces: al volver de CW el VNA no siempre toma todo a la primera
    for _ in range(2):
        configurar_busqueda(pna, centro, span, puntos, ifbw, potencia)
    m = f"CALC:MEAS{num}:MARK{marcador}"
    pna.write(f"{m}:STAT ON")
    pna.write(f"{m}:FUNC:SEL MIN")
    pna.write(f"{m}:FUNC:TRAC ON")              # busca el mínimo en cada barrido
    est = leer_estado(pna)
    disparar_single(pna, est)
    pna.write(f"{m}:FUNC:EXEC MIN")             # por si el tracking todavía no actualizó
    pna.ask("*OPC?")
    f_min = float(pna.ask(f"{m}:X?"))
    s21_min = float(pna.ask(f"{m}:Y?").split(",")[0])
    # traza con la que se buscó el mínimo (compleja, sin electrical delay)
    activa = pna.active_trace()
    pna.active_trace(num)
    d = _leer_binario(pna, "CALC:DATA? SDATA")
    pna.active_trace(activa)
    z = d[0::2] + 1j * d[1::2]
    f = eje_x(pna, est)
    # el tracking se apaga antes de pasar a CW (en CW buscaría el mínimo en el eje de tiempo)
    pna.write(f"{m}:FUNC:TRAC OFF")
    f_obj = f_min + float(desplazamiento)
    if desplazamiento:
        pna.write(f"{m}:X {f_obj:.3f}")
    pna.write(f"{m}:SET CWFR")                  # marcador -> frecuencia CW (no cambia el sweep type)
    pna.ask("*OPC?")
    return dict(f_min=f_min, s21_min_db=s21_min, f_cw_objetivo=f_obj, f=f, z=z, estado=est)


def configurar_cw(pna, num, puntos, ifbw, f_esperada, tol=2.0):
    """VNA en CW (la f CW ya la puso el marcador). Verifica que coincida con el mínimo.

    Deja seleccionada la traza `num` para leer SDATA. Devuelve (f_cw, estado).
    """
    pna.sweep_type("CW")
    pna.points(int(puntos))
    pna.if_bandwidth(float(ifbw))
    pna.averages_enabled(False)
    f_cw = float(pna.cw())
    if abs(f_cw - f_esperada) > tol:
        print(f"  aviso: la f CW del VNA ({f_cw/1e9:.6f} GHz) no es la del marcador "
              f"({f_esperada/1e9:.6f} GHz); la seteo directamente")
        pna.cw(float(round(f_esperada)))
        f_cw = float(pna.cw())
    pna.active_trace(num)
    pna.ask("*OPC?")
    return f_cw, leer_estado(pna)


def leer_cw(pna, est):
    """Un barrido CW (N puntos) de la traza seleccionada -> (promedio complejo, error del promedio)."""
    disparar_single(pna, est)
    d = _leer_binario(pna, "CALC:DATA? SDATA")
    z = d[0::2] + 1j * d[1::2]
    err = np.std(z) / np.sqrt(len(z)) if len(z) > 1 else np.nan
    return complex(np.mean(z)), float(err)


# =============================================================================
# 2) Figura en vivo (matplotlib tk)
# =============================================================================
class _Vivo:
    """Ventana que se actualiza mientras se mide.

    Arriba: traza de búsqueda con el mínimo / f CW y referencia sin drive vs corriente.
    Medio : Δ|S21| y Δfase vs f del drive (respecto de la referencia sin drive) de la corriente
            y potencia que se están midiendo.
    Abajo : mapas de Δ|S21| y Δfase. Eje x = corriente de la bobina si hay más de una
            (para la potencia p_mapa), si no la potencia del drive; eje y = f del drive.
    """

    def __init__(self, corrientes, p_drive, f_drive, p_mapa=0, cada=0.3):
        import matplotlib.pyplot as plt
        self.plt = plt
        self.cada = cada
        self.t_ult = 0.0
        self.corrientes = np.asarray(corrientes, float)
        self.p_drive = np.asarray(p_drive, float)
        self.f_drive = np.asarray(f_drive, float)
        self.p_mapa = p_mapa
        self.por_corriente = len(self.corrientes) > 1
        ny = len(self.corrientes) if self.por_corriente else len(self.p_drive)
        self.yv = self.corrientes * 1e6 if self.por_corriente else self.p_drive
        self.orden = np.argsort(self.yv)
        self.fila = np.empty(ny, int)
        self.fila[self.orden] = np.arange(ny)            # índice medido -> fila del mapa
        self.Mm = np.full((ny, len(self.f_drive)), np.nan)
        self.Mp = np.full((ny, len(self.f_drive)), np.nan)
        self.I_hist, self.fcw_hist = [], []
        self.Iref_hist, self.ref_hist = [], []

        plt.ion()
        fig = plt.figure(figsize=(13, 9.5))
        gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 1.3])
        self.ax_b = fig.add_subplot(gs[0, 0])
        self.ax_c = fig.add_subplot(gs[0, 1])
        self.ax_m = fig.add_subplot(gs[1, 0])
        self.ax_p = fig.add_subplot(gs[1, 1], sharex=self.ax_m)
        self.ax_Mm = fig.add_subplot(gs[2, 0])
        self.ax_Mp = fig.add_subplot(gs[2, 1], sharex=self.ax_Mm, sharey=self.ax_Mm)
        self.fig = fig

        self.ln_b, = self.ax_b.plot([], [], lw=1)
        self.pt_b, = self.ax_b.plot([], [], "rv", ms=7, label="mínimo (marcador)")
        self.vl_b = self.ax_b.axvline(np.nan, color="r", ls=":", lw=1, label="f CW")
        self.ax_b.set_xlabel("f VNA [GHz]"); self.ax_b.set_ylabel("|S21| [dB]")
        self.ax_b.legend(fontsize=8, loc="lower right")
        self.ln_c, = self.ax_c.plot([], [], "o-", ms=4, label="f CW (búsquedas)")
        self.ax_c.set_xlabel("I bobina [uA]"); self.ax_c.set_ylabel("f CW − primera [kHz]", color="C0")
        self.ax_c2 = self.ax_c.twinx()
        self.ln_r, = self.ax_c2.plot([], [], ".-", ms=3, lw=0.8, color="C1")
        self.ax_c2.set_ylabel("ref sin drive |S21| [dB]", color="C1")
        self.ax_c.set_title("f CW y referencia sin drive vs corriente", fontsize=10)

        fg = self.f_drive / 1e9
        self.ln_m, = self.ax_m.plot([], [], lw=1)
        self.ln_p, = self.ax_p.plot([], [], lw=1, color="C3")
        self.ax_m.set_ylabel("Δ|S21| [dB]"); self.ax_p.set_ylabel("Δfase [°]")
        for a in (self.ax_m, self.ax_p):
            a.axhline(0, color="0.6", lw=0.8)
            a.set_xlim(fg.min(), fg.max())
            a.grid(alpha=0.3)

        # mapas: x = corriente (o potencia), y = f del drive -> se dibuja M.T
        x0, x1 = self.yv[self.orden[0]], self.yv[self.orden[-1]]
        if x1 == x0:
            x0, x1 = x0 - 0.5, x0 + 0.5
        ext = [x0, x1, fg.min(), fg.max()]
        cm = plt.get_cmap("RdBu_r").copy(); cm.set_bad("0.85")
        self.im_m = self.ax_Mm.imshow(self.Mm.T, aspect="auto", origin="lower", extent=ext, cmap=cm,
                                      interpolation="nearest")
        self.im_p = self.ax_Mp.imshow(self.Mp.T, aspect="auto", origin="lower", extent=ext, cmap=cm,
                                      interpolation="nearest")
        fig.colorbar(self.im_m, ax=self.ax_Mm, label="Δ|S21| [dB]")
        fig.colorbar(self.im_p, ax=self.ax_Mp, label="Δfase [°]")
        xlab = "I bobina [uA]" if self.por_corriente else "P drive [dBm]"
        for a in (self.ax_Mm, self.ax_Mp):
            a.set_xlabel(xlab); a.set_ylabel("f drive [GHz]")
        if self.por_corriente and len(self.p_drive) > 1:
            self.ax_Mm.set_title(f"mapa para P drive = {self.p_drive[p_mapa]:g} dBm", fontsize=10)
        self.titulo = fig.suptitle("")
        fig.tight_layout()
        plt.show(block=False)
        self._dibujar()

    def _dibujar(self):
        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()
        self.t_ult = time.time()

    def busqueda(self, I, b, f_cw):
        """Nueva corriente: traza de búsqueda, mínimo y f CW."""
        mag = 20 * np.log10(np.abs(b["z"]))
        self.ln_b.set_data(b["f"] / 1e9, mag)
        self.pt_b.set_data([b["f_min"] / 1e9], [b["s21_min_db"]])
        self.vl_b.set_xdata([f_cw / 1e9, f_cw / 1e9])
        self.ax_b.relim(); self.ax_b.autoscale_view()
        self.ax_b.set_title(f"búsqueda en I = {I*1e6:.2f} uA | f CW = {f_cw/1e9:.6f} GHz", fontsize=10)
        self.I_hist.append(I * 1e6); self.fcw_hist.append(f_cw)
        self.ln_c.set_data(self.I_hist, (np.array(self.fcw_hist) - self.fcw_hist[0]) / 1e3)
        self.ax_c.relim(); self.ax_c.autoscale_view()
        self._dibujar()

    def referencia(self, I, z_ref):
        """Referencia sin drive de cada corriente (con f CW fija: el resonador visto a un tono)."""
        self.Iref_hist.append(I * 1e6); self.ref_hist.append(20 * np.log10(abs(z_ref)))
        self.ln_r.set_data(self.Iref_hist, self.ref_hist)
        self.ax_c2.relim(); self.ax_c2.autoscale_view(scalex=False)
        xs = self.Iref_hist + self.I_hist
        x0, x1 = min(xs), max(xs)
        d = 0.03 * (x1 - x0) or 1e-3
        self.ax_c.set_xlim(x0 - d, x1 + d)

    def punto(self, kI, kP, fila_m, fila_p, titulo, forzar=False):
        """Actualiza las curvas de la corriente/potencia en curso y el mapa (con límite de tasa)."""
        fg = self.f_drive / 1e9
        self.ln_m.set_data(fg, fila_m)
        self.ln_p.set_data(fg, fila_p)
        k = kI if self.por_corriente else kP
        if (not self.por_corriente) or kP == self.p_mapa:
            self.Mm[self.fila[k]] = fila_m
            self.Mp[self.fila[k]] = fila_p
        if not forzar and time.time() - self.t_ult < self.cada:
            return
        for a in (self.ax_m, self.ax_p):
            a.relim(); a.autoscale_view(scalex=False)
        for im, M in ((self.im_m, self.Mm), (self.im_p, self.Mp)):
            im.set_data(M.T)
            if np.isfinite(M).any():
                v = np.nanpercentile(np.abs(M), 99) or 1e-3
                im.set_clim(-v, v)
        self.titulo.set_text(titulo)
        self._dibujar()


# =============================================================================
# 3) Medición
# =============================================================================
def _medicion(nombre):
    from qcodes.dataset import Measurement
    m = Measurement(name=nombre)
    m.write_period = 5
    m.register_custom_parameter("p_drive", label="potencia del drive (SG)", unit="dBm")
    m.register_custom_parameter("f_drive", label="frecuencia del drive (SG)", unit="Hz")
    sp = ("p_drive", "f_drive")
    m.register_custom_parameter("s21_re", label="Re S21 (CW)", unit="", setpoints=sp)
    m.register_custom_parameter("s21_im", label="Im S21 (CW)", unit="", setpoints=sp)
    m.register_custom_parameter("s21_mag", label="|S21| (CW)", unit="dB", setpoints=sp)
    m.register_custom_parameter("s21_phase", label="fase S21 (CW)", unit="deg", setpoints=sp)
    m.register_custom_parameter("s21_err", label="error del promedio de S21", unit="", setpoints=sp)
    m.register_custom_parameter("f_busqueda", label="f VNA (búsqueda del mínimo)", unit="Hz")
    sb = ("f_busqueda",)
    m.register_custom_parameter("busqueda_re", label="Re S21 búsqueda", unit="", setpoints=sb)
    m.register_custom_parameter("busqueda_im", label="Im S21 búsqueda", unit="", setpoints=sb)
    m.register_custom_parameter("busqueda_mag", label="|S21| búsqueda", unit="dB", setpoints=sb)
    return m


def _n_busquedas(n_corrientes, buscar_cada=0):
    """Cuántas veces se busca el mínimo en un barrido de n_corrientes."""
    if buscar_cada and buscar_cada > 0:
        return int(np.ceil(n_corrientes / int(buscar_cada)))
    return 1 if n_corrientes else 0


def tiempo_estimado(n_corrientes, n_p, n_f, puntos_cw, ifbw_cw, espera_drive,
                    puntos_busqueda, ifbw_busqueda, espera_I, n_ref=5, overhead=0.05, buscar_cada=0):
    """Duración aproximada [s] (sin contar la rampa de la bobina).

    buscar_cada: como en medir_dos_tonos (0 = una sola búsqueda al principio; 1 = en cada corriente).
    """
    t_pt = puntos_cw / ifbw_cw + espera_drive + overhead
    t_I = espera_I + n_ref * (puntos_cw / ifbw_cw + overhead)
    t_b = puntos_busqueda / ifbw_busqueda + 1.0
    return (n_corrientes * (t_I + n_p * n_f * t_pt)
            + _n_busquedas(n_corrientes, buscar_cada) * t_b), t_pt


def medir_dos_tonos(pna, sg, fuente, corrientes, f_drive, p_drive,
                    centro_busqueda, span_busqueda, puntos_busqueda=201, ifbw_busqueda=100,
                    p_vna=-10, puntos_cw=10, ifbw_cw=100, desplazamiento_cw=0.0,
                    buscar_cada=0, seguir_minimo=True, marcador=1,
                    espera_drive=0.01, excluir_cw=20e6, n_ref=5,
                    paso_rampa=1e-6, espera_rampa=0.05, espera_I=1.0, preparar_desde=None,
                    apagar_drive_al_final=True, volver_a=None,
                    vivo=True, p_mapa=0, refresco=0.5,
                    nombre="dos_tonos", tag=None, **extra_metadata):
    """Espectroscopía de dos tonos: busca el mínimo del resonador (marcador del VNA -> f CW) en la
    primera corriente y, para cada corriente de la bobina, toma la referencia sin drive en CW y
    barre el drive del SG en potencia y frecuencia.

    corrientes [A]  : corrientes de la bobina, en el orden en que se miden.
    f_drive [Hz]    : frecuencias del SG (en el orden en que se barren).
    p_drive [dBm]   : potencias del SG (una o varias; lazo externo a f_drive).
    centro_busqueda, span_busqueda [Hz], puntos_busqueda, ifbw_busqueda [Hz]: barrido LIN con el
                      que el marcador busca el mínimo de |S21| (con el drive apagado).
    p_vna [dBm]     : potencia del VNA (la misma para la búsqueda y para la lectura CW).
    puntos_cw, ifbw_cw: puntos CW que se promedian en cada frecuencia del drive y su IFBW.
    desplazamiento_cw [Hz]: f CW = f del mínimo + desplazamiento (0 = en el mínimo).
    buscar_cada     : 0 = busca el mínimo solo en la primera corriente y la f CW queda fija (el VNA
                      se queda en CW todo el barrido); n = vuelve a buscar cada n corrientes
                      (1 = en todas, como antes).
    seguir_minimo   : True = la próxima búsqueda se centra en el mínimo anterior.
    excluir_cw [Hz] : no mide frecuencias del drive a menos de esto de la f CW (None = mide todo).
    n_ref           : barridos CW promediados para la referencia sin drive en cada corriente.
    preparar_desde [A]: si no es None, antes de la primera corriente rampa hasta acá (histéresis).
    volver_a [A]    : corriente final de la bobina (None = queda en la última).
    vivo            : ventana de matplotlib (usar %matplotlib tk) que se actualiza mientras mide.
    p_mapa          : índice de p_drive que se muestra en el mapa cuando hay varias corrientes.
    extra_metadata  : se guarda en cada run (resonador, att_vna, att_drive, notas, ...).

    Devuelve dict(tag, runs=[{run_id, current_uA, f_min, f_cw}]).
    """
    corrientes = np.atleast_1d(np.asarray(corrientes, float))
    f_drive = np.atleast_1d(np.asarray(f_drive, float))
    p_drive = np.atleast_1d(np.asarray(p_drive, float))
    if np.any(np.abs(corrientes) > fuente.I_max):
        raise ValueError(f"hay corrientes por encima de I_max = {fuente.I_max*1e6:.1f} uA")
    for p in p_drive:                          # valida con los límites del driver antes de empezar
        sg.power.validate(float(p))
    for f in (f_drive.min(), f_drive.max()):
        sg.frequency.validate(float(f))

    tag = tag or time.strftime("%Y%m%d_%H%M%S")
    num = traza_s21_mlog(pna)
    buscar_cada = int(buscar_cada or 0)
    t_tot, t_pt = tiempo_estimado(len(corrientes), len(p_drive), len(f_drive), puntos_cw, ifbw_cw,
                                  espera_drive, puntos_busqueda, ifbw_busqueda, espera_I, n_ref,
                                  buscar_cada=buscar_cada)
    nb = _n_busquedas(len(corrientes), buscar_cada)
    print(f"dos tonos {tag}: {len(corrientes)} corriente(s) x {len(p_drive)} potencia(s) x "
          f"{len(f_drive)} frecuencias del drive, {nb} búsqueda(s) del mínimo | "
          f"~{t_pt*1e3:.0f} ms por punto -> ~{t_tot/60:.1f} min (sin la rampa de la bobina)")

    def toca_buscar(k):
        return k == 0 or (buscar_cada > 0 and k % buscar_cada == 0)

    md_comun = dict(
        dos_tonos=tag, n_corrientes=len(corrientes),
        corrientes_uA=json.dumps([round(c * 1e6, 4) for c in corrientes]),
        f_drive_start=float(f_drive[0]), f_drive_stop=float(f_drive[-1]), n_f_drive=len(f_drive),
        f_drive_step=float(np.median(np.diff(f_drive))) if len(f_drive) > 1 else 0.0,
        p_drive=json.dumps([float(p) for p in p_drive]),
        centro_busqueda=float(centro_busqueda), span_busqueda=float(span_busqueda),
        puntos_busqueda=int(puntos_busqueda), ifbw_busqueda=float(ifbw_busqueda),
        pow_vna=float(p_vna), N_points_cw=int(puntos_cw), bw_vna=float(ifbw_cw),
        desplazamiento_cw=float(desplazamiento_cw), seguir_minimo=int(bool(seguir_minimo)),
        buscar_cada=buscar_cada, marcador=int(marcador), traza_marcador=int(num),
        espera_drive=float(espera_drive), excluir_cw=float(excluir_cw or 0), n_ref=int(n_ref),
        paso_rampa=float(paso_rampa or 0), espera_I=float(espera_I),
        preparar_desde_uA=float(preparar_desde * 1e6) if preparar_desde is not None else "",
    )
    md_comun.update(extra_metadata)

    v = _Vivo(corrientes, p_drive, f_drive, p_mapa=p_mapa, cada=refresco) if vivo else None
    runs = []
    centro = float(centro_busqueda)
    t0 = time.time()
    hechos, n_tot = 0, len(corrientes) * len(p_drive) * len(f_drive)
    try:
        sg.rf_output(False)
        if preparar_desde is not None:
            print(f"preparando: rampa a {preparar_desde*1e6:.2f} uA")
            fuente.rampa(preparar_desde, paso=paso_rampa, espera=espera_rampa)
            time.sleep(espera_I)

        b = None
        for kI, I in enumerate(corrientes):
            # --- nueva corriente de la bobina ---
            fuente.rampa(I, paso=paso_rampa, espera=espera_rampa)
            time.sleep(espera_I)
            sg.rf_output(False)
            busco = toca_buscar(kI)
            if busco:
                # búsqueda del mínimo (nunca al cambiar el drive); después el VNA queda en CW
                b = buscar_minimo(pna, num, centro, span_busqueda, puntos_busqueda, ifbw_busqueda,
                                  potencia=p_vna, marcador=marcador, desplazamiento=desplazamiento_cw)
                f_cw, est_cw = configurar_cw(pna, num, puntos_cw, ifbw_cw, b["f_cw_objetivo"])
                if seguir_minimo:
                    centro = b["f_min"]
            # referencia sin drive (el VNA ya está en CW en f_cw)
            zs = [leer_cw(pna, est_cw)[0] for _ in range(max(int(n_ref), 1))]
            z_ref = complex(np.mean(zs))
            if busco:
                print(f"I = {I*1e6:8.3f} uA | mínimo {b['f_min']/1e9:.6f} GHz ({b['s21_min_db']:.2f} dB) | "
                      f"f CW {f_cw/1e9:.6f} GHz | ref {20*np.log10(abs(z_ref)):.2f} dB")
            else:
                print(f"I = {I*1e6:8.3f} uA | f CW fija {f_cw/1e9:.6f} GHz | "
                      f"ref {20*np.log10(abs(z_ref)):.2f} dB, {np.degrees(np.angle(z_ref)):+.2f}°")
            if v:
                if busco:
                    v.busqueda(I, b, f_cw)
                v.referencia(I, z_ref)

            usar = np.ones(len(f_drive), bool)
            if excluir_cw:
                usar = np.abs(f_drive - f_cw) >= excluir_cw

            with _medicion(f"{nombre}_I{I*1e6:+.2f}uA").run() as saver:
                ds = saver.dataset
                rid = ds.run_id
                # f_min_marcador / s21_min_db son de la última búsqueda (busqueda_aqui = 0 si no se
                # buscó en esta corriente, y entonces el run no tiene traza de búsqueda)
                md = dict(md_comun, dos_tonos_idx=kI, current_uA=float(I * 1e6),
                          busqueda_aqui=int(busco),
                          f_min_marcador=b["f_min"], s21_min_db=b["s21_min_db"], f_cw=f_cw,
                          ref_re=z_ref.real, ref_im=z_ref.imag,
                          ref_mag_db=float(20 * np.log10(abs(z_ref))),
                          ref_phase_deg=float(np.degrees(np.angle(z_ref))),
                          timestamp_inicio=time.strftime("%Y-%m-%d %H:%M:%S"),
                          vna_cw_state=json.dumps({k: x for k, x in est_cw.items() if k != "trazas"},
                                                  default=str),
                          vna_trazas=json.dumps(est_cw["trazas"], default=str))
                try:
                    md["sg_idn"] = sg.ask("*IDN?").strip()
                except Exception:
                    pass
                for k, x in md.items():
                    ds.add_metadata(k, x)
                if busco:
                    zb = b["z"]
                    saver.add_result(("f_busqueda", b["f"]), ("busqueda_re", zb.real),
                                     ("busqueda_im", zb.imag),
                                     ("busqueda_mag", 20 * np.log10(np.abs(zb))))
                runs.append(dict(run_id=rid, current_uA=float(I * 1e6), f_min=b["f_min"], f_cw=f_cw))

                for kP, P in enumerate(p_drive):
                    sg.power(float(P))
                    fila_m = np.full(len(f_drive), np.nan)
                    fila_p = np.full(len(f_drive), np.nan)
                    sg.frequency(float(f_drive[np.argmax(usar)] if usar.any() else f_drive[0]))
                    sg.rf_output(True)
                    for kf, fd in enumerate(f_drive):
                        hechos += 1
                        if not usar[kf]:
                            continue
                        sg.frequency(float(fd))
                        if espera_drive:
                            time.sleep(espera_drive)
                        z, err = leer_cw(pna, est_cw)
                        saver.add_result(("p_drive", float(P)), ("f_drive", float(fd)),
                                         ("s21_re", z.real), ("s21_im", z.imag),
                                         ("s21_mag", 20 * np.log10(abs(z))),
                                         ("s21_phase", np.degrees(np.angle(z))), ("s21_err", err))
                        fila_m[kf] = 20 * np.log10(abs(z) / abs(z_ref))
                        fila_p[kf] = np.degrees(np.angle(z / z_ref))
                        if v:
                            el = time.time() - t0
                            falta = el / hechos * (n_tot - hechos)
                            tit = (f"run {rid} | I = {I*1e6:.2f} uA ({kI+1}/{len(corrientes)}) | "
                                   f"P drive = {P:g} dBm ({kP+1}/{len(p_drive)}) | "
                                   f"f drive = {fd/1e9:.4f} GHz | {el/60:.1f} min, faltan ~{falta/60:.0f} min")
                            v.punto(kI, kP, fila_m, fila_p, tit)
                    sg.rf_output(False)
                    if v and usar.any():        # al terminar cada potencia, dibuja sí o sí
                        v.punto(kI, kP, fila_m, fila_p, tit, forzar=True)
                ds.add_metadata("timestamp_fin", time.strftime("%Y-%m-%d %H:%M:%S"))
            print(f"  run {rid} listo")
    except KeyboardInterrupt:
        print("interrumpido por el usuario: lo medido hasta acá quedó guardado")
    finally:
        if apagar_drive_al_final:
            try:
                sg.rf_output(False)
            except Exception as e:
                print("no pude apagar el SG:", e)
        try:
            # el VNA vuelve a barrer el resonador en continuo, para verlo en pantalla
            configurar_busqueda(pna, centro, span_busqueda, puntos_busqueda, ifbw_busqueda)
            pna.sweep_mode("CONT")
        except Exception as e:
            print("no pude devolver el VNA a barrido continuo:", e)
        if volver_a is not None:
            fuente.rampa(volver_a, paso=paso_rampa, espera=espera_rampa)
    if runs:
        print(f"runs {runs[0]['run_id']} a {runs[-1]['run_id']} (dos_tonos = '{tag}')")
    return dict(tag=tag, runs=runs)


# =============================================================================
# 3b) Un tono: barrido de corriente de la bobina leyendo en CW, sin drive
# =============================================================================
class _VivoUnTono:
    """Ventana en vivo del barrido a un tono.

    Arriba: última traza de búsqueda con el mínimo / f del mínimo vs corriente.
    Abajo : |S21| y fase en CW vs corriente de la bobina.
    """

    def __init__(self, corrientes, cada=0.3):
        import matplotlib.pyplot as plt
        self.plt = plt
        self.cada = cada
        self.t_ult = 0.0
        self.I, self.mag, self.fase, self.Ib, self.fmin = [], [], [], [], []
        plt.ion()
        fig, ax = plt.subplots(2, 2, figsize=(13, 8))
        self.fig = fig
        self.ax_b, self.ax_c, self.ax_m, self.ax_p = ax.ravel()
        self.ln_b, = self.ax_b.plot([], [], lw=1)
        self.pt_b, = self.ax_b.plot([], [], "rv", ms=7, label="mínimo (marcador)")
        self.vl_b = self.ax_b.axvline(np.nan, color="r", ls=":", lw=1, label="f CW")
        self.ax_b.set_xlabel("f VNA [GHz]"); self.ax_b.set_ylabel("|S21| [dB]")
        self.ax_b.legend(fontsize=8, loc="lower right")
        self.ln_c, = self.ax_c.plot([], [], "o-", ms=3, lw=0.8)
        self.ax_c.set_ylabel("f mínimo − primera [kHz]")
        self.ln_m, = self.ax_m.plot([], [], ".-", ms=3, lw=0.6)
        self.ln_p, = self.ax_p.plot([], [], ".-", ms=3, lw=0.6, color="C3")
        self.ax_m.set_ylabel("|S21| CW [dB]"); self.ax_p.set_ylabel("fase CW [°]")
        I = np.asarray(corrientes, float) * 1e6
        lo, hi = I.min(), I.max()
        if hi == lo:
            lo, hi = lo - 0.5, hi + 0.5
        for a in (self.ax_c, self.ax_m, self.ax_p):
            a.set_xlabel("I bobina [uA]"); a.set_xlim(lo, hi); a.grid(alpha=0.3)
        self.titulo = fig.suptitle("")
        fig.tight_layout()
        plt.show(block=False)
        self._dibujar()

    def _dibujar(self):
        self.fig.canvas.draw_idle()
        self.fig.canvas.flush_events()
        self.t_ult = time.time()

    def busqueda(self, I, b, f_cw):
        self.ln_b.set_data(b["f"] / 1e9, 20 * np.log10(np.abs(b["z"])))
        self.pt_b.set_data([b["f_min"] / 1e9], [b["s21_min_db"]])
        self.vl_b.set_xdata([f_cw / 1e9, f_cw / 1e9])
        self.ax_b.relim(); self.ax_b.autoscale_view()
        self.ax_b.set_title(f"búsqueda en I = {I*1e6:.3f} uA | f CW = {f_cw/1e9:.6f} GHz", fontsize=10)
        self.Ib.append(I * 1e6); self.fmin.append(b["f_min"])
        self.ln_c.set_data(self.Ib, (np.array(self.fmin) - self.fmin[0]) / 1e3)
        self.ax_c.relim(); self.ax_c.autoscale_view(scalex=False)

    def punto(self, I, z, titulo, forzar=False):
        self.I.append(I * 1e6)
        self.mag.append(20 * np.log10(abs(z)))
        self.fase.append(np.angle(z))
        if not forzar and time.time() - self.t_ult < self.cada:
            return
        self.ln_m.set_data(self.I, self.mag)
        self.ln_p.set_data(self.I, np.degrees(np.unwrap(self.fase)))
        for a in (self.ax_m, self.ax_p):
            a.relim(); a.autoscale_view(scalex=False)
        self.titulo.set_text(titulo)
        self._dibujar()


def _medicion_un_tono(nombre, guardar_busqueda):
    from qcodes.dataset import Measurement
    m = Measurement(name=nombre)
    m.write_period = 5
    m.register_custom_parameter("I_bobina", label="corriente de la bobina", unit="A")
    sp = ("I_bobina",)
    for p, lab, u in (("s21_re", "Re S21 (CW)", ""), ("s21_im", "Im S21 (CW)", ""),
                      ("s21_mag", "|S21| (CW)", "dB"), ("s21_phase", "fase S21 (CW)", "deg"),
                      ("s21_err", "error del promedio de S21", ""),
                      ("f_cw", "f CW del VNA", "Hz"),
                      ("f_min", "f del mínimo (marcador; NaN si no se buscó)", "Hz"),
                      ("t", "tiempo desde el inicio", "s")):
        m.register_custom_parameter(p, label=lab, unit=u, setpoints=sp)
    if guardar_busqueda:
        m.register_custom_parameter("I_busqueda", label="corriente de la bobina (búsqueda)", unit="A")
        m.register_custom_parameter("f_busqueda", label="f VNA (búsqueda del mínimo)", unit="Hz")
        sb = ("I_busqueda", "f_busqueda")
        m.register_custom_parameter("busqueda_re", label="Re S21 búsqueda", unit="", setpoints=sb)
        m.register_custom_parameter("busqueda_im", label="Im S21 búsqueda", unit="", setpoints=sb)
        m.register_custom_parameter("busqueda_mag", label="|S21| búsqueda", unit="dB", setpoints=sb)
    return m


def tiempo_estimado_un_tono(n_corrientes, puntos_cw, ifbw_cw, espera_I, n_busquedas,
                            puntos_busqueda, ifbw_busqueda, overhead=0.1):
    """Duración aproximada [s] (sin contar la rampa de la bobina)."""
    return (n_corrientes * (puntos_cw / ifbw_cw + espera_I + overhead)
            + n_busquedas * (puntos_busqueda / ifbw_busqueda + 1.0))


def medir_un_tono(pna, fuente, corrientes, centro_busqueda, span_busqueda,
                  puntos_busqueda=201, ifbw_busqueda=100, p_vna=-10,
                  puntos_cw=10, ifbw_cw=100, f_cw=None, buscar_cada=0,
                  desplazamiento_cw=0.0, seguir_minimo=True, marcador=1, guardar_busqueda=True,
                  paso_rampa=1e-6, espera_rampa=0.05, espera_I=0.5, preparar_desde=None,
                  volver_a=None, sg=None, vivo=True, refresco=0.5,
                  nombre="un_tono", tag=None, **extra_metadata):
    """Un tono: barre la corriente de la bobina y en cada corriente lee S21 en CW (sin drive).

    corrientes [A]  : corrientes de la bobina, en el orden en que se miden.
    f_cw [Hz]       : si se da, la f CW queda fija ahí y no se busca el mínimo nunca.
    buscar_cada     : (solo si f_cw es None) cada cuántas corrientes se busca el mínimo y se mueve
                      la f CW a él. 0 = solo en la primera corriente (f CW fija todo el barrido);
                      1 = en todas (la f CW sigue al resonador); n = cada n corrientes.
    centro_busqueda, span_busqueda, puntos_busqueda, ifbw_busqueda: barrido LIN de la búsqueda.
    p_vna [dBm]     : potencia del VNA (búsqueda y CW).
    puntos_cw, ifbw_cw: puntos CW promediados en cada corriente y su IFBW.
    desplazamiento_cw [Hz]: f CW = f del mínimo + desplazamiento.
    seguir_minimo   : la próxima búsqueda se centra en el último mínimo.
    guardar_busqueda: guarda las trazas de búsqueda en el run (2D: I_busqueda x f_busqueda).
    preparar_desde [A]: antes de la primera corriente rampa hasta acá (histéresis).
    volver_a [A]    : corriente final de la bobina (None = queda en la última).
    sg              : opcional; si se pasa, se apaga su salida al empezar (por las dudas).
    extra_metadata  : se guarda en el run (resonador, att_vna, notas, ...).

    Todo el barrido es UN run. Devuelve dict(tag, run_id).
    """
    corrientes = np.atleast_1d(np.asarray(corrientes, float))
    if np.any(np.abs(corrientes) > fuente.I_max):
        raise ValueError(f"hay corrientes por encima de I_max = {fuente.I_max*1e6:.1f} uA")
    buscar_cada = int(buscar_cada or 0)
    if f_cw is not None and buscar_cada:
        print("aviso: con f_cw dado no se busca el mínimo (buscar_cada se ignora)")
    nI = len(corrientes)

    def toca_buscar(k):
        if f_cw is not None:
            return False
        return k == 0 or (buscar_cada > 0 and k % buscar_cada == 0)

    n_busq = sum(toca_buscar(k) for k in range(nI))
    tag = tag or time.strftime("%Y%m%d_%H%M%S")
    num = traza_s21_mlog(pna)
    t_tot = tiempo_estimado_un_tono(nI, puntos_cw, ifbw_cw, espera_I, n_busq,
                                    puntos_busqueda, ifbw_busqueda)
    modo = (f"f CW fija en {f_cw/1e9:.6f} GHz" if f_cw is not None else
            "f CW fija (mínimo buscado en la primera corriente)" if buscar_cada == 0 else
            f"f CW sigue al mínimo cada {buscar_cada} corriente(s)")
    print(f"un tono {tag}: {nI} corrientes, {modo}, {n_busq} búsqueda(s) -> ~{t_tot/60:.1f} min "
          f"(sin la rampa de la bobina)")

    md = dict(
        un_tono=tag, n_corrientes=nI,
        corrientes_uA=json.dumps([round(c * 1e6, 5) for c in corrientes]),
        centro_busqueda=float(centro_busqueda), span_busqueda=float(span_busqueda),
        puntos_busqueda=int(puntos_busqueda), ifbw_busqueda=float(ifbw_busqueda),
        pow_vna=float(p_vna), N_points_cw=int(puntos_cw), bw_vna=float(ifbw_cw),
        f_cw_fija=float(f_cw) if f_cw is not None else "", buscar_cada=buscar_cada,
        desplazamiento_cw=float(desplazamiento_cw), seguir_minimo=int(bool(seguir_minimo)),
        marcador=int(marcador), traza_marcador=int(num), guardar_busqueda=int(bool(guardar_busqueda)),
        paso_rampa=float(paso_rampa or 0), espera_I=float(espera_I),
        preparar_desde_uA=float(preparar_desde * 1e6) if preparar_desde is not None else "",
    )
    md.update(extra_metadata)

    v = _VivoUnTono(corrientes, cada=refresco) if vivo else None
    centro = float(centro_busqueda)
    rid = None
    t0 = time.time()
    try:
        if sg is not None:
            sg.rf_output(False)
        if preparar_desde is not None:
            print(f"preparando: rampa a {preparar_desde*1e6:.2f} uA")
            fuente.rampa(preparar_desde, paso=paso_rampa, espera=espera_rampa)
            time.sleep(espera_I)
        est_cw, f_act = None, np.nan
        if f_cw is not None:
            pna.power(float(p_vna))
            pna.cw(float(round(f_cw)))
            f_act, est_cw = configurar_cw(pna, num, puntos_cw, ifbw_cw, float(f_cw))

        with _medicion_un_tono(nombre, guardar_busqueda).run() as saver:
            ds = saver.dataset
            rid = ds.run_id
            md["timestamp_inicio"] = time.strftime("%Y-%m-%d %H:%M:%S")
            for k, x in md.items():
                ds.add_metadata(k, x)
            aviso = max(nI // 20, 1)
            for kI, I in enumerate(corrientes):
                fuente.rampa(I, paso=paso_rampa, espera=espera_rampa)
                time.sleep(espera_I)
                f_min = np.nan
                if toca_buscar(kI):
                    b = buscar_minimo(pna, num, centro, span_busqueda, puntos_busqueda, ifbw_busqueda,
                                      potencia=p_vna, marcador=marcador,
                                      desplazamiento=desplazamiento_cw)
                    f_act, est_cw = configurar_cw(pna, num, puntos_cw, ifbw_cw, b["f_cw_objetivo"])
                    f_min = b["f_min"]
                    if seguir_minimo:
                        centro = f_min
                    if guardar_busqueda:
                        zb = b["z"]
                        saver.add_result(("I_busqueda", np.full(len(b["f"]), I)),
                                         ("f_busqueda", b["f"]), ("busqueda_re", zb.real),
                                         ("busqueda_im", zb.imag),
                                         ("busqueda_mag", 20 * np.log10(np.abs(zb))))
                    if kI == 0:
                        ds.add_metadata("f_cw_inicial", float(f_act))
                    print(f"I = {I*1e6:9.4f} uA | mínimo {f_min/1e9:.6f} GHz ({b['s21_min_db']:.2f} dB) "
                          f"| f CW {f_act/1e9:.6f} GHz")
                    if v:
                        v.busqueda(I, b, f_act)
                z, err = leer_cw(pna, est_cw)
                el = time.time() - t0
                saver.add_result(("I_bobina", float(I)), ("s21_re", z.real), ("s21_im", z.imag),
                                 ("s21_mag", 20 * np.log10(abs(z))),
                                 ("s21_phase", np.degrees(np.angle(z))), ("s21_err", err),
                                 ("f_cw", float(f_act)), ("f_min", float(f_min)), ("t", el))
                falta = el / (kI + 1) * (nI - kI - 1)
                if (kI + 1) % aviso == 0 or kI == nI - 1:
                    print(f"  {kI+1}/{nI}  I = {I*1e6:9.4f} uA  |S21| {20*np.log10(abs(z)):7.2f} dB  "
                          f"fase {np.degrees(np.angle(z)):7.2f}°  | {el/60:.1f} min, faltan ~{falta/60:.1f} min")
                if v:
                    tit = (f"run {rid} | I = {I*1e6:.4f} uA ({kI+1}/{nI}) | f CW = {f_act/1e9:.6f} GHz | "
                           f"{el/60:.1f} min, faltan ~{falta/60:.1f} min")
                    v.punto(I, z, tit, forzar=(kI == nI - 1))
            ds.add_metadata("timestamp_fin", time.strftime("%Y-%m-%d %H:%M:%S"))
        print(f"run {rid} listo (un_tono = '{tag}')")
    except KeyboardInterrupt:
        print("interrumpido por el usuario: lo medido hasta acá quedó guardado")
    finally:
        try:
            configurar_busqueda(pna, centro, span_busqueda, puntos_busqueda, ifbw_busqueda)
            pna.sweep_mode("CONT")
        except Exception as e:
            print("no pude devolver el VNA a barrido continuo:", e)
        if volver_a is not None:
            fuente.rampa(volver_a, paso=paso_rampa, espera=espera_rampa)
    return dict(tag=tag, run_id=rid)


# =============================================================================
# 4) Carga
# =============================================================================
def cargar_dos_tonos(run_id):
    """Run de dos tonos -> dict(p_drive, f_drive, z, err, z_ref, f_busqueda, z_busqueda, metadata).

    p_drive y f_drive son las columnas en el orden en que se midieron (una fila por punto).
    """
    from qcodes.dataset import load_by_id
    ds = load_by_id(run_id)
    d = ds.get_parameter_data()
    md = ds.metadata
    r = d["s21_re"]
    out = dict(
        p_drive=np.ravel(r["p_drive"]).astype(float), f_drive=np.ravel(r["f_drive"]).astype(float),
        z=np.ravel(r["s21_re"]) + 1j * np.ravel(d["s21_im"]["s21_im"]),
        err=np.ravel(d["s21_err"]["s21_err"]),
        z_ref=float(md.get("ref_re", np.nan)) + 1j * float(md.get("ref_im", np.nan)),
        metadata=md,
    )
    # con buscar_cada = 0 solo el primer run tiene traza de búsqueda; los demás quedan vacíos
    if "busqueda_re" in d and len(np.ravel(d["busqueda_re"]["busqueda_re"])):
        out["f_busqueda"] = np.ravel(d["busqueda_re"]["f_busqueda"]).astype(float)
        out["z_busqueda"] = (np.ravel(d["busqueda_re"]["busqueda_re"])
                             + 1j * np.ravel(d["busqueda_im"]["busqueda_im"]))
    else:
        out["f_busqueda"] = np.array([])
        out["z_busqueda"] = np.array([], complex)
    return out


def cargar_un_tono(run_id):
    """Run de un tono -> dict(I [A], z, err, f_cw, f_min, t, metadata) y, si se guardaron,
    I_busqueda [n_b], f_busqueda [n_b x puntos], z_busqueda [n_b x puntos] (una fila por búsqueda).
    """
    from qcodes.dataset import load_by_id
    ds = load_by_id(run_id)
    d = ds.get_parameter_data()
    md = ds.metadata
    r = d["s21_re"]
    col = lambda p: np.ravel(d[p][p]).astype(float)
    out = dict(I=np.ravel(r["I_bobina"]).astype(float),
               z=np.ravel(r["s21_re"]) + 1j * np.ravel(d["s21_im"]["s21_im"]),
               err=col("s21_err"), f_cw=col("f_cw"), f_min=col("f_min"), t=col("t"), metadata=md)
    if "busqueda_re" in d:
        Ib = np.ravel(d["busqueda_re"]["I_busqueda"]).astype(float)
        fb = np.ravel(d["busqueda_re"]["f_busqueda"]).astype(float)
        zb = np.ravel(d["busqueda_re"]["busqueda_re"]) + 1j * np.ravel(d["busqueda_im"]["busqueda_im"])
        n = int(md.get("puntos_busqueda", 0) or 0)
        if n and len(fb) % n == 0:
            out.update(I_busqueda=Ib.reshape(-1, n)[:, 0], f_busqueda=fb.reshape(-1, n),
                       z_busqueda=zb.reshape(-1, n))
        else:
            out.update(I_busqueda=Ib, f_busqueda=fb, z_busqueda=zb)
    return out


def buscar_dos_tonos(tag):
    """Runs de una medición de dos tonos (por su tag), ordenados como se midieron."""
    from qcodes.dataset import experiments
    runs = []
    for e in experiments():
        for ds in e.data_sets():
            md = ds.metadata
            if md.get("dos_tonos") == tag:
                runs.append(dict(run_id=ds.run_id, current_uA=md.get("current_uA"),
                                 f_cw=md.get("f_cw"), _i=md.get("dos_tonos_idx", 0)))
    runs.sort(key=lambda x: x.pop("_i"))
    return dict(tag=tag, runs=runs)
