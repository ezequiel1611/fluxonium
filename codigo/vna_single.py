r"""
Lectura "single" del PNA (Keysight N5245A, driver QCoDeS) tal cual está configurado en pantalla.

NO toca la configuración del VNA (start/stop/puntos/potencia/IFBW/promedios/formato de trazas):
lee el estado, dispara UN barrido (o lee lo que ya está en pantalla) y guarda.

Qué se guarda por cada medición:
    - f_vna [Hz]                         eje de frecuencia del barrido
    - s21_mag [dB], s21_phase [deg]      de los datos COMPLEJOS (SDATA) de la traza S21
      s21_re, s21_im                     -> formato que ya reconoce fit_resonancia.cargar_traza
                                         (SIN electrical delay: la fase incluye el cable)
    - trN_<FORMATO>                      cada traza en el formato que tiene en pantalla (FDATA),
                                         p.ej. tr1_MLOG, tr2_UPH  -> "exactamente lo que ves"
    - metadata: todo el estado del VNA (pow_vna, bw_vna, averages, electrical delay, ...)

Uso típico (en el notebook, con `pna` ya creado):

    from vna_single import medir_single, guardar_npz

    r = medir_single(pna, att_vna=-62)            # barre una vez, guarda en la base activa
    print(r['run_id'])
    guardar_npz(r, 'data')                        # opcional: copia en .npz

    # si el VNA ya está en HOLD y querés lo que está en pantalla sin re-barrer:
    r = medir_single(pna, disparar=False)

Postprocesado:
    from fit_resonancia import ajustar_resonancia
    res = ajustar_resonancia(db, r['run_id'], att_total=-112)
"""

import json
import os
import time

import numpy as np

# Formatos de CALC:FORM que devuelven 2 valores por punto (re, im)
_FORMATOS_COMPLEJOS = {"SMIT", "SADM", "POL", "MLIN_PHAS"}


# =============================================================================
# 1) Estado del VNA
# =============================================================================
def _catalogo(pna):
    """[(nombre_medicion, 'S21'), ...] de las trazas del canal 1."""
    cat = pna.ask("CALC:PAR:CAT:EXT?").strip().strip('"').split(",")
    return list(zip(cat[::2], cat[1::2]))


def _info_trazas(pna):
    """Lista de dicts {nombre, sparam, num, formato, edelay} para cada traza del canal.

    El electrical delay es POR TRAZA en el PNA (p.ej. 0 en la de magnitud y 62 ns en la
    de fase), por eso se lee para cada una.
    """
    activa = pna.active_trace()
    trazas = []
    for nombre, sparam in _catalogo(pna):
        pna.write(f"CALC:PAR:SEL '{nombre}'")
        num = int(pna.ask("CALC:PAR:MNUM?"))
        formato = pna.ask("CALC:FORM?").strip()
        edelay = float(pna.ask("CALC:CORR:EDEL:TIME?"))
        trazas.append(dict(nombre=nombre, sparam=sparam.upper(), num=num, formato=formato,
                           edelay=edelay))
    pna.active_trace(activa)
    return trazas


def leer_estado(pna):
    """Lee (sin modificar) la configuración actual del VNA."""
    est = dict(
        sweep_type=pna.sweep_type(),
        sweep_mode=pna.sweep_mode(),
        start=pna.start(),
        stop=pna.stop(),
        center=pna.center(),
        span=pna.span(),
        points=int(pna.points()),
        if_bandwidth=pna.if_bandwidth(),
        power=pna.power(),
        averages_enabled=bool(pna.averages_enabled()),
        averages=int(pna.averages()),
        sweep_time=pna.sweep_time(),
        trigger_source=pna.trigger_source(),
    )
    if est["sweep_type"] == "CW":
        est["cw"] = pna.cw()
    est["trazas"] = _info_trazas(pna)
    return est


def resumen(est):
    """Una línea con la configuración leída del VNA."""
    prom = est["averages"] if est["averages_enabled"] else 1
    return (f"{est['sweep_type']} {est['start']/1e9:.6f}-{est['stop']/1e9:.6f} GHz "
            f"(centro {est['center']/1e9:.6f} GHz, span {est['span']/1e6:g} MHz) | "
            f"{est['points']} pts | IFBW {est['if_bandwidth']:g} Hz | P {est['power']:g} dBm | "
            f"prom {prom} | trazas: "
            + ", ".join(f"{t['sparam']}/{t['formato']}"
                        + (f" (delay {t['edelay']*1e9:.3f} ns)" if t["edelay"] else "")
                        for t in est["trazas"]))


# =============================================================================
# 1b) Configuración del VNA desde el notebook
# =============================================================================
def configurar_vna(pna, f_center=None, span=None, power=None, if_bandwidth=None, points=None,
                   averages=None, edelay_fase=None, sweep_type="LIN", continuo=True,
                   autoescala=True, verbose=True):
    """Configura el VNA. Cualquier argumento en None NO se toca (queda como en el panel).

    f_center, span [Hz] | power [dBm] | if_bandwidth [Hz] | points | averages (1 = apagados)
    edelay_fase [s]: electrical delay de las trazas S21 en formato de fase (PHAS/UPH/PPH).
    continuo=True  : deja el VNA barriendo en continuo para que veas el resultado en pantalla.
    autoescala=True: autoescala el eje y de todas las trazas.

    Devuelve el estado leído del VNA después de configurar (el VNA puede redondear valores,
    p.ej. el IFBW a los valores permitidos).
    """
    if sweep_type is not None:
        pna.sweep_type(sweep_type)

    if f_center is not None or span is not None:
        c = pna.center() if f_center is None else float(f_center)
        s = pna.span() if span is None else float(span)
        if s <= 0:
            raise ValueError("span tiene que ser > 0")
        if c - s / 2 < pna.min_freq or c + s / 2 > pna.max_freq:
            raise ValueError(f"{(c - s/2)/1e9:.6f}-{(c + s/2)/1e9:.6f} GHz fuera del rango del "
                             f"VNA ({pna.min_freq/1e9:g}-{pna.max_freq/1e9:g} GHz)")
        # directo por SCPI: el validador del driver no deja spans < 10 MHz
        pna.write(f"SENS:FREQ:CENT {c:.6f}")
        pna.write(f"SENS:FREQ:SPAN {s:.6f}")

    if points is not None:
        pna.points(int(points))
    if if_bandwidth is not None:
        pna.if_bandwidth(float(if_bandwidth))
    if power is not None:
        pna.power(float(power))
    if averages is not None:
        n = int(averages)
        if n > 1:
            pna.averages(n)
            pna.averages_enabled(True)
        else:
            pna.averages_enabled(False)

    trazas = _info_trazas(pna)
    if edelay_fase is not None:
        activa = pna.active_trace()
        for t in trazas:
            if t["sparam"] == "S21" and t["formato"] in _FORMATOS_FASE:
                pna.active_trace(t["num"])
                pna.write(f"CALC:CORR:EDEL:TIME {float(edelay_fase):.6e}")
        pna.active_trace(activa)

    if continuo:
        pna.sweep_mode("CONT")
    if autoescala:
        for t in trazas:
            pna.write(f"DISP:MEAS{t['num']}:Y:AUTO")

    pna.ask("*OPC?")                       # espera a que el VNA aplique todo
    est = leer_estado(pna)
    if verbose:
        print(resumen(est))
        n_barr = est["averages"] if est["averages_enabled"] else 1
        print(f"tiempo de barrido ~{est['sweep_time']:.2f} s"
              + (f" x {n_barr} promedios = {est['sweep_time'] * n_barr:.1f} s" if n_barr > 1 else ""))
    return est


def eje_x(pna, est):
    """Eje x del barrido. Frecuencia en Hz (LIN/LOG/SEGM) o tiempo en s (CW)."""
    n, tipo = est["points"], est["sweep_type"]
    if tipo == "LIN":
        return np.linspace(est["start"], est["stop"], n)
    if tipo == "LOG":
        return np.geomspace(est["start"], est["stop"], n)
    if tipo == "CW":
        return np.linspace(0, est["sweep_time"], n)
    # SEGM u otros: se lo pedimos al VNA en ASCII (REAL,32 pierde resolución en Hz)
    pna.write("FORM ASC,0")
    try:
        return np.array(pna.ask("SENS:X?").split(","), dtype=float)
    finally:
        pna.write("FORM REAL,32")
        pna.write("FORM:BORD NORM")


# =============================================================================
# 2) Disparo de un barrido y lectura de trazas
# =============================================================================
def disparar_single(pna, est, timeout=None):
    """Un barrido completo (o un grupo = n° de promedios si están activados).

    Igual que traces.tr1.run_sweep() del driver: deja el VNA en HOLD.
    Devuelve el sweep_mode previo para poder restaurarlo.
    """
    previo = pna.sweep_mode()
    if est["averages_enabled"]:
        pna.reset_averages()
        pna.group_trigger_count(est["averages"])
        pna.sweep_mode("GRO")
    else:
        pna.sweep_mode("SING")

    n_barridos = est["averages"] if est["averages_enabled"] else 1
    t_est = est["sweep_time"] * n_barridos
    if timeout is None:
        timeout = 3 * t_est + 30
    t0 = time.time()
    while pna.sweep_mode() != "HOLD":
        if time.time() - t0 > timeout:
            raise TimeoutError(f"el VNA no terminó el barrido en {timeout:.0f} s "
                               f"(trigger_source={est['trigger_source']})")
        time.sleep(min(0.1, max(t_est / 50, 0.02)))
    return previo


def _leer_binario(pna, cmd):
    return np.array(pna.visa_handle.query_binary_values(cmd, datatype="f", is_big_endian=True))


def leer_trazas(pna, est):
    """Lee cada traza: formateada (FDATA, lo de pantalla) y compleja (SDATA). No barre."""
    activa = pna.active_trace()
    out = []
    try:
        for t in est["trazas"]:
            pna.active_trace(t["num"])
            fdata = _leer_binario(pna, "CALC:DATA? FDATA")
            sdata = _leer_binario(pna, "CALC:DATA? SDATA")
            if t["formato"] in _FORMATOS_COMPLEJOS or fdata.size == 2 * est["points"]:
                fdata = fdata[0::2] + 1j * fdata[1::2]
            out.append(dict(t, fdata=fdata, sdata=sdata[0::2] + 1j * sdata[1::2]))
    finally:
        pna.active_trace(activa)
    return out


def adquirir(pna, disparar=True, restaurar=True, timeout=None):
    """Lee estado + (opcional) dispara un barrido + lee todas las trazas. No guarda nada.

    disparar=True : hace UN barrido nuevo (SING o GRO con los promedios configurados).
    disparar=False: lee lo que hay en pantalla (útil si el VNA ya está en HOLD).
    restaurar=True: al final devuelve el VNA al sweep_mode que tenía (p.ej. CONT).
    """
    est = leer_estado(pna)
    x = eje_x(pna, est)
    previo = None
    t_ini = time.strftime("%Y-%m-%d %H:%M:%S")
    try:
        if disparar:
            previo = disparar_single(pna, est, timeout=timeout)
        trazas = leer_trazas(pna, est)
    finally:
        if restaurar and previo is not None:
            pna.sweep_mode(previo)
    return dict(estado=est, x=x, trazas=trazas, timestamp=t_ini, disparado=disparar)


def traza_s21(r, sparam="S21"):
    """Primera traza que mide `sparam` (por defecto S21)."""
    for t in r["trazas"]:
        if t["sparam"] == sparam.upper():
            return t
    raise ValueError(f"no hay ninguna traza {sparam} en el VNA: "
                     f"{[(t['nombre'], t['sparam']) for t in r['trazas']]}")


_FORMATOS_FASE = ("PHAS", "UPH", "PPH")


def fase_pantalla(r, sparam="S21"):
    """Fase [deg] tal como la muestra el VNA (con el electrical delay de esa traza).

    Si hay una traza `sparam` en formato de fase (PHAS/UPH/PPH) usa sus datos de pantalla.
    Si no, la calcula de los datos complejos aplicando el electrical delay de la traza
    `sparam` y la deja entre ±180° (como PHAS).
    Devuelve (fase_deg, edelay_s, origen).
    """
    for t in r["trazas"]:
        if t["sparam"] == sparam.upper() and t["formato"] in _FORMATOS_FASE:
            return np.asarray(t["fdata"], float), t["edelay"], f"tr{t['num']} {t['formato']}"
    t = traza_s21(r, sparam)
    z = t["sdata"] * np.exp(2j * np.pi * r["x"] * t["edelay"])
    return np.rad2deg(np.angle(z)), t["edelay"], f"SDATA + delay de tr{t['num']}"


# =============================================================================
# 3) Guardado
# =============================================================================
def _metadata(r, extra):
    est = r["estado"]
    md = dict(
        pow_vna=est["power"], bw_vna=est["if_bandwidth"], N_points=est["points"],
        f_start=est["start"], f_stop=est["stop"], f0_vna=est["center"], span_vna=est["span"],
        averages=est["averages"] if est["averages_enabled"] else 1,
        sweep_type=est["sweep_type"],
        sweep_time=est["sweep_time"], disparado=int(r["disparado"]), timestamp=r["timestamp"],
        trazas=json.dumps([{k: t[k] for k in ("nombre", "sparam", "num", "formato", "edelay")}
                           for t in r["trazas"]]),
        vna_state=json.dumps({k: v for k, v in est.items() if k != "trazas"}),
    )
    try:                                   # delay de la traza de fase S21 (lo que ves en pantalla)
        md["electrical_delay"] = fase_pantalla(r)[1]
    except ValueError:
        pass
    md.update(extra or {})
    return md


def guardar_db(r, sparam="S21", nombre="vna_single", **extra_metadata):
    """Guarda en la base QCoDeS activa (initialise_or_create_database_at + experimento).

    extra_metadata: lo que quieras agregar, p.ej. att_vna=-62, current_uA=18.3, muestra='r3'.
    Devuelve el run_id.
    """
    from qcodes.dataset import Measurement

    est = r["estado"]
    if est["sweep_type"] not in ("LIN", "LOG", "SEGM"):
        raise ValueError(f"sweep_type={est['sweep_type']}: el eje x no es frecuencia, "
                         "usá adquirir() y guardá a mano (o guardar_npz)")
    s21 = traza_s21(r, sparam)
    z = s21["sdata"]

    meas = Measurement(name=nombre)
    meas.register_custom_parameter("f_vna", label="VNA frequency", unit="Hz")
    sp = ("f_vna",)
    meas.register_custom_parameter("s21_mag", label=f"|{sparam}|", unit="dB", setpoints=sp)
    meas.register_custom_parameter("s21_phase", label=f"arg {sparam}", unit="deg", setpoints=sp)
    meas.register_custom_parameter("s21_re", label=f"Re {sparam}", unit="", setpoints=sp)
    meas.register_custom_parameter("s21_im", label=f"Im {sparam}", unit="", setpoints=sp)

    resultados = [
        ("f_vna", r["x"]),
        ("s21_mag", 20 * np.log10(np.abs(z))),
        ("s21_phase", np.rad2deg(np.angle(z))),
        ("s21_re", z.real),
        ("s21_im", z.imag),
    ]
    # cada traza como se ve en pantalla (formatos reales; los complejos -> _re/_im)
    for t in r["trazas"]:
        base = f"tr{t['num']}_{t['sparam']}_{t['formato']}"
        lab = f"{t['nombre']} ({t['sparam']}, {t['formato']})"
        if np.iscomplexobj(t["fdata"]):
            for suf, val in (("re", t["fdata"].real), ("im", t["fdata"].imag)):
                meas.register_custom_parameter(f"{base}_{suf}", label=f"{lab} {suf}", setpoints=sp)
                resultados.append((f"{base}_{suf}", val))
        else:
            meas.register_custom_parameter(base, label=lab, setpoints=sp)
            resultados.append((base, t["fdata"]))

    with meas.run() as datasaver:
        for k, v in _metadata(r, extra_metadata).items():
            datasaver.dataset.add_metadata(k, v)
        datasaver.add_result(*resultados)
        run_id = datasaver.run_id
    return run_id


def guardar_npz(r, carpeta="data", nombre=None, sparam="S21", **extra_metadata):
    """Copia en .npz: f, z (complejo S21), cada traza formateada y el estado (json)."""
    os.makedirs(carpeta, exist_ok=True)
    if nombre is None:
        nombre = time.strftime("%Y%m%d_%H%M%S") + "_vna_single"
        if r.get("run_id") is not None:
            nombre += f"_run{r['run_id']}"
    path = os.path.join(carpeta, nombre + ".npz")
    arrays = {f"tr{t['num']}_{t['sparam']}_{t['formato']}": t["fdata"] for t in r["trazas"]}
    try:
        arrays["z"] = traza_s21(r, sparam)["sdata"]
    except ValueError:
        pass
    np.savez(path, f=r["x"], **arrays,
             metadata=json.dumps(_metadata(r, extra_metadata), default=str))
    return path


# =============================================================================
# 4) Todo junto
# =============================================================================
def medir_single(pna, disparar=True, restaurar=True, guardar=True, sparam="S21",
                 nombre="vna_single", graficar=True, **extra_metadata):
    """Adquiere lo que está configurado en el VNA y lo guarda en la base activa.

    Devuelve el dict de adquisición con 'run_id' agregado.
    """
    r = adquirir(pna, disparar=disparar, restaurar=restaurar)
    est = r["estado"]
    print(resumen(est))
    r["run_id"] = guardar_db(r, sparam=sparam, nombre=nombre, **extra_metadata) if guardar else None
    if r["run_id"] is not None:
        print("run_id:", r["run_id"])
    if graficar:
        graficar_adquisicion(r)
    return r


def graficar_adquisicion(r):
    """Un subplot por traza, en el formato de pantalla."""
    import matplotlib.pyplot as plt

    trazas = r["trazas"]
    x = r["x"] / (1e9 if r["estado"]["sweep_type"] != "CW" else 1)
    fig, axs = plt.subplots(len(trazas), 1, sharex=True, figsize=(8, 2.6 * len(trazas)), squeeze=False)
    for ax, t in zip(axs[:, 0], trazas):
        y = t["fdata"]
        if np.iscomplexobj(y):
            ax.plot(x, y.real, label="re"); ax.plot(x, y.imag, label="im"); ax.legend()
        else:
            ax.plot(x, y)
        ax.set_ylabel(f"{t['sparam']} {t['formato']}")
    axs[-1, 0].set_xlabel("f [GHz]" if r["estado"]["sweep_type"] != "CW" else "t [s]")
    titulo = f"run {r['run_id']}" if r.get("run_id") is not None else r["timestamp"]
    axs[0, 0].set_title(titulo)
    fig.tight_layout()
    plt.show()
    return fig


# =============================================================================
# 5) Barrido en potencia (un run por potencia, todos con el mismo tag 'barrido')
# =============================================================================
def _valor_para(x, p):
    """None -> None | función -> x(p) | número -> x."""
    if x is None:
        return None
    return x(p) if callable(x) else x


def barrido_potencia(pna, potencias, ifbw=None, averages=None, resonador="", nombre=None,
                     sparam="S21", graficar=True, tag=None, verbose=True, **extra_metadata):
    """Mide un barrido de frecuencia (tal como está configurado) para cada potencia del VNA.

    potencias: lista de potencias del VNA [dBm]. Conviene ir de menor a mayor
               (a potencia alta el resonador puede ser no lineal / con histéresis).
    ifbw, averages: None (no tocar), un número fijo, o una función de la potencia,
               p.ej. ifbw=lambda p: 10 if p < -30 else 100   (menos IFBW a baja potencia).
    Cada potencia se guarda como un run propio (compatible con fit_resonancia), con
    metadata barrido=<tag>, barrido_idx, resonador y lo que pases en extra_metadata.
    Al final devuelve el VNA a la potencia/IFBW/promedios/modo que tenía.

    tag: etiqueta del barrido (None = fecha y hora). verbose=False: no imprime cada potencia.

    Devuelve dict(barrido=tag, runs=[{run_id, power, ifbw, averages}, ...]).
    """
    potencias = [float(p) for p in np.atleast_1d(potencias)]
    previo = dict(power=pna.power(), ifbw=pna.if_bandwidth(),
                  avg_on=bool(pna.averages_enabled()), avg=int(pna.averages()),
                  modo=pna.sweep_mode())
    tag = tag or time.strftime("%Y%m%d_%H%M%S")
    nombre = nombre or (f"barrido_pot_{resonador}" if resonador else "barrido_pot")
    runs, curvas = [], []
    if verbose:
        print(f"barrido {tag}: {len(potencias)} potencias de {potencias[0]:g} a {potencias[-1]:g} dBm")
    t0 = time.time()
    try:
        for i, p in enumerate(potencias):
            est = configurar_vna(pna, power=p, if_bandwidth=_valor_para(ifbw, p),
                                 averages=_valor_para(averages, p), sweep_type=None,
                                 continuo=False, autoescala=False, verbose=False)
            r = adquirir(pna, disparar=True, restaurar=False)
            rid = guardar_db(r, sparam=sparam, nombre=nombre, barrido=tag, barrido_idx=i,
                             resonador=resonador, **extra_metadata)
            n_prom = est["averages"] if est["averages_enabled"] else 1
            runs.append(dict(run_id=rid, power=est["power"], ifbw=est["if_bandwidth"],
                             averages=n_prom))
            z = traza_s21(r, sparam)["sdata"]
            curvas.append((est["power"], r["x"], 20 * np.log10(np.abs(z))))

            el = time.time() - t0
            falta = el / (i + 1) * (len(potencias) - i - 1)
            if verbose:
                print(f"  {i + 1}/{len(potencias)} | P = {est['power']:6.1f} dBm | IFBW {est['if_bandwidth']:g} Hz"
                      f" | prom {n_prom} | run {rid} | {el / 60:.1f} min, faltan ~{falta / 60:.1f} min")
    finally:
        pna.power(previo["power"])
        pna.if_bandwidth(previo["ifbw"])
        pna.averages(previo["avg"])
        pna.averages_enabled(previo["avg_on"])
        pna.sweep_mode(previo["modo"])

    if runs and verbose:
        print(f"listo: runs {runs[0]['run_id']} a {runs[-1]['run_id']} (barrido = '{tag}')")
    if graficar and curvas:
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(figsize=(8, 4.5))
        cmap = plt.get_cmap("viridis")
        pmin, pmax = min(c[0] for c in curvas), max(c[0] for c in curvas)
        for p, x, mag in curvas:
            ax.plot(x / 1e9, mag, lw=1, color=cmap((p - pmin) / ((pmax - pmin) or 1)))
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(pmin, pmax))
        fig.colorbar(sm, ax=ax, label="P VNA [dBm]")
        ax.set_xlabel("f [GHz]"); ax.set_ylabel(f"|{sparam}| [dB]")
        ax.set_title(f"barrido en potencia {tag} {resonador}")
        fig.tight_layout()
        plt.show()
    return dict(barrido=tag, runs=runs)


def buscar_barrido(barrido):
    """Recupera los runs de un barrido (por su tag) desde la base activa, p.ej. tras reiniciar."""
    from qcodes.dataset import experiments

    runs = []
    for e in experiments():
        for ds in e.data_sets():
            md = ds.metadata
            if md.get("barrido") == barrido:
                runs.append(dict(run_id=ds.run_id, power=md.get("pow_vna"), ifbw=md.get("bw_vna"),
                                 averages=md.get("averages"), _i=md.get("barrido_idx", 0)))
    runs.sort(key=lambda x: x.pop("_i"))
    return dict(barrido=barrido, runs=runs)


def buscar_barrido_corriente(barrido_corriente):
    """Recupera un barrido en potencia x corriente (por su tag) desde la base activa.

    Devuelve dict(barrido_corriente=tag, barridos=[{current_uA, barrido, runs}, ...]) ordenado
    como se midió; cada 'barridos[k]' tiene el mismo formato que devuelve barrido_potencia.
    """
    from qcodes.dataset import experiments

    grupos = {}
    for e in experiments():
        for ds in e.data_sets():
            md = ds.metadata
            if md.get("barrido_corriente") != barrido_corriente:
                continue
            g = grupos.setdefault(md["barrido"], dict(current_uA=md.get("current_uA"),
                                                      barrido=md["barrido"], runs=[],
                                                      _k=md.get("barrido_corriente_idx", 0)))
            g["runs"].append(dict(run_id=ds.run_id, power=md.get("pow_vna"), ifbw=md.get("bw_vna"),
                                  averages=md.get("averages"), _i=md.get("barrido_idx", 0)))
    barridos = sorted(grupos.values(), key=lambda g: g.pop("_k"))
    for g in barridos:
        g["runs"].sort(key=lambda x: x.pop("_i"))
    return dict(barrido_corriente=barrido_corriente, barridos=barridos)
