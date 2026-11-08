r"""
Barrido en corriente de la bobina (flujo) con el PNA y búsqueda de punch-out.

Lo usa punch_out.ipynb:
    1) barrido_corriente(...)  -> UN run 2D por sentido: |S21|(f, I) con la config. actual del VNA
    2) ajustar_barrido_corriente(...) -> fr, Qi, Qc, κ vs corriente (ajuste notch por traza)
    3) en la corriente elegida: vna_single.barrido_potencia(...) y ajustar_runs(...) -> fr, Qi vs P

Formato del run 2D (compatible con fit_resonancia.ajustar_resonancia(db, run, corte={'current': I_uA})):
    current [uA], f_vna [Hz] | s21_mag [dB], s21_phase [deg] (SDATA, sin electrical delay), s21_re, s21_im
"""

import time

import numpy as np

from vna_single import (_metadata, adquirir, eje_x, leer_estado, traza_s21)


# =============================================================================
# 1) Fuente de corriente de la bobina (Keithley, mismos comandos que readout.ipynb)
# =============================================================================
class FuenteCorriente:
    """Keithley de la bobina. I_max [A]: límite de seguridad, no deja setear más que eso."""

    def __init__(self, inst, I_max=500e-6):
        self.inst = inst
        self.I_max = float(I_max)
        self.ultima = None                  # última corriente seteada [A]

    def set_current(self, I):
        I = float(I)
        if abs(I) > self.I_max:
            raise ValueError(f"|I| = {abs(I)*1e6:.1f} uA supera I_max = {self.I_max*1e6:.1f} uA")
        self.inst.write(f"I{I:.4e}X")
        self.ultima = I

    def get_current(self):
        return float(self.inst.query("CURR?").split(",")[0][4:])

    def actual(self):
        """Última corriente seteada (o la leída de la fuente si todavía no se seteó ninguna)."""
        if self.ultima is None:
            try:
                self.ultima = self.get_current()
            except Exception:
                self.ultima = 0.0
        return self.ultima

    def rampa(self, objetivo, paso=1e-6, espera=0.05):
        """Lleva la corriente al objetivo en pasos de `paso` [A] (None = directo)."""
        objetivo = float(objetivo)
        if abs(objetivo) > self.I_max:
            raise ValueError(f"|I| = {abs(objetivo)*1e6:.1f} uA supera I_max = {self.I_max*1e6:.1f} uA")
        if paso is None:
            self.set_current(objetivo)
            return
        paso = abs(float(paso))                 # el sentido lo da objetivo - actual
        if paso == 0:
            raise ValueError("paso = 0: usá un paso > 0 o paso=None para ir directo")
        actual = self.actual()
        n = max(int(np.ceil(abs(objetivo - actual) / paso)), 1)
        for i_ in np.linspace(actual, objetivo, n + 1)[1:]:
            self.set_current(i_)
            time.sleep(espera)

    def ir_a(self, objetivo, desde_abajo=True, margen=10e-6, paso=1e-6, espera=0.05):
        """Llega a `objetivo` siempre desde el mismo lado (histéresis por vórtices / flujo atrapado).

        desde_abajo=True: si la corriente actual es mayor que el objetivo, primero baja a
        objetivo - margen y después sube. desde_abajo=False: lo mismo desde arriba. None: directo.
        """
        actual = self.actual()
        if desde_abajo is True and actual > objetivo:
            self.rampa(objetivo - margen, paso, espera)
        elif desde_abajo is False and actual < objetivo:
            self.rampa(objetivo + margen, paso, espera)
        self.rampa(objetivo, paso, espera)


def conectar_fuente(direccion="GPIB0::15::INSTR", I_max=500e-6, rm=None):
    import pyvisa
    rm = rm or pyvisa.ResourceManager()
    return FuenteCorriente(rm.open_resource(direccion), I_max=I_max)


# =============================================================================
# 2) Barrido en corriente: un run 2D |S21|(f, I)
# =============================================================================
def _medicion_2d(nombre):
    from qcodes.dataset import Measurement
    m = Measurement(name=nombre)
    m.write_period = 2
    m.register_custom_parameter("current", label="corriente bobina", unit="uA")
    m.register_custom_parameter("f_vna", label="frecuencia VNA", unit="Hz")
    sp = ("current", "f_vna")
    m.register_custom_parameter("s21_mag", label="|S21|", unit="dB", setpoints=sp)
    m.register_custom_parameter("s21_phase", label="fase S21 (sin delay)", unit="deg", setpoints=sp)
    m.register_custom_parameter("s21_re", label="Re S21", unit="", setpoints=sp)
    m.register_custom_parameter("s21_im", label="Im S21", unit="", setpoints=sp)
    return m


def _mapa_vivo(corrientes, x):
    """Figura que se actualiza en el notebook mientras se mide."""
    import matplotlib.pyplot as plt
    from IPython.display import display
    plt.ioff()
    fig, (ax_t, ax_m) = plt.subplots(1, 2, figsize=(13, 4.5), gridspec_kw=dict(width_ratios=[1, 1.4]))
    ln, = ax_t.plot([], [], lw=1)
    ax_t.set_xlabel("f [GHz]"); ax_t.set_ylabel("|S21| [dB]")
    orden = np.argsort(corrientes)
    M = np.full((len(corrientes), len(x)), np.nan)
    Is = np.sort(corrientes) * 1e6
    im = ax_m.imshow(M, aspect="auto", origin="lower", cmap="viridis", interpolation="nearest",
                     extent=[x[0] / 1e9, x[-1] / 1e9, Is[0], Is[-1] if len(Is) > 1 else Is[0] + 1])
    fig.colorbar(im, ax=ax_m, label="|S21| − mediana [dB]")
    ax_m.set_xlabel("f [GHz]"); ax_m.set_ylabel("I bobina [uA]")
    fig.tight_layout()
    h = display(fig, display_id=True)
    plt.ion()

    def actualizar(k, I, mag, titulo):
        ln.set_data(x / 1e9, mag)
        ax_t.relim(); ax_t.autoscale_view()
        ax_t.set_title(f"I = {I*1e6:.2f} uA", fontsize=10)
        M[k] = mag - np.median(mag)
        Z = M[orden]
        if np.isfinite(Z).any():
            im.set_data(Z)
            im.set_clim(np.nanpercentile(Z, 1), np.nanpercentile(Z, 99))
        ax_m.set_title(titulo, fontsize=10)
        if h is not None:
            h.update(fig)

    def cerrar():
        plt.close(fig)

    return actualizar, cerrar


def barrido_corriente(pna, fuente, corrientes, espera=0.5, paso_rampa=1e-6, espera_rampa=0.05,
                      resonador="", nombre=None, sparam="S21", ida_y_vuelta=False, volver_a=None,
                      vivo=True, actualizar_cada=1, tag=None, **extra_metadata):
    """Para cada corriente: rampa, espera y un barrido del VNA (tal como está configurado).

    corrientes [A] en el orden en que se miden. ida_y_vuelta=True: después mide las mismas
    corrientes en orden inverso en OTRO run (para ver histéresis).
    volver_a: corriente final [A] (None = queda en la última).
    Devuelve dict(tag, runs=[{run_id, sentido}], corrientes) — runs[0] ida, runs[1] vuelta.
    """
    corrientes = np.asarray(corrientes, float)
    if np.any(np.abs(corrientes) > fuente.I_max):
        raise ValueError(f"hay corrientes por encima de I_max = {fuente.I_max*1e6:.1f} uA")
    est = leer_estado(pna)
    if est["sweep_type"] not in ("LIN", "LOG", "SEGM"):
        raise ValueError(f"sweep_type = {est['sweep_type']}: poné el VNA en barrido de frecuencia")
    x = eje_x(pna, est)
    tag = tag or time.strftime("%Y%m%d_%H%M%S")
    nombre = nombre or (f"barrido_I_{resonador}" if resonador else "barrido_I")
    sentidos = [("ida", corrientes)] + ([("vuelta", corrientes[::-1])] if ida_y_vuelta else [])

    n_prom = est["averages"] if est["averages_enabled"] else 1
    t_traza = est["sweep_time"] * n_prom + espera
    n_tot = sum(len(c) for _, c in sentidos)
    print(f"barrido en corriente {tag}: {n_tot} trazas ({' + '.join(s for s, _ in sentidos)}), "
          f"~{t_traza:.1f} s c/u -> ~{n_tot * t_traza / 60:.0f} min (sin contar la rampa)")

    previo = pna.sweep_mode()
    runs = []
    t0, hechas = time.time(), 0
    try:
        for sentido, cs in sentidos:
            actualizar, cerrar = _mapa_vivo(cs, x) if vivo else (None, None)
            with _medicion_2d(f"{nombre}_{sentido}").run() as saver:
                ds = saver.dataset
                rid = ds.run_id
                runs.append(dict(run_id=rid, sentido=sentido))
                for k, I in enumerate(cs):
                    fuente.rampa(I, paso=paso_rampa, espera=espera_rampa)
                    time.sleep(espera)
                    r = adquirir(pna, disparar=True, restaurar=False)
                    z = traza_s21(r, sparam)["sdata"]
                    if k == 0:
                        md = _metadata(r, dict(extra_metadata, barrido_I=tag, sentido=sentido,
                                               resonador=resonador, espera_I=espera,
                                               I_start_uA=cs[0] * 1e6, I_stop_uA=cs[-1] * 1e6,
                                               n_corrientes=len(cs)))
                        md.pop("timestamp", None)
                        for key, v in md.items():
                            ds.add_metadata(key, v)
                        ds.add_metadata("timestamp_inicio", r["timestamp"])
                    mag = 20 * np.log10(np.abs(z))
                    saver.add_result(("current", np.full(len(x), I * 1e6)), ("f_vna", x),
                                     ("s21_mag", mag), ("s21_phase", np.rad2deg(np.angle(z))),
                                     ("s21_re", z.real), ("s21_im", z.imag))
                    hechas += 1
                    el = time.time() - t0
                    falta = el / hechas * (n_tot - hechas)
                    titulo = (f"run {rid} ({sentido}) | {k + 1}/{len(cs)} | {el/60:.1f} min, "
                              f"faltan ~{falta/60:.0f} min")
                    if vivo and ((k + 1) % actualizar_cada == 0 or k == len(cs) - 1):
                        actualizar(k, I, mag, titulo)
                    elif not vivo and ((k + 1) % 10 == 0 or k == len(cs) - 1):
                        print("  " + titulo)
            if vivo:
                cerrar()
            print(f"run {rid} ({sentido}) listo")
    finally:
        pna.sweep_mode(previo)
        if volver_a is not None:
            fuente.rampa(volver_a, paso=paso_rampa, espera=espera_rampa)
    return dict(tag=tag, runs=runs, corrientes=corrientes)


# =============================================================================
# 3) Carga y ajuste
# =============================================================================
def cargar_barrido_corriente(run_id):
    """Run 2D -> dict(I_uA (orden de medición), f [Hz], Z (n_I x n_f, complejo), metadata)."""
    from qcodes.dataset import load_by_id
    ds = load_by_id(run_id)
    d = ds.get_parameter_data()
    I = np.ravel(d["s21_re"]["current"]).astype(float)
    f = np.ravel(d["s21_re"]["f_vna"]).astype(float)
    z = np.ravel(d["s21_re"]["s21_re"]) + 1j * np.ravel(d["s21_im"]["s21_im"])
    # bloques consecutivos de la misma corriente (en el orden en que se midieron)
    cortes = np.r_[0, np.flatnonzero(np.diff(I) != 0) + 1, len(I)]
    n_f = cortes[1] - cortes[0]               # descarta la última traza si quedó incompleta
    completos = [(a, b) for a, b in zip(cortes[:-1], cortes[1:]) if b - a == n_f]
    I_u = np.array([I[a] for a, _ in completos])
    Z = np.array([z[a:b] for a, b in completos])
    return dict(run_id=run_id, I_uA=I_u, f=f[:n_f], Z=Z, metadata=dict(ds.metadata), nombre=ds.name)


def _ajustar_trazas(f, Z, fcrop="auto", n_anchos=4):
    from fit_resonancia import ajustar_notch
    claves = ("fr", "fr_err", "Ql", "Qi", "Qi_err", "Qc", "Qc_err", "kappa", "caida_db", "f_dip")
    out = {k: np.full(len(Z), np.nan) for k in claves}
    out["ok"] = np.zeros(len(Z), bool)
    for i, z in enumerate(Z):
        try:
            r = ajustar_notch(f, z, fcrop=fcrop, n_anchos=n_anchos)
        except Exception:
            continue
        for k in claves:
            out[k][i] = r[k]
        out["ok"][i] = r["ok"]
    return out


def ajustar_barrido_corriente(datos, fcrop="auto", n_anchos=4):
    """Ajuste notch de cada traza del run 2D. Devuelve dict de arrays (uno por corriente)."""
    out = _ajustar_trazas(datos["f"], datos["Z"], fcrop, n_anchos)
    out["I_uA"] = datos["I_uA"]
    return out


def ajustar_runs(runs, fcrop="auto", n_anchos=4):
    """Ajuste notch de una lista de runs 1D (p.ej. barrido_potencia(...)['runs'])."""
    from qcodes.dataset import load_by_id
    f0, Z, P = None, [], []
    for x in runs:
        ds = load_by_id(x["run_id"])
        d = ds.get_parameter_data()
        f = np.ravel(d["s21_re"]["f_vna"]).astype(float)
        z = np.ravel(d["s21_re"]["s21_re"]) + 1j * np.ravel(d["s21_im"]["s21_im"])
        f0 = f if f0 is None else f0
        Z.append(np.interp(f0, f, z.real) + 1j * np.interp(f0, f, z.imag) if len(f) != len(f0) else z)
        P.append(float(ds.metadata.get("pow_vna", np.nan)))
    out = _ajustar_trazas(f0, np.array(Z), fcrop, n_anchos)
    out["P"] = np.array(P)
    out["f"], out["Z"] = f0, np.array(Z)
    return out


def saltos(I_uA, fr, k=6, umbral_Hz=None, n_sigma=3.0, minimo_Hz=30e3, sep=None,
           signo="dominante", periodico=True, n_sigma_periodico=2.5, n_ajuste=8, min_signif=2.0):
    """Detecta saltos de fr vs corriente.

    Para cada posición j ajusta, en una ventana de k puntos a cada lado, una recta común más
    un escalón en j:   fr = a + b (I - I_j) + c · [I >= I_j]
    c es el tamaño del salto. Al tener la recta en el modelo, la pendiente dentro de cada rama
    del diente de sierra no se confunde con un salto (con datos poco ruidosos, como a potencia
    alta, la diferencia de medias sí la confundía).
    Cerca de los extremos del barrido la ventana se achica (mínimo 4 puntos por lado).
    Se aceptan los máximos locales de |c| por encima del umbral, separados más de `sep` puntos.

    umbral_Hz : umbral fijo [Hz]; None = max(n_sigma × error de c, minimo_Hz), con el error de c
                calculado a partir del ruido punto a punto, que se estima robustamente con la MAD
                de las segundas diferencias (no la afectan ni la pendiente ni la curvatura suave).
    sep       : separación mínima entre saltos [puntos]; None = k.
    signo     : 'dominante' = se queda solo con los saltos del signo más frecuente;
                '+' / '-' = solo ese signo (en el orden de medición); None = todos.
    periodico : si hay >= 3 saltos, ajusta su período y busca saltos que falten en las posiciones
                esperadas (± período/6) con un umbral más bajo (n_sigma_periodico), con el mismo
                signo y una ventana más ancha (0.4 períodos por lado, sin cruzar otros saltos).
                Sirve para los saltos chicos en datos ruidosos.
    n_ajuste  : puntos a cada lado para medir el tamaño del salto con dos rectas locales
                (evaluadas en la posición del salto, sin cruzar otros saltos).
    min_signif: descarta saltos cuyo tamaño ajustado sea < min_signif × su error (o de signo
                contrario); None = no filtra.

    Devuelve una lista de dicts (en el orden de medición) con
    j (el salto ocurre entre los índices j-1 y j), I (posición [uA]), df y df_err [Hz] (tamaño
    del salto en el orden de medición, de las rectas locales) y score [Hz] (c del ajuste con escalón).
    """
    I_uA, fr = np.asarray(I_uA, float), np.asarray(fr, float)
    idx = np.flatnonzero(np.isfinite(fr) & np.isfinite(I_uA))      # ignora ajustes fallidos
    I, y = I_uA[idx], fr[idx]
    n = len(y)
    sep = k if sep is None else sep
    if n < 2 * k + 1:
        return []
    d2 = y[2:] - 2 * y[1:-1] + y[:-2]
    sig = 1.4826 * np.median(np.abs(d2 - np.median(d2))) / np.sqrt(6)    # ruido de un punto

    def escalon(j, kl, kr):
        """Ajuste recta + escalón en j con kl puntos a la izquierda y kr a la derecha: (c, error)."""
        x = I[j - kl:j + kr] - (I[j - 1] + I[j]) / 2
        X = np.column_stack([np.ones(kl + kr), x, np.r_[np.zeros(kl), np.ones(kr)]])
        coef, *_ = np.linalg.lstsq(X, y[j - kl:j + kr], rcond=None)
        return coef[2], sig * np.sqrt(np.linalg.inv(X.T @ X)[2, 2])

    sc, err_c = np.zeros(n), np.full(n, np.inf)
    for j in range(4, n - 3):
        sc[j], err_c[j] = escalon(j, min(k, j), min(k, n - j))
    umbral = np.maximum(n_sigma * err_c, minimo_Hz) if umbral_Hz is None else np.full(n, umbral_Hz)
    elegidos = []
    for j in np.argsort(-np.abs(sc)):
        if abs(sc[j]) == 0:
            break
        if abs(sc[j]) < umbral[j]:
            continue
        if all(abs(j - e) > sep for e in elegidos):
            elegidos.append(int(j))
    elegidos.sort()
    if signo is not None and elegidos:
        s = np.sign(sc[elegidos])
        quiero = (1 if s.sum() >= 0 else -1) if signo == "dominante" else (1 if signo == "+" else -1)
        elegidos = [j for j in elegidos if np.sign(sc[j]) == quiero]

    if periodico and len(elegidos) >= 3:
        Im = (I[:-1] + I[1:]) / 2                      # Im[j-1] = posición de un salto en j
        pos = Im[np.array(elegidos) - 1]
        o = np.argsort(pos); pos = pos[o]
        T0 = np.median(np.diff(pos))
        m_ = np.round((pos - pos[0]) / T0)
        T, x0 = np.polyfit(m_, pos, 1)
        sgn = np.sign(np.sum(np.sign(sc[elegidos])))
        lo, hi = min(Im[3], Im[n - 5]), max(Im[3], Im[n - 5])
        k2 = max(k, int(round(0.4 * abs(T) / np.median(np.abs(np.diff(I))))))
        for m in range(int(np.floor((lo - x0) / T)) - 1, int(np.ceil((hi - x0) / T)) + 2):
            esperado = x0 + m * T
            if not lo <= esperado <= hi or np.min(np.abs(pos - esperado)) < T / 4:
                continue
            cand = []
            for j in range(4, n - 3):
                if abs(Im[j - 1] - esperado) > abs(T) / 6 or not all(abs(j - e) > sep for e in elegidos):
                    continue
                izq = max([e for e in elegidos if e < j], default=0)
                der = min([e for e in elegidos if e > j], default=n)
                c, e_c = escalon(j, min(k2, j - izq), min(k2, der - j))
                if np.sign(c) == sgn and abs(c) >= max(n_sigma_periodico * e_c, minimo_Hz):
                    cand.append((abs(c) / e_c, j, c))
            if cand:
                _, j, c = max(cand)
                sc[j] = c
                elegidos.append(j)
        elegidos.sort()

    out = []
    for m, j in enumerate(elegidos):
        Ij = (I[j - 1] + I[j]) / 2
        a = max(j - n_ajuste, elegidos[m - 1] if m > 0 else 0)
        b = min(j + n_ajuste, elegidos[m + 1] if m + 1 < len(elegidos) else n)
        df, df_err = sc[j], np.nan
        if j - a >= 3 and b - j >= 3:
            L, cL = np.polyfit(I[a:j], y[a:j], 1, cov=True)
            R, cR = np.polyfit(I[j:b], y[j:b], 1, cov=True)
            v = np.array([Ij, 1.0])
            df = np.polyval(R, Ij) - np.polyval(L, Ij)
            df_err = np.sqrt(v @ cL @ v + v @ cR @ v)
        if min_signif is not None and np.isfinite(df_err) and df * np.sign(sc[j]) < min_signif * df_err:
            continue
        out.append(dict(j=int(idx[j]), I=float(Ij), df=float(df), df_err=float(df_err), score=float(sc[j])))
    return out
