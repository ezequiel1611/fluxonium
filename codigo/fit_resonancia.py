r"""
Ajuste de una resonancia notch (resonador acoplado a una línea) desde una base QCoDeS.
No necesita repos externos: solo numpy, scipy, pandas, matplotlib y qcodes.

Modelo (Probst et al., Rev. Sci. Instrum. 86, 024706 (2015)) + fondo lineal complejo:

    S21(f) = (A0 + A1·x) · e^{-2πiτ(f - fm)} · [ 1 - (Ql/|Qc|) e^{iφ} / (1 + 2iQl (f - fr)/fr) ]

    1/Qi = 1/Ql - cos(φ)/|Qc|        Qc = |Qc| / cos(φ)        κ/2π = fr/Ql

Número medio de fotones en el resonador:

    <n> = 2 P_in Ql² / (Qc ħ ω0²),    ω0 = 2π fr,    P_in [W] = 1e-3 · 10^((P_VNA + att_total)/10)
    (resonador "hanger" acoplado a una línea de paso: Probst 2015, Bruno 2015)

Uso:
    res = ajustar_resonancia(r"C:\...\fluxonium_database\resonadores", 30, att_total=-112)

    # run con varias corrientes / frecuencias de SG: elegir una traza
    res = ajustar_resonancia(db, 129, corte={'current': 50}, att_total=-112)

Formatos de run que reconoce:
    - barrido simple del PNA:  pna_tr1_magnitude / pna_tr2_(unwrapped_)phase vs pna_frequency_axis
                               (+ 'pow' guardado como parámetro)
    - custom:                  s21_mag / s21_phase vs f_vna   (o s21_mag_res / s21_phase_res)
                               con o sin una segunda variable (current, f_sg, ...) -> usar `corte`
La potencia del VNA se toma de: argumento `potencia`, parámetro 'pow', o metadata 'pow_vna'/'pow_cw'.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from scipy.constants import hbar


# =============================================================================
# 1) Carga desde la base de datos
# =============================================================================
def cargar_traza(db, run_id, corte=None):
    """Devuelve dict con freq [Hz], z (S21 complejo), potencia [dBm o None] y metadata.

    corte: dict {variable: valor} para elegir UNA traza en runs 2D,
           p.ej. {'current': 50} (µA) o {'f_sg': 13.6} (GHz). Usa el valor medido más cercano.
    """
    from qcodes.dataset import load_by_id
    from qcodes.dataset.sqlite.database import connect

    if not os.path.isfile(db):
        raise FileNotFoundError(f'no existe la base: {db}')
    ds = load_by_id(run_id, conn=connect(db))
    data = ds.get_parameter_data()

    # --- magnitud: el parámetro cuyo eje es una frecuencia del VNA ---
    def es_eje_frecuencia(k):
        return k in ('f_vna', 'pna_frequency_axis') or 'frequency_axis' in k

    cand = [k for k in data if ('mag' in k.lower()) and any(es_eje_frecuencia(s) for s in data[k])]
    if not cand:
        raise ValueError(f'run {run_id}: no encontré magnitud vs frecuencia del VNA. '
                         f'Parámetros: {list(data)}')
    k_mag = sorted(cand, key=lambda k: ('_res' not in k, k))[0]       # prioriza *_res
    k_pha = next((k for k in data if 'phase' in k.lower()
                  and k.replace('phase', 'mag') == k_mag.replace('magnitude', 'mag')), None) \
        or next(k for k in data if 'phase' in k.lower()
                and any(es_eje_frecuencia(s) for s in data[k]))
    k_f = next(s for s in data[k_mag] if es_eje_frecuencia(s))
    otros = [s for s in data[k_mag] if s not in (k_mag, k_f)]

    f = np.ravel(data[k_mag][k_f]).astype(float)
    mag = np.ravel(data[k_mag][k_mag]).astype(float)
    pha = np.ravel(data[k_pha][k_pha]).astype(float)

    # --- elegir una traza si el run es 2D ---
    elegido = {}
    if otros:
        if not corte:
            vals = {o: np.unique(np.ravel(data[k_mag][o])) for o in otros}
            raise ValueError(f'run {run_id} tiene varias trazas ({", ".join(otros)}). '
                             f'Pasá corte={{"{otros[0]}": valor}}. '
                             f'Rango: {vals[otros[0]].min():g} a {vals[otros[0]].max():g}')
        m = np.ones_like(f, bool)
        for o in otros:
            if o not in corte:
                raise ValueError(f'falta {o!r} en corte (variables: {otros})')
            v = np.ravel(data[k_mag][o]).astype(float)
            u = np.unique(v)
            real = u[np.argmin(np.abs(u - corte[o]))]
            m &= v == real
            elegido[o] = real
        f, mag, pha = f[m], mag[m], pha[m]

    # --- unidades ---
    if np.nanmax(f) < 1e6:                     # guardado en GHz
        f = f * 1e9
    unidad = (ds.paramspecs[k_pha].unit or 'deg').lower()
    pha_rad = pha if unidad.startswith('rad') else np.deg2rad(pha)
    o = np.argsort(f)
    f, mag, pha_rad = f[o], mag[o], np.unwrap(pha_rad[o])
    z = 10 ** (mag / 20) * np.exp(1j * pha_rad)

    # --- potencia del VNA ---
    pot = None
    if 'pow' in data:
        pot = float(np.ravel(data['pow']['pow'])[0])
    else:
        for k in ('pow_vna', 'pow_cw', 'power_vna'):
            if k in ds.metadata:
                pot = float(ds.metadata[k]); break

    return dict(run_id=run_id, freq=f, z=z, mag_db=mag, potencia=pot, corte=elegido,
                metadata=dict(ds.metadata), nombre=ds.name)


# =============================================================================
# 2) Ajuste
# =============================================================================
def _estimacion_inicial(f, z):
    """Delay (pendiente de fase en los bordes), fr, FWHM y profundidad del dip."""
    fm = f.mean()
    b = max(3, len(f) // 5)
    e = np.r_[:b, len(f) - b:len(f)]
    ph = np.unwrap(np.angle(z))
    tau = -np.polyfit(f[e] - fm, ph[e], 1)[0] / (2 * np.pi)
    zc = z * np.exp(2j * np.pi * tau * (f - fm))
    zn = np.abs(zc / np.mean(zc[e]))
    zs = np.convolve(zn, np.ones(3) / 3, 'same')
    i0 = np.argmin(zs[1:-1]) + 1
    d = min(zn[i0], 0.95)
    medio = np.sqrt((1 + d ** 2) / 2)
    l = r = i0
    while l > 0 and zn[l] < medio:
        l -= 1
    while r < len(f) - 1 and zn[r] < medio:
        r += 1
    fwhm = max(f[r] - f[l], 3 * np.median(np.diff(f)))
    return tau, f[i0], fwhm, d


def ajustar_notch(f, z, fcrop='auto', n_anchos=4):
    """Ajusta el modelo notch a S21 complejo. Devuelve dict con parámetros y errores."""
    tau, fr, fwhm, d = _estimacion_inicial(f, z)
    if isinstance(fcrop, str) and fcrop == 'auto':
        fcrop = (fr - n_anchos * fwhm, fr + n_anchos * fwhm)
    if fcrop is not None:
        m = (f >= fcrop[0]) & (f <= fcrop[1])
        f, z = f[m], z[m]
    if len(f) < 10:
        raise ValueError('muy pocos puntos en la ventana de ajuste')

    Ql, Qc = fr / fwhm, fr / fwhm / (1 - d)
    fm, span = f.mean(), f[-1] - f[0]
    x, w0 = (f - fm) / span, fr / Ql

    def forma(fr_, ql, qc, p):
        return 1 - (ql / qc) * np.exp(1j * p) / (1 + 2j * ql * (f - fr_) / fr_)

    g = forma(fr, Ql, Qc, 0) * np.exp(-2j * np.pi * tau * (f - fm))
    (A0, A1), *_ = np.linalg.lstsq(np.c_[g, g * x], z, rcond=None)
    s = abs(A0)

    def unpack(u):
        return (fr + w0 * u[0], np.exp(u[1]), np.exp(u[2]), u[3],
                tau + u[4] / (2 * np.pi * span),
                s * (u[5] + 1j * u[6]), s * (u[7] + 1j * u[8]))

    def modelo(u):
        fr_, ql, qc, p, t, a0, a1 = unpack(u)
        return (a0 + a1 * x) * np.exp(-2j * np.pi * t * (f - fm)) * forma(fr_, ql, qc, p)

    def resid(u):
        dz = (modelo(u) - z) / s
        return np.r_[dz.real, dz.imag]

    u0 = [0, np.log(Ql), np.log(Qc), 0, 0, A0.real / s, A0.imag / s, A1.real / s, A1.imag / s]
    with np.errstate(all='ignore'):
        sol = least_squares(resid, u0, method='lm', max_nfev=5000)
    u = sol.x
    s2 = 2 * sol.cost / max(2 * len(f) - len(u), 1)
    try:
        cov = np.linalg.inv(sol.jac.T @ sol.jac) * s2
    except np.linalg.LinAlgError:
        cov = np.full((len(u), len(u)), np.nan)

    fr_, ql, qc, p, t, a0, a1 = unpack(u)
    qi = 1 / (1 / ql - np.cos(p) / qc)
    g_qi = np.zeros(len(u))
    g_qi[1], g_qi[2], g_qi[3] = qi ** 2 / ql, -qi ** 2 * np.cos(p) / qc, -qi ** 2 * np.sin(p) / qc
    Qc_corr = qc / np.cos(p)
    g_qc = np.zeros(len(u))
    g_qc[2], g_qc[3] = Qc_corr, Qc_corr * np.tan(p)

    zs = modelo(u)
    ruido = np.std(np.diff(z, 2)) / np.sqrt(6)
    entorno = (a0 + a1 * x) * np.exp(-2j * np.pi * t * (f - fm))

    # profundidad del dip (modelo normalizado)
    ff = np.linspace(f[0], f[-1], 20001)
    zz = np.abs(1 - ql / qc * np.exp(1j * p) / (1 + 2j * ql * (ff - fr_) / fr_))

    return dict(
        fr=fr_, fr_err=w0 * np.sqrt(cov[0, 0]),
        Ql=ql, Ql_err=ql * np.sqrt(cov[1, 1]),
        Qi=qi, Qi_err=np.sqrt(g_qi @ cov @ g_qi),
        Qc=Qc_corr, Qc_err=np.sqrt(g_qc @ cov @ g_qc),
        absQc=qc, phi0=p, delay=t,
        kappa=fr_ / ql, kappa_i=fr_ / qi, kappa_c=fr_ / Qc_corr,
        f_dip=ff[np.argmin(zz)], caida_db=-20 * np.log10(zz.min()),
        res_noise=np.sqrt(np.mean(np.abs(z - zs) ** 2)) / ruido if ruido > 0 else np.nan,
        ok=bool(sol.success and ql > 0 and qi > 0 and f[0] <= fr_ <= f[-1]),
        _f=f, _z=z, _zs=zs, _zn=z / entorno, _zsn=zs / entorno,
    )


# =============================================================================
# 3) Número de fotones
# =============================================================================
def numero_fotones(P_vna_dBm, att_total_dB, fr, Ql, Qc):
    """<n> = 2 P_in Ql² / (Qc ħ ω0²).  P_in = P_VNA + att_total (dBm) -> W.

    Factor 2 (no 4): en un hanger la señal entra por un lado y el resonador pierde hacia
    los dos, así que el acoplamiento de entrada es κc/2  ->  <n> = 2 κc P / (ħ ω κ²).
    """
    P_in_W = 1e-3 * 10 ** ((P_vna_dBm + att_total_dB) / 10)
    return 2 * P_in_W * Ql ** 2 / (Qc * hbar * (2 * np.pi * fr) ** 2)


# =============================================================================
# 4) Función principal
# =============================================================================
def ajustar_resonancia(db, run_id, att_total=None, potencia=None, corte=None,
                       fcrop='auto', n_anchos=4, graficar=True, imprimir=True):
    """Carga el run `run_id` de la base `db`, ajusta la resonancia y calcula <n>.

    att_total: atenuación total hasta el resonador [dB, negativa], p.ej. -112.
               Si es None no se calcula <n>.
    potencia:  potencia del VNA [dBm]; None = la guardada en el run.
    corte:     {'current': valor} o {'f_sg': valor} si el run tiene varias trazas.
    fcrop:     'auto' (dip ± n_anchos·FWHM), None (todo) o (f1, f2) en Hz.
    """
    d = cargar_traza(db, run_id, corte)
    r = ajustar_notch(d['freq'], d['z'], fcrop=fcrop, n_anchos=n_anchos)

    P = potencia if potencia is not None else d['potencia']
    r.update(run_id=run_id, corte=d['corte'], P_vna=P, att_total=att_total)
    if att_total is not None and P is not None:
        r['P_in_dBm'] = P + att_total
        r['n_fotones'] = numero_fotones(P, att_total, r['fr'], r['Ql'], r['Qc'])
        r['P_un_foton_dBm'] = r['P_in_dBm'] - 10 * np.log10(r['n_fotones'])   # P_in con <n> = 1
    elif att_total is not None:
        print('⚠ no encontré la potencia del VNA en el run: pasala con potencia=...')

    if imprimir:
        imprimir_resultado(r)
    if graficar:
        graficar_ajuste(r)
    return {k: v for k, v in r.items() if not k.startswith('_')}


def imprimir_resultado(r):
    corte = ''.join(f', {k} = {v:g}' for k, v in r['corte'].items())
    print(f"run {r['run_id']}{corte}")
    print(f"  fr     = {r['fr']/1e9:.7f} GHz ± {r['fr_err']/1e3:.2f} kHz")
    print(f"  κ/2π   = {r['kappa']/1e3:.2f} kHz   (κi/2π = {r['kappa_i']/1e3:.2f}, κc/2π = {r['kappa_c']/1e3:.2f} kHz)")
    print(f"  Ql     = {r['Ql']:.4g} ± {r['Ql_err']:.2g}")
    print(f"  Qi     = {r['Qi']:.4g} ± {r['Qi_err']:.2g}")
    print(f"  Qc     = {r['Qc']:.4g} ± {r['Qc_err']:.2g}   (|Qc| = {r['absQc']:.4g}, φ0 = {r['phi0']:.3f} rad)")
    print(f"  caída  = {r['caida_db']:.2f} dB,  delay = {r['delay']*1e9:.2f} ns")
    print(f"  res/noise = {r['res_noise']:.2f}  {'(ok)' if r['ok'] else '⚠ ajuste dudoso'}")
    if 'n_fotones' in r:
        print(f"  P_VNA = {r['P_vna']:.1f} dBm, att = {r['att_total']:.1f} dB -> P_in = {r['P_in_dBm']:.1f} dBm")
        print(f"  <n>    = {r['n_fotones']:.3g} fotones   (<n> = 1 en P_in ≈ {r['P_un_foton_dBm']:.1f} dBm)")


def graficar_ajuste(r):
    f = r['_f'] / 1e9
    z, zs, zn, zsn = r['_z'], r['_zs'], r['_zn'], r['_zsn']
    corr = np.exp(2j * np.pi * r['_f'] * r['delay'])

    fig = plt.figure(figsize=(11, 6))
    gs = fig.add_gridspec(2, 2, width_ratios=[1.4, 1])
    ax_m = fig.add_subplot(gs[0, 0])
    ax_p = fig.add_subplot(gs[1, 0], sharex=ax_m)
    ax_c = fig.add_subplot(gs[:, 1])

    ax_m.plot(f, 20 * np.log10(np.abs(z)), '.', ms=3, label='datos')
    ax_m.plot(f, 20 * np.log10(np.abs(zs)), '-', lw=1.5, label='ajuste')
    ax_m.axvline(r['fr'] / 1e9, color='gray', ls=':', lw=1, label='fr')
    ax_m.set_ylabel('|S21| [dB]'); ax_m.legend()

    ax_p.plot(f, np.unwrap(np.angle(z * corr)), '.', ms=3)
    ax_p.plot(f, np.unwrap(np.angle(zs * corr)), '-', lw=1.5)
    ax_p.axvline(r['fr'] / 1e9, color='gray', ls=':', lw=1)
    ax_p.set_ylabel('fase − 2πfτ [rad]'); ax_p.set_xlabel('frecuencia [GHz]')

    ax_c.plot(zn.real, zn.imag, '.', ms=3, label='datos norm.')
    ax_c.plot(zsn.real, zsn.imag, '-', lw=1.5, label='ajuste')
    ax_c.plot([1], [0], 'k+', ms=10)
    ax_c.set_aspect('equal', 'datalim')
    ax_c.set_xlabel('Re S21 (norm.)'); ax_c.set_ylabel('Im S21 (norm.)')
    ax_c.legend(loc='upper left')

    txt = (f"fr = {r['fr']/1e9:.6f} GHz\n"
           f"κ/2π = {r['kappa']/1e3:.1f} kHz\n"
           f"Ql = {r['Ql']:.3g}\nQi = {r['Qi']:.3g} ± {r['Qi_err']:.2g}\n"
           f"Qc = {r['Qc']:.3g}\nres/noise = {r['res_noise']:.2f}")
    if 'n_fotones' in r:
        txt += f"\n⟨n⟩ = {r['n_fotones']:.3g}"
    ax_c.text(0.02, 0.02, txt, transform=ax_c.transAxes, fontsize=9, va='bottom',
              family='monospace', bbox=dict(boxstyle='round', fc='white', alpha=0.85))

    corte = ''.join(f', {k} = {v:g}' for k, v in r['corte'].items())
    fig.suptitle(f"run {r['run_id']}{corte}")
    fig.tight_layout()
    plt.show()
    return fig