r"""
Número de fotones en un resonador a partir de UN run de QCoDeS, usando fit_resonancia.py.

    <n> = 4 P_in Ql² / (Qc ħ ω0²),   ω0 = 2π fr,   P_in [W] = 1e-3 · 10^((P_VNA + att_ext + att_cryo)/10)

Por defecto att_ext = -62 dB (los dos atenuadores afuera) y att_cryo = -50 dB (adentro del criostato):
con P_VNA = 0 dBm -> P_in = -112 dBm = 6.3e-15 W.

"Puede que no sea lineal": a potencias altas el resonador se vuelve no lineal (Kerr): el dip se
inclina hacia un lado y se ensancha, y la Lorentziana ajusta mal. Por eso se ajustan DOS modelos:

  lineal     S21 = entorno · [1 - (Ql/|Qc|) e^{iφ} / (1 + 2i y0)],     y0 = Ql (f - fr)/fr
  no lineal  igual pero con y en lugar de y0, donde y resuelve  y = y0 + a / (1 + 4y²)
             (Swenson et al., J. Appl. Phys. 113, 104501 (2013)).
             a = 0 es lineal; a > 0 inclina el dip hacia frecuencias bajas;
             |a| > 4/(3√3) ≈ 0.77 -> bifurcación (el dip "salta").

Con modelo='auto' se usa el no lineal solo si mejora el ajuste de forma significativa.
En ese caso fr es la frecuencia a baja potencia (sin el corrimiento Kerr), y Ql, Qc son los
lineales, que es lo que pide la fórmula de <n>.

Uso:
    import fotones_resonador as fot
    r = fot.fotones(30, db=r"C:\...\fluxonium_database\resonadores")          # potencia del run
    r = fot.fotones(30, db=..., potencia=0)                                     # forzar 0 dBm
    r = fot.fotones(129, db=..., corte={'current': 50})                         # run 2D
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

import fit_resonancia as frs            # cargar_traza, ajustar_notch, numero_fotones

A_BIFURCACION = 4 / (3 * np.sqrt(3))    # ≈ 0.77


# =============================================================================
# Modelo no lineal (Kerr)
# =============================================================================
def _y_kerr(y0, a):
    """Resuelve y = y0 + a/(1+4y²)  <=>  y³ - y0 y² + y/4 - (y0+a)/4 = 0  para cada punto.
    Si hay 3 raíces reales (bifurcación) toma la rama de un barrido hacia frecuencias crecientes."""
    n = len(y0)
    C = np.zeros((n, 3, 3))
    C[:, 0, 0] = y0                      # -c2,  c2 = -y0
    C[:, 0, 1] = -0.25                   # -c1,  c1 = 1/4
    C[:, 0, 2] = (y0 + a) / 4            # -c0,  c0 = -(y0+a)/4
    C[:, 1, 0] = 1.0
    C[:, 2, 1] = 1.0
    raices = np.linalg.eigvals(C)
    reales = np.where(np.abs(raices.imag) < 1e-7 * (1 + np.abs(raices.real)), raices.real, np.nan)
    return np.nanmax(reales, axis=1) if a >= 0 else np.nanmin(reales, axis=1)


def ajustar_notch_kerr(f, z, lin):
    """Ajuste no lineal arrancando del ajuste lineal `lin` (salida de frs.ajustar_notch).
    Usa la misma ventana de frecuencias que el lineal."""
    fm, span = f.mean(), f[-1] - f[0]
    x = (f - fm) / span
    fr0, Ql0, Qc0, p0, tau0 = lin['fr'], lin['Ql'], lin['absQc'], lin['phi0'], lin['delay']
    w0 = fr0 / Ql0

    def forma(fr_, ql, qc, p, a):
        y = _y_kerr(ql * (f - fr_) / fr_, a)
        return 1 - (ql / qc) * np.exp(1j * p) / (1 + 2j * y)

    g = forma(fr0, Ql0, Qc0, p0, 0.0) * np.exp(-2j * np.pi * tau0 * (f - fm))
    (A0, A1), *_ = np.linalg.lstsq(np.c_[g, g * x], z, rcond=None)
    s = abs(A0)

    def unpack(u):
        return (fr0 + w0 * u[0], Ql0 * np.exp(u[1]), Qc0 * np.exp(u[2]), u[3], u[4],
                tau0 + u[5] / (2 * np.pi * span),
                s * (u[6] + 1j * u[7]), s * (u[8] + 1j * u[9]))

    def modelo(u):
        fr_, ql, qc, p, a, t, a0, a1 = unpack(u)
        return (a0 + a1 * x) * np.exp(-2j * np.pi * t * (f - fm)) * forma(fr_, ql, qc, p, a)

    def resid(u):
        r = (modelo(u) - z) / s
        r = np.r_[r.real, r.imag]
        return np.where(np.isfinite(r), r, 1e3)

    mejor = None
    for a_ini in (0.05, 0.3, -0.05, -0.3, 0.7):          # varios arranques: el cúbico tiene ramas
        u0 = [0, 0, 0, p0, a_ini, 0, A0.real / s, A0.imag / s, A1.real / s, A1.imag / s]
        with np.errstate(all='ignore'):
            sol = least_squares(resid, u0, method='lm', max_nfev=4000)
        if mejor is None or sol.cost < mejor.cost:
            mejor = sol
    u = mejor.x
    s2 = 2 * mejor.cost / max(2 * len(f) - len(u), 1)
    try:
        cov = np.linalg.inv(mejor.jac.T @ mejor.jac) * s2
    except np.linalg.LinAlgError:
        cov = np.full((len(u), len(u)), np.nan)

    fr_, ql, qc, p, a, t, a0, a1 = unpack(u)
    qi = 1 / (1 / ql - np.cos(p) / qc)
    zs = modelo(u)
    ruido = np.std(np.diff(z, 2)) / np.sqrt(6)
    return dict(fr=fr_, fr_err=w0 * np.sqrt(cov[0, 0]), Ql=ql, Ql_err=ql * np.sqrt(cov[1, 1]),
                Qi=qi, Qc=qc / np.cos(p), absQc=qc, phi0=p, a=a, a_err=np.sqrt(cov[4, 4]),
                delay=t, kappa=fr_ / ql, chi2=s2,
                res_noise=np.sqrt(np.mean(np.abs(z - zs) ** 2)) / ruido if ruido > 0 else np.nan,
                ok=bool(mejor.success and ql > 0 and qi > 0), _zs=zs)


# =============================================================================
# Función principal
# =============================================================================
def fotones(run_id, db=None, potencia=None, att_ext=-62.0, att_cryo=-50.0, corte=None,
            modelo='auto', fcrop='auto', n_anchos=6, graficar=True):
    """Ajusta la resonancia del run `run_id` y calcula <n>.

    db:        path de la base (None = la configurada en qc.config)
    potencia:  P_VNA [dBm] a la entrada; None = la guardada en el run
    att_ext:   atenuación fuera del criostato [dB]   (default -62)
    att_cryo:  atenuación dentro del criostato [dB]  (default -50)
    corte:     {'current': valor} / {'f_sg': valor} para runs con varias trazas
    modelo:    'auto' | 'lineal' | 'kerr'
    """
    import qcodes as qc
    db = db or qc.config['core']['db_location']

    d = frs.cargar_traza(db, run_id, corte)
    P = potencia if potencia is not None else d['potencia']
    if P is None:
        raise ValueError('no encontré la potencia del VNA en el run: pasala con potencia=...')

    # ---- ajuste lineal (fit_resonancia) ----
    lin = frs.ajustar_notch(d['freq'], d['z'], fcrop=fcrop, n_anchos=n_anchos)
    f, z = lin['_f'], lin['_z']
    chi2_lin = np.mean(np.abs(z - lin['_zs']) ** 2)

    # ---- ajuste no lineal (Kerr) ----
    kerr = None
    if modelo in ('auto', 'kerr'):
        try:
            kerr = ajustar_notch_kerr(f, z, lin)
        except Exception as e:  # noqa: BLE001
            print('ajuste no lineal falló:', e)

    usar_kerr = False
    if kerr is not None:
        chi2_k = np.mean(np.abs(z - kerr['_zs']) ** 2)
        mejora = chi2_lin / chi2_k if chi2_k > 0 else 1
        significativo = abs(kerr['a']) > 3 * (kerr['a_err'] if np.isfinite(kerr['a_err']) else np.inf)
        usar_kerr = modelo == 'kerr' or (kerr['ok'] and mejora > 1.2 and significativo)
        kerr['mejora_chi2'] = mejora
    elegido = kerr if usar_kerr else lin

    # ---- número de fotones ----
    att = att_ext + att_cryo
    P_in_dBm = P + att
    P_in_W = 1e-3 * 10 ** (P_in_dBm / 10)
    n = frs.numero_fotones(P, att, elegido['fr'], elegido['Ql'], elegido['Qc'])

    r = dict(run_id=run_id, corte=d['corte'], modelo='kerr' if usar_kerr else 'lineal',
             P_vna_dBm=P, att_ext=att_ext, att_cryo=att_cryo, P_in_dBm=P_in_dBm, P_in_W=P_in_W,
             fr=elegido['fr'], Ql=elegido['Ql'], Qi=elegido['Qi'], Qc=elegido['Qc'],
             kappa=elegido['fr'] / elegido['Ql'], res_noise=elegido['res_noise'],
             n_fotones=n, P_in_un_foton_dBm=P_in_dBm - 10 * np.log10(n),
             P_vna_un_foton_dBm=P - 10 * np.log10(n),
             n_lineal=frs.numero_fotones(P, att, lin['fr'], lin['Ql'], lin['Qc']))
    if kerr is not None:
        r.update(a_kerr=kerr['a'], a_kerr_err=kerr['a_err'], mejora_chi2=kerr['mejora_chi2'],
                 n_kerr=frs.numero_fotones(P, att, kerr['fr'], kerr['Ql'], kerr['Qc']))

    _imprimir(r, lin, kerr)
    if graficar:
        _graficar(r, f, z, lin, kerr)
    return r


def _imprimir(r, lin, kerr):
    corte = ''.join(f', {k} = {v:g}' for k, v in r['corte'].items())
    print(f"run {r['run_id']}{corte}")
    print(f"  P_VNA = {r['P_vna_dBm']:.1f} dBm, att = {r['att_ext']:.0f} + {r['att_cryo']:.0f} dB"
          f"  ->  P_in = {r['P_in_dBm']:.1f} dBm = {r['P_in_W']:.3g} W")
    print(f"  lineal : fr = {lin['fr']/1e9:.7f} GHz, Ql = {lin['Ql']:.4g}, Qc = {lin['Qc']:.4g}, "
          f"Qi = {lin['Qi']:.4g}, res/noise = {lin['res_noise']:.2f}  ->  <n> = {r['n_lineal']:.3g}")
    if kerr is not None:
        print(f"  Kerr   : fr = {kerr['fr']/1e9:.7f} GHz, Ql = {kerr['Ql']:.4g}, Qc = {kerr['Qc']:.4g}, "
              f"Qi = {kerr['Qi']:.4g}, a = {kerr['a']:.3f} ± {kerr['a_err']:.2g}, "
              f"res/noise = {kerr['res_noise']:.2f}  ->  <n> = {r['n_kerr']:.3g}")
    print(f"  modelo usado: {r['modelo']}")
    print(f"  <n> = {r['n_fotones']:.3g} fotones   "
          f"(<n> = 1 en P_in ≈ {r['P_in_un_foton_dBm']:.1f} dBm, P_VNA ≈ {r['P_vna_un_foton_dBm']:.1f} dBm)")
    if kerr is not None and abs(kerr['a']) > A_BIFURCACION:
        print(f"  ⚠ |a| > {A_BIFURCACION:.2f}: el resonador está bifurcado; <n> es solo una estimación")
    elif r['modelo'] == 'kerr':
        print('  ⚠ el resonador es no lineal a esta potencia: <n> usa fr, Ql, Qc del modelo Kerr')
    if r['res_noise'] > 3:
        print('  ⚠ los residuos son grandes: revisar ventana (fcrop) o la forma del dip')


def _graficar(r, f, z, lin, kerr):
    fg = f / 1e9
    corr = np.exp(2j * np.pi * f * lin['delay'])
    fig, ax = plt.subplots(2, 1, figsize=(9, 6.5), sharex=True)
    ax[0].plot(fg, 20 * np.log10(np.abs(z)), '.', ms=3, color='k', label='datos')
    ax[0].plot(fg, 20 * np.log10(np.abs(lin['_zs'])), '-', lw=1.4,
               label=f"lineal (res/noise {lin['res_noise']:.2f})")
    ax[1].plot(fg, np.unwrap(np.angle(z * corr)), '.', ms=3, color='k')
    ax[1].plot(fg, np.unwrap(np.angle(lin['_zs'] * corr)), '-', lw=1.4)
    if kerr is not None:
        ax[0].plot(fg, 20 * np.log10(np.abs(kerr['_zs'])), '--', lw=1.4,
                   label=f"Kerr a={kerr['a']:.2f} (res/noise {kerr['res_noise']:.2f})")
        ax[1].plot(fg, np.unwrap(np.angle(kerr['_zs'] * corr)), '--', lw=1.4)
    ax[0].axvline(r['fr'] / 1e9, color='gray', ls=':', lw=1)
    ax[1].axvline(r['fr'] / 1e9, color='gray', ls=':', lw=1)
    ax[0].set_ylabel('|S21| [dB]'); ax[0].legend(fontsize=8)
    ax[1].set_ylabel('fase − 2πfτ [rad]'); ax[1].set_xlabel('frecuencia [GHz]')
    corte = ''.join(f', {k} = {v:g}' for k, v in r['corte'].items())
    ax[0].set_title(f"run {r['run_id']}{corte}  |  P_in = {r['P_in_dBm']:.1f} dBm  |  "
                    f"<n> = {r['n_fotones']:.3g}  ({r['modelo']})")
    fig.tight_layout()
    plt.show()
    return fig