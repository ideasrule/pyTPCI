"""Initial condition for pyTPCI: a Parker-type wind on top of the cold hydrostatic base.

  * Mass-loss rate: energy-limited, Mdot = eff * 4 pi Rp^3 F_XUV cos(angle) / (G Mp), where F_XUV is the
    flux below 912 A in the Cloudy spectrum file and the 4 pi matches the 1D spherical geometry
    (every direction receives the slanted flux F cos(angle)).
  * Wind: steady transonic solution of the momentum equation for a prescribed theta(r) = P/rho = kT/(mu mH)
    in the PLUTO gravity (planet + stellar tide), integrated in both directions from the critical point;
    rho = Mdot / (4 pi r^2 v).
  * Base: theta goes from the base value (init_temp, init_mu) to the wind's with a tanh of width `width`
    centred on r_j, the one radius for which the pressure at r = 1 equals the boundary value.
  * theta of the wind: the first guess is isothermal (T_wind, mu_wind). A single temperature does not work
    for a deep potential: where the foot of the wind is dense and neutral, 10^4 K is far above thermal
    balance (Lya cooling >> heating), and a colder wind has too high a foot pressure to join the base.
    So, if a `cloudy_eq(r, rho, tag)` callback is given, T and mu are iterated to Cloudy's thermal
    equilibrium for the current density (damped in ln theta), with the wind rebuilt each time.
  * free_mdot=True (second stage): the energy-limited Mdot fixes where the cold layer ends, and that
    can be inconsistent with the heating (HD 209458b at eff = 0.1: Mdot 6x too low, the cold layer is
    put out to 1.16 Rp, where Cloudy heats it 30-250x faster than it cools, and it blows off at once).
    So theta(r) is then relaxed to max(Cloudy equilibrium, base value) everywhere, without the imposed
    tanh, and Mdot follows from the base pressure: the wind's v(r) does not depend on Mdot, so
    Mdot = 4 pi Rp^2 P_surf v(1) / theta(1). The damping is reduced when Mdot would jump by more than 3x.
Everything is in cgs here; build() writes ic.dat in PLUTO code units: r  rho  v  P  mu
"""
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
from scipy.ndimage import gaussian_filter1d

# PLUTO's constants (pluto.h), so that the base is hydrostatic in PLUTO's own gravity
G, kB, mH, Mearth, Msun, au, Rearth = 6.6726e-8, 1.3806505e-16, 1.6733e-24, 5.9736e27, 2.e33, 1.49597892e13, 6.378136e8

def f_xuv(spec_file):
    """Flux below 912 A [erg/cm2/s] of a Cloudy 'interpolate' table normalised by an 'intensity' line"""
    lognu, logf, total, lo, hi = [], [], None, 0.0, np.inf
    for line in open(spec_file):
        p = line.split()
        if p and p[0] == "continue":
            lognu.append(float(p[1])); logf.append(float(p[2]))
        elif p and p[0] == "intensity":
            total = float(p[1])
            assert "linear" in line and "range" in line and p[-1] == "Ryd"
            lo, hi = float(p[p.index("range") + 1]), float(p[p.index("to") + 1])
    nu, f = 10**np.array(lognu), 10**np.array(logf)
    ryd = 3.2898e15
    seg = lambda a, b: np.trapezoid(np.where((nu >= a) & (nu <= b), f, 0.0), nu)
    return total * seg(ryd, hi * ryd) / seg(lo * ryd, hi * ryd)

def wind(r, theta, g, mdot):
    """Transonic steady wind for theta(r): (v - theta/v) dv/dr = 2 theta/r - dtheta/dr + g =: N"""
    thp = np.gradient(theta, r)
    N = 2 * theta / r - thp + g
    up = np.where((N[:-1] < 0) & (N[1:] >= 0))[0]
    assert len(up), "no critical point in the domain"
    i = up[-1]                                    # outermost: N > 0 everywhere beyond it
    rc = r[i] - N[i] * (r[i + 1] - r[i]) / (N[i + 1] - N[i])
    th_c, thp_c = np.interp(rc, r, theta), np.interp(rc, r, thp)
    Np_c = (N[i + 1] - N[i]) / (r[i + 1] - r[i])
    vc = np.sqrt(th_c)
    slope = (thp_c + np.sqrt(thp_c**2 + 8 * vc**2 * Np_c)) / (4 * vc)
    rhs = lambda rr, y: np.interp(rr, r, N) / (np.exp(2 * y) - np.interp(rr, r, theta))   # y = ln v
    d = 3 * (r[i + 1] - r[i])
    lnv = np.empty_like(r)
    out, inn = r > rc + d, r < rc - d
    kw = dict(method="LSODA", rtol=1e-8, atol=1e-10)
    lnv[out] = solve_ivp(rhs, (rc + d, r[-1]), [np.log(vc + slope * d)], t_eval=r[out], **kw).y[0]
    lnv[inn] = solve_ivp(rhs, (rc - d, r[0]), [np.log(vc - slope * d)], t_eval=r[inn][::-1], **kw).y[0][::-1]
    mid = ~(out | inn)
    lnv[mid] = np.log(vc + slope * (r[mid] - rc))
    v = np.exp(lnv)
    return v, mdot / (4 * np.pi * r**2 * v), rc

def build(Rp, Mp, M_star, a, init_temp, init_mu, n_surf, spec_file, angle_deg, unit_density, unit_velocity,
          r_out=15.0, T_wind=1e4, mu_wind=0.64, eff=0.1, width=0.003, cloudy_eq=None, n_iter=12, tol=0.03,
          damp=0.5, free_mdot=True, n_iter_free=30, output_file="ic.dat"):
    # same expressions as params.h / init.c
    Rp_c = Rp / 6.378136e8 * Rearth; Mp_c = Mp / 5.97e27 * Mearth; Ms_c = M_star / 2e33 * Msun; a_c = a / 1.5e13 * au
    l_cm = a_c / (1 + Mp_c / Ms_c)
    x = np.concatenate([np.linspace(1.0, 1.5, 5001)[:-1], np.geomspace(1.5, r_out * 1.02, 3000)])
    r = x * Rp        # tpci.py and PLUTO lengths differ by 1e-4; r/Rp is what both use as the coordinate
    g = -G * Mp_c / r**2 + G * Ms_c / (a_c - r)**2 - G * (Ms_c + Mp_c) / a_c**3 * (l_cm - r)
    F = f_xuv(spec_file)
    mdot = eff * 4 * np.pi * Rp_c**3 * F * np.cos(np.radians(angle_deg)) / (G * Mp_c)
    th_c = kB * init_temp / (init_mu * mH)
    P_surf = n_surf * kB * init_temp
    print(f"Parker IC: F_XUV = {F:.4e} erg/cm2/s, energy-limited Mdot = {mdot:.4e} g/s (eff = {eff})", flush=True)

    def solve(th_hot, mu_hot):
        """Find r_j and return the full profile for a given wind theta(r) (>= th_c)"""
        def profile(xj):
            s = 0.5 * (1 + np.tanh((x - xj) / width))
            th_up = np.interp(np.maximum(x, xj + 3 * width), x, th_hot)   # flat below the transition
            theta = th_c + (th_up - th_c) * s
            v, rho, rc = wind(r, theta, g, mdot)
            return rho, v, rho * theta, init_mu + (mu_hot - init_mu) * s, rc
        resid = lambda xj: np.log(profile(xj)[2][0] / P_surf)
        lo = 1.0 + 2 * width
        assert resid(lo) < 0, "the wind's foot pressure exceeds the base pressure: no cold layer fits"
        hi = lo
        while resid(hi) < 0:
            lo, hi = hi, hi + 2 * width
            assert hi < 1.5, "no join radius below 1.5 Rp"
        xj = brentq(resid, lo, hi, xtol=1e-7)
        return (xj,) + profile(xj)

    th_hot = np.full(len(x), kB * T_wind / (mu_wind * mH))
    mu_hot = np.full(len(x), mu_wind)
    xj, rho, v, P, mu, rc = solve(th_hot, mu_hot)
    print(f"Parker IC: isothermal guess T = {T_wind:.0f} K, mu = {mu_wind}: r_join = {xj:.4f}, r_sonic = {rc / Rp:.3f}", flush=True)
    if cloudy_eq is not None:
        lx = np.linspace(0, np.log(x[-1]), 20000)                    # uniform ln r grid for the smoothing
        for it in range(n_iter):
            r_cl, T_cl, mu_cl = cloudy_eq(x, rho, f"ic{it}")
            # (static equilibrium has no expansion cooling, so T runs to 2-5e4 K far out. Capping it at
            # T_wind was tried: the foot then needs a higher density for the same Mdot, gets more neutral
            # and colder, and the iteration runs away until no cold layer fits.)
            th_cl = np.maximum(kB * np.interp(x, r_cl, T_cl) / (np.interp(x, r_cl, mu_cl) * mH), th_c)
            sm = gaussian_filter1d(np.interp(lx, np.log(x), np.log(th_cl)), 0.01 / (lx[1] - lx[0]), mode="nearest")
            th_new = np.exp(np.interp(np.log(x), lx, sm))
            above = x > xj + 3 * width
            err = np.abs(np.log(th_new / th_hot))[above]
            print(f"Parker IC: iteration {it}: max |dln theta| = {err.max():.3f} at r = {x[above][np.argmax(err)]:.3f}; "
                  f"r_join = {xj:.4f}, r_sonic = {rc / Rp:.3f}", flush=True)
            if err.max() < tol:
                break
            th_try = np.exp((1 - damp) * np.log(th_hot) + damp * np.log(th_new))
            mu_try = (1 - damp) * mu_hot + damp * np.interp(x, r_cl, mu_cl)
            try:
                xj, rho, v, P, mu, rc = solve(th_try, mu_try)
            except AssertionError as e:
                print("Parker IC: iteration stopped, keeping the last profile:", e, flush=True)
                break
            th_hot, mu_hot = th_try, mu_try
    if cloudy_eq is not None and free_mdot:
        theta = P / rho
        smooth = lambda th: np.exp(np.interp(np.log(x), lx, gaussian_filter1d(
            np.interp(lx, np.log(x), np.log(th)), 0.01 / (lx[1] - lx[0]), mode="nearest")))
        for it in range(n_iter_free):
            r_cl, T_cl, mu_cl = cloudy_eq(x, rho, f"icf{it}")
            mu = np.interp(x, r_cl, mu_cl)
            th_new = np.maximum(smooth(np.maximum(kB * np.interp(x, r_cl, T_cl) / (mu * mH), th_c)), th_c)
            err = np.abs(np.log(th_new / theta))
            d = damp
            while True:
                th_try = np.exp((1 - d) * np.log(theta) + d * np.log(th_new))
                v_try, _, rc = wind(r, th_try, g, 1.0)
                mdot_try = 4 * np.pi * r[0]**2 * P_surf * v_try[0] / th_try[0]
                if abs(np.log(mdot_try / mdot)) < np.log(3.0) or d < 0.02:
                    break
                d *= 0.5
            print(f"Parker IC: free-Mdot iteration {it}: max |dln theta| = {err.max():.3f} at r = {x[np.argmax(err)]:.3f}; "
                  f"damping {d:.3f}; Mdot {mdot:.3e} -> {mdot_try:.3e}; r_sonic = {rc / Rp:.3f}", flush=True)
            done = err.max() < 2 * tol and abs(np.log(mdot_try / mdot)) < 2 * tol
            theta, v, mdot = th_try, v_try, mdot_try
            rho = mdot / (4 * np.pi * r**2 * v); P = rho * theta
            if done:
                break
    np.savetxt(output_file, np.column_stack([x, rho / unit_density, v / unit_velocity,
                                              P / (unit_density * unit_velocity**2), mu]), fmt="%.10e")
    for xx in (1.0, 1.01, 1.03, 1.05, 1.1, 1.2, 1.5, 2, 3, 5, 10, 15):
        i = np.argmin(abs(x - xx))
        print(f"  r={xx:<5} rho={rho[i]:.3e} v={v[i]/1e5:.3e} km/s theta={P[i]/rho[i]/1e10:.1f} T={P[i]/rho[i]*mu[i]*mH/kB:.0f} mu={mu[i]:.2f}", flush=True)
    return dict(F_xuv=F, mdot=mdot, r_sonic=rc / Rp, r_join=xj)
