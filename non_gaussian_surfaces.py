import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib import cm

from scipy.stats import johnsonsu, johnsonsb, lognorm, norm, probplot

from j_johnson_M import f_johnson_M


# ----------------------------------------------------------------------
# Johnson-distribution helpers
# ----------------------------------------------------------------------

def _johnson_dist(coef, j_type):
    """
    Return a frozen scipy.stats distribution matching the Johnson curve
    described by `coef` = (gamma, delta, xi, lam) and `j_type`, or None for
    the degenerate 'ST' (two-point) boundary case, which isn't continuous.

    Reparametrizations used (Z ~ N(0,1) is the underlying standard normal):
      SU (unbounded): Y = xi + lam*sinh((Z-gamma)/delta)         -> johnsonsu
      SB (bounded):   Y = xi + lam/(1+exp(-(Z-gamma)/delta))     -> johnsonsb
      SL (lognormal): Y = xi + lam*exp((Z-gamma)/delta)          -> lognorm(s=1/delta, scale=lam*exp(-gamma/delta))
      SN (normal):    Y = (Z-gamma)/delta                        -> norm(loc=-gamma/delta, scale=1/delta)
    """
    gamma, delta, xi, lam = coef

    if j_type == 'SU':
        return johnsonsu(a=gamma, b=delta, loc=xi, scale=lam)
    elif j_type == 'SB':
        return johnsonsb(a=gamma, b=delta, loc=xi, scale=lam)
    elif j_type == 'SL':
        return lognorm(s=1.0/delta, loc=xi, scale=lam*np.exp(-gamma/delta))
    elif j_type == 'SN':
        return norm(loc=-gamma/delta, scale=1.0/delta)
    elif j_type == 'ST':
        return None
    else:
        raise ValueError(f'Unknown Johnson type: {j_type!r}')


def _sample_johnson(coef, j_type, size, rng=None):
    """
    Draw random samples from the Johnson distribution described by
    `coef` = (gamma, delta, xi, lam) and `j_type` (one of 'SN','SL','SB','SU','ST'),
    as returned by f_johnson_M().

    f_johnson_M can return ANY of the five Johnson families depending on the
    requested (skewness, kurtosis) pair -- which family you get is not a free
    choice, it is determined by where (skew, kurt) falls relative to the
    Johnson-system boundary curves. Sampling always with `johnsonsu` regardless
    of the fitted type (as an earlier version of this code did) silently
    applies the wrong inverse transform whenever the fit is SB, SL, SN or ST,
    and produces a surface whose actual moments can be wildly different from
    what was asked for.
    """
    if rng is None:
        rng = np.random.default_rng()

    if j_type == 'ST':
        # Degenerate two-point ("boundary") distribution: happens right on
        # the kurt = 1 + skew^2 boundary line, where the only distribution
        # with those moments is a two-valued (Bernoulli-like) variable.
        # f_johnson_M packs it as: delta = P(value == lam), xi/lam = the two values.
        _, delta, xi, lam = coef
        p = delta
        u = rng.random(size=size)
        return np.where(u < p, lam, xi)

    return _johnson_dist(coef, j_type).rvs(size=size, random_state=rng)


# ----------------------------------------------------------------------
# Small numerical helpers
# ----------------------------------------------------------------------

def _hhcf_1d(z, N):
    """1-D height-height correlation function G(r) = sqrt(<(z(x+r)-z(x))^2>)."""
    hhcf1d = np.zeros(N // 2)
    for ndif in range(N // 2):
        surf1 = z[:, :N - ndif]
        surf2 = z[:, ndif:]
        hhcf1d[ndif] = np.sqrt(np.mean((surf1 - surf2) ** 2))
    return hhcf1d


def _moments(v):
    """mean, std, skewness, excess-free (Pearson) kurtosis of a flattened array."""
    from scipy.stats import skew, kurtosis
    return np.mean(v), np.std(v), skew(v), kurtosis(v, fisher=False)


# ----------------------------------------------------------------------
# Visualization
# ----------------------------------------------------------------------

_ACCENT = '#e8590c'
_ACCENT2 = '#1971c2'
_INK = '#1a1a1a'
_PANEL_BG = '#fafafa'


def plot_surface_gallery(z_gs, z_ngs, coef, j_type, err, requested, N, corr=True,
                          cmap='inferno', threshold=None):
    """
    Build a single "gallery + report" figure summarizing the generated
    surface: 3D views, a top-down heightmap, a binarized mask, a histogram
    with the fitted Johnson PDF overlaid, a Q-Q plot against that fit, the
    height-height correlation function, and a text report card comparing
    requested vs. achieved statistics.

    Returns the matplotlib Figure (does not call plt.show()).
    """
    mu_req, rms_req, skew_req, kurt_req, corlength_x, corlength_y, alpha = requested
    v_ngs = z_ngs.flatten()
    mean_a, std_a, skew_a, kurt_a = _moments(v_ngs)
    if threshold is None:
        threshold = np.mean(z_ngs)

    fig = plt.figure(figsize=(16, 11), facecolor='white')
    gs = GridSpec(3, 3, figure=fig, height_ratios=[1.15, 1, 0.9],
                  hspace=0.5, wspace=0.35, top=0.90, bottom=0.05, left=0.05, right=0.97)
    Xg, Yg = np.meshgrid(np.arange(N), np.arange(N))

    # --- Row 0: 3D Gaussian | 3D non-Gaussian | top-down heatmap --- #
    ax0 = fig.add_subplot(gs[0, 0], projection='3d')
    ax0.plot_surface(Xg, Yg, z_gs, cmap=cmap, edgecolor='none', antialiased=True)
    ax0.set_title('Gaussian surface\n(spatial structure only)', fontsize=10.5, color=_INK, pad=14)
    ax0.view_init(35, -60)
    ax0.set_xticklabels([]); ax0.set_yticklabels([])

    ax1 = fig.add_subplot(gs[0, 1], projection='3d')
    ax1.plot_surface(Xg, Yg, z_ngs, cmap=cmap, edgecolor='none', antialiased=True)
    ax1.set_title(f'Non-Gaussian surface\n(type {j_type})', fontsize=10.5, color=_INK, pad=14)
    ax1.view_init(35, -60)
    ax1.set_xticklabels([]); ax1.set_yticklabels([])

    ax2 = fig.add_subplot(gs[0, 2])
    im = ax2.imshow(z_ngs, cmap=cmap, origin='lower')
    ax2.set_title('Top-down heightmap', fontsize=10.5, color=_INK)
    ax2.set_xticks([]); ax2.set_yticks([])
    cb = fig.colorbar(im, ax=ax2, fraction=0.046, pad=0.04)
    cb.set_label('height', fontsize=9)

    # --- Row 1: binarized mask | histogram+pdf | Q-Q plot --- #
    ax3 = fig.add_subplot(gs[1, 0])
    ax3.imshow(z_ngs > threshold, cmap='gray', origin='lower')
    ax3.set_title(f'Binarized (z > {threshold:.2f})', fontsize=10.5, color=_INK)
    ax3.set_xticks([]); ax3.set_yticks([])

    ax4 = fig.add_subplot(gs[1, 1])
    ax4.set_facecolor(_PANEL_BG)
    counts, bins, _ = ax4.hist(v_ngs, bins=120, density=True, color=_ACCENT2,
                                alpha=0.75, edgecolor='none', label='sampled surface')
    xs = np.linspace(bins[0], bins[-1], 400)
    dist = _johnson_dist(coef, j_type)
    if dist is not None:
        ax4.plot(xs, dist.pdf(xs), color=_ACCENT, lw=2.2, label=f'fitted Johnson-{j_type} PDF')
    else:
        # ST: discrete two-point distribution, draw stems instead of a PDF
        _, p, x_lo, x_hi = coef
        ax4.vlines([x_lo, x_hi], 0, [1 - p, p], color=_ACCENT, lw=3, label='fitted two-point mass')
    ax4.axvline(mean_a, color=_INK, ls='--', lw=1, alpha=0.6)
    ax4.set_title('Height distribution', fontsize=10.5, color=_INK)
    ax4.set_xlabel('height'); ax4.set_ylabel('density')
    ax4.legend(fontsize=8, frameon=False)

    ax5 = fig.add_subplot(gs[1, 2])
    ax5.set_facecolor(_PANEL_BG)
    if dist is not None:
        (osm, osr), (slope, intercept, r) = probplot(v_ngs, dist=dist, fit=True)
        ax5.scatter(osm, osr, s=4, color=_ACCENT2, alpha=0.5)
        ax5.plot(osm, slope * osm + intercept, color=_ACCENT, lw=2)
        ax5.set_title(f'Q-Q plot vs. fit  (R={r:.4f})', fontsize=10.5, color=_INK)
    else:
        ax5.text(0.5, 0.5, 'Q-Q plot not applicable\n(ST is a discrete\ntwo-point distribution)',
                  ha='center', va='center', fontsize=9.5, color=_INK, transform=ax5.transAxes)
        ax5.set_xticks([]); ax5.set_yticks([])
        ax5.set_title('Q-Q plot', fontsize=10.5, color=_INK)
    ax5.set_xlabel('theoretical quantiles'); ax5.set_ylabel('sample quantiles')

    # --- Row 2: HHCF | statistics report card --- #
    ax6 = fig.add_subplot(gs[2, 0])
    ax6.set_facecolor(_PANEL_BG)
    if corr:
        hhcf1d = _hhcf_1d(z_ngs, N)
        r_axis = np.arange(1, N // 2)
        ax6.loglog(r_axis, hhcf1d[1:], color=_ACCENT2, lw=1.8)
        ax6.set_xlabel('r'); ax6.set_ylabel('G(r)')
        ax6.set_title('Height-height correlation', fontsize=10.5, color=_INK)
        ax6.grid(True, which='both', alpha=0.25)
    else:
        ax6.text(0.5, 0.5, 'Set corr=True to compute\nthe correlation function',
                  ha='center', va='center', fontsize=9.5, color=_INK, transform=ax6.transAxes)
        ax6.set_xticks([]); ax6.set_yticks([])
        ax6.set_title('Height-height correlation', fontsize=10.5, color=_INK)

    ax7 = fig.add_subplot(gs[2, 1:])
    ax7.axis('off')
    rows = [
        ('mean',      mu_req,   mean_a),
        ('rms (std)', rms_req,  std_a),
        ('skewness',  skew_req, skew_a),
        ('kurtosis',  kurt_req, kurt_a),
    ]

    def _pct_err(req, ach):
        if abs(req) < 1e-9:
            return None
        return 100.0 * (ach - req) / abs(req)

    header = f"{'metric':<11}{'requested':>12}{'achieved':>12}{'error':>10}"
    lines = [header, '-' * len(header)]
    for name, req, ach in rows:
        e = _pct_err(req, ach)
        e_str = f"{e:+.1f}%" if e is not None else "  n/a"
        lines.append(f"{name:<11}{req:>12.3f}{ach:>12.3f}{e_str:>10}")
    lines.append('')
    lines.append('* mean is n/a: surface is shifted so min height = 0 (rms/skew/kurt')
    lines.append('  are shift-invariant, so those three stay directly comparable)')
    lines.append(f"Johnson type fitted : {j_type}" + (f"   ({err})" if err else ""))
    lines.append(f"Correlation lengths : x={corlength_x}, y={corlength_y}   roughness alpha={alpha}")
    lines.append(f"Grid size           : {N} x {N}")

    ax7.text(0.0, 1.0, 'STATISTICAL REPORT', fontsize=13.5, fontweight='bold',
             color=_INK, transform=ax7.transAxes, va='top')
    ax7.text(0.0, 0.86, '\n'.join(lines), fontsize=10.5, family='monospace',
             color=_INK, transform=ax7.transAxes, va='top')

    fig.suptitle('Non-Gaussian Self-Affine Surface -- Summary', fontsize=15, fontweight='bold', y=0.98)
    return fig


def print_report(coef, j_type, err, requested, z_ngs):
    """Console version of the statistics report card (no plotting)."""
    mu_req, rms_req, skew_req, kurt_req, corlength_x, corlength_y, alpha = requested
    mean_a, std_a, skew_a, kurt_a = _moments(z_ngs.flatten())

    def _pct_err(req, ach):
        if abs(req) < 1e-9:
            return None
        return 100.0 * (ach - req) / abs(req)

    rows = [('mean', mu_req, mean_a), ('rms (std)', rms_req, std_a),
            ('skewness', skew_req, skew_a), ('kurtosis', kurt_req, kurt_a)]

    width = 48
    print('+' + '-' * width + '+')
    print('| {:^{w}} |'.format('NON-GAUSSIAN SURFACE -- STATISTICAL REPORT', w=width - 2))
    print('+' + '-' * width + '+')
    print(f"| {'metric':<11}{'requested':>12}{'achieved':>12}{'error':>10} |")
    print('|' + '-' * width + '|')
    for name, req, ach in rows:
        e = _pct_err(req, ach)
        e_str = f"{e:+.1f}%" if e is not None else "n/a"
        print(f"| {name:<11}{req:>12.3f}{ach:>12.3f}{e_str:>10} |")
    print('+' + '-' * width + '+')
    print(f"Johnson type fitted : {j_type}" + (f"  ({err})" if err else ""))
    print(f"Correlation lengths : x={corlength_x}, y={corlength_y}   roughness alpha={alpha}")
    print("* mean is n/a: surface is shifted so min height = 0 (rms/skew/kurt are shift-invariant)")


# ----------------------------------------------------------------------
# Main entry point
# ----------------------------------------------------------------------

def SAimage_fft_2(N=500, rms=3, skewness=1, kurtosis=3, corlength_x=10, corlength_y=10, alpha=0.9,
                  threshold=None, non_Gauss=True, corr=True, invert=False, seed=None,
                  show=True, clip_percentiles=None, plot_style='rich'):
    """
    Generate a (non-)Gaussian self-affine surface, following Yang et al.,
    CMES vol.103 no.4 pp.251-279 (2014): a Gaussian surface with the requested
    correlation lengths/roughness exponent is generated first, then an
    independent 1-D noise sample matching (mean, rms, skewness, kurtosis) is
    generated via the Johnson system and remapped onto the Gaussian surface's
    rank order, so the final surface keeps the Gaussian surface's spatial
    correlation but the requested marginal (height) distribution.

    Parameters
    ----------
    invert : bool, default False
        If True, flips the sign of the final non-Gaussian surface before
        shifting it to be non-negative. NOTE: this also flips the sign of its
        skewness relative to what you requested. Off by default so the
        surface's measured skewness matches `skewness`.
    seed : int or None
        Seed for reproducibility.
    show : bool
        If False, build the figure(s) but don't call plt.show() (useful for
        headless/batch runs and tests). The figure(s) are still returned.
    clip_percentiles : (low, high) or None, default None
        If given (e.g. (0.1, 99.9)), clips the sampled noise to these
        percentiles before renormalizing. Off by default: for high-kurtosis
        / SU-type requests the tails ARE the kurtosis, so clipping them
        systematically pulls the surface's measured kurtosis below what you
        asked for. Only enable this if you need to bound the height range
        and can accept some undershoot on kurtosis.
    plot_style : {'rich', 'classic', 'none'}, default 'rich'
        'rich'    -- one combined figure: 3D views, heightmap, histogram with
                     fitted Johnson PDF, Q-Q plot, correlation function, and a
                     text statistics report card (requested vs. achieved).
        'classic' -- the original three separate figures (3D surface,
                     binarized mask, optional loglog correlation function).
        'none'    -- skip plotting entirely, just return the height array.

    Returns
    -------
    z_ngs : ndarray, shape (N, N)
        The generated surface (non-negative). If `non_Gauss=False`, the
        Gaussian surface is returned instead.
    """
    rng = np.random.default_rng(seed)

    # --- Gaussian Surface Generation --- #
    txmin = tymin = -N/2
    txmax = tymax = N/2
    dtx = dty = (txmax - txmin) / N
    tx = np.arange(txmin, txmax, dtx)
    ty = np.arange(tymin, tymax, dty)

    R = np.zeros((N+1, N+1))
    TX, TY = np.meshgrid(tx, ty, indexing='ij')
    r = np.sqrt((TX/corlength_x)**2 + (TY/corlength_y)**2)
    R[:N, :N] = rms**2 * np.exp(-np.abs(r)**(2*alpha))

    FR = np.fft.fft2(R, s=[N, N])
    AMPR = np.sqrt(dtx**2 + dty**2) * np.abs(FR)

    X = rng.random((N, N))
    X = (X - np.mean(X)) / np.std(X)
    XF = np.fft.fft2(X, s=[N, N])

    YF = XF * np.sqrt(AMPR)
    z = np.real(np.fft.ifft2(YF, s=[N, N]))
    z = (z - np.mean(z)) * rms / np.std(z)
    z_gs = z.copy()

    if threshold is None:
        threshold = np.mean(z_gs)

    if not non_Gauss:
        if plot_style != 'none':
            fig = plt.figure()
            ax = fig.add_subplot(111, projection='3d')
            Xg, Yg = np.meshgrid(np.arange(N), np.arange(N))
            ax.plot_surface(Xg, Yg, z_gs, cmap='viridis', edgecolor='none')
            ax.view_init(30, 150)
            plt.title('Gaussian Surface')
            if show:
                plt.show()

            plt.figure()
            plt.imshow(z_gs < threshold, cmap='gray')
            plt.title('Binarized Gaussian Surface')
            if show:
                plt.show()
        return z_gs - np.min(z_gs)  # Return surface shifted to be positive

    # --- Non-Gaussian Noise --- #
    coef, j_type, err = f_johnson_M(0, rms, skewness, kurtosis)
    gamma, delta, xi, lam = coef
    print(f"Johnson type={j_type} params: gamma={gamma}, delta={delta}, xi={xi}, lambda={lam}"
          + (f"  [note: {err}]" if err else ""))

    # IMPORTANT: dispatch on the fitted type -- see _sample_johnson docstring.
    # Sampling always with johnsonsu regardless of the fitted type (as an
    # earlier version of this code did) silently applies the wrong inverse
    # transform whenever the fit isn't 'SU', giving order-of-magnitude wrong
    # moments.
    z_ngn = _sample_johnson(coef, j_type, size=(N, N), rng=rng)

    if clip_percentiles is not None:
        z_ngn = np.clip(z_ngn, *np.percentile(z_ngn, clip_percentiles))
    z_ngn = (z_ngn - np.mean(z_ngn)) / np.std(z_ngn) * rms

    # --- Rank-based Mapping --- #
    v_gs = z_gs.flatten()
    v_ngn = z_ngn.flatten()

    Igs = np.argsort(v_gs)
    vs_ngn_sorted = np.sort(v_ngn)

    v_ngs = np.zeros_like(v_gs)
    v_ngs[Igs] = vs_ngn_sorted
    z_ngs = v_ngs.reshape((N, N))

    if invert:
        z_ngs = -z_ngs
    z_ngs = z_ngs - np.min(z_ngs)

    if threshold is not None and threshold < np.min(z_ngs):
        # threshold was computed on the pre-shift Gaussian surface; re-center
        # it on the actual non-Gaussian height range for the binarized view
        threshold = np.mean(z_ngs)

    # --- Plotting --- #
    if plot_style == 'rich':
        requested = (0, rms, skewness, kurtosis, corlength_x, corlength_y, alpha)
        fig = plot_surface_gallery(z_gs, z_ngs, coef, j_type, err, requested, N,
                                    corr=corr, threshold=threshold)
        print_report(coef, j_type, err, requested, z_ngs)
        if show:
            plt.show()

    elif plot_style == 'classic':
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        Xg, Yg = np.meshgrid(np.arange(N), np.arange(N))
        ax.plot_surface(Xg, Yg, z_ngs, cmap='hot', edgecolor='black')
        ax.set_xlim(0, N)
        ax.set_ylim(0, N)
        ax.set_zlim(np.min(z_ngs), np.max(z_ngs))
        ax.view_init(30, 150)
        plt.title('Non-Gaussian Surface')
        if show:
            plt.show()

        plt.figure()
        plt.imshow(z_ngs > threshold, cmap='gray')
        plt.title('Binarized Non-Gaussian Surface')
        if show:
            plt.show()

        if corr:
            hhcf1d = _hhcf_1d(z_ngs, N)
            plt.figure()
            plt.loglog(np.arange(N//2), hhcf1d)
            plt.grid()
            plt.xlabel('log(r(nm))')
            plt.ylabel('log(G(r) (nm))')
            plt.title('1-D height-height correlation function')
            if show:
                plt.show()

    # plot_style == 'none' -> skip plotting entirely

    return z_ngs