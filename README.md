# Non-Gaussian Surface Simulations

**Author:** Sebastian Korsak

# 🌄 Gaussian & Non-Gaussian Self-Affine Surface Generator

This Python module generates **2D random surfaces** — rough, textured height maps like the ones you'd measure with an AFM, a profilometer, or see on a fracture or corrosion surface — with a prescribed **spatial correlation** (how "wiggly" the surface is and over what length scale) and a prescribed **height distribution** (how the heights themselves are spread out: symmetric or lopsided, gently rounded or spiky).

It implements a spectral (FFT-based) synthesis method for the spatial part, and the **Johnson system of distributions** for the statistical part, following:

> Yang et al., *CMES*, vol. 103, no. 4, pp. 251–279, 2014

---

## 💭 What are we actually simulating, and why?

Real rough surfaces are almost never perfectly Gaussian. Think of a corroded metal plate: most of it sits close to some average height, but a few deep pits go much further down than a Gaussian bell curve would predict — that's a *skewed*, *heavy-tailed* surface, not a symmetric bumpy one. A machined surface, by contrast, might have its peaks gently truncated and its valleys smoothed — a different, *bounded* shape.

To simulate that faithfully, this module asks you for exactly the four numbers a surface scientist usually measures or specifies:

- **mean** — the average height,
- **rms (root-mean-square roughness)** — how far heights typically wander from the mean,
- **skewness** — whether the surface is lopsided (tall spiky peaks vs. deep narrow valleys, or the reverse),
- **kurtosis** — how "spiky vs. flat" the height distribution is (occasional extreme excursions vs. everything hovering near the same range),

and combines them with a **correlation length** and a **roughness (Hurst) exponent** that describe how far apart two points on the surface need to be before their heights stop "knowing about" each other, and how rough the surface looks at each scale in between. The goal is a surface that looks and behaves statistically like the real thing it's meant to stand in for — not just "random bumps," but bumps with the right personality.

---

## 🧠 The Math, and the Intuition Behind It

### 1. Spatial correlation: the Gaussian scaffold

We first build a **Gaussian** random surface with the right spatial structure — this is the "scaffold" that the real, non-Gaussian heights will later be poured into.

The spatial autocorrelation function is modeled as a stretched exponential:

```
R(tx, ty) = rms² · exp( − ( sqrt((tx/ξx)² + (ty/ξy)²) )^(2α) )
```

- `ξx`, `ξy` — correlation lengths along x and y. Beyond this distance, two points on the surface are essentially statistically independent.
- `α` — the roughness (Hurst) exponent. Small `α` (→0) gives a jagged, fractal-like surface; `α = 0.5` corresponds to Brownian-motion-like roughness; `α → 1` gives smoother, more rounded undulations.

By the **Wiener–Khinchin theorem**, the power spectrum of the surface is just the Fourier transform of `R`. So the recipe is:

1. Build `R(tx, ty)` on a grid, take its 2D FFT to get the power spectrum `|FR|`.
2. Generate white noise, take *its* FFT.
3. Multiply the noise spectrum by `sqrt(|FR|)` — this "colors" the noise with the correlation structure we want.
4. Inverse FFT back to real space, and rescale to the requested `rms`.

The result, `z_gs`, is a **Gaussian** surface (skewness 0, kurtosis 3) with exactly the correlation lengths and roughness exponent you asked for. If that's all you need, this is the whole story — see `non_Gauss=False`.

### 2. Shaping the marginal distribution: the Johnson system

Gaussian heights are the easy case. But if you want skewness ≠ 0 or kurtosis ≠ 3, there is no single closed-form distribution that lets you dial in mean, rms, skewness, and kurtosis independently — you need a *family* of distributions flexible enough to cover the whole feasible (skewness, kurtosis) plane.

That's what the **Johnson system** does. It says: take a standard normal variable `Z`, and pass it through one of four transformations to get a variable `Y` with the moments you want:

| Type | Name | Transform `Y = ...` | Shape |
|------|------|----------------------|-------|
| `SN` | Normal | `(Z − γ)/δ` | plain Gaussian (skew=0, kurt=3) |
| `SL` | Lognormal | `ξ + λ·exp((Z − γ)/δ)` | one-sided, unbounded on one side |
| `SB` | Bounded | `ξ + λ / (1 + exp(−(Z − γ)/δ))` | bounded above **and** below |
| `SU` | Unbounded | `ξ + λ·sinh((Z − γ)/δ)` | unbounded both sides, heavy tails |

plus a degenerate boundary case `ST`, a two-point (Bernoulli-like) distribution that shows up exactly on the line `kurtosis = 1 + skewness²`, where no continuous distribution can match the requested moments.

**Given your (mean, rms, skewness, kurtosis), which of these five you get is not your choice** — it's determined entirely by where that (skewness, kurtosis) point falls relative to the Johnson system's boundary curves. A request for `(skewness=1, kurtosis=3)`, for instance, lands you in **SB** territory, not SU. The module fits the right `(γ, δ, ξ, λ)` parameters for whichever family applies, and samples from *that* family — this dispatch is the key correctness fix described below.

### 3. Marrying the two: rank-order mapping

We now have two independent things: a Gaussian field `z_gs` with the right *spatial correlation*, and an independent noise sample `z_ngn` with the right *marginal distribution* (mean/rms/skew/kurtosis) but no spatial structure at all.

The trick — a **rank-order (copula-style) mapping** — combines them for free: sort both arrays, and reassign values so that whichever pixel had the *k*-th smallest value in `z_gs` now gets the *k*-th smallest value of `z_ngn`. This preserves the *spatial pattern* of the Gaussian field (where the hills and valleys are) while giving every pixel a height drawn from the exact target distribution. The final surface is then shifted so its minimum height is zero, matching how physical height maps (e.g. AFM scans) are usually reported.

---

## 🔧 Function

```python
SAimage_fft_2(N=500, rms=3, skewness=1, kurtosis=3, corlength_x=10, corlength_y=10,
              alpha=0.9, threshold=None, non_Gauss=True, corr=True, invert=False,
              seed=None, show=True, clip_percentiles=None, plot_style='rich')
```

---

## 📥 Parameters

| Parameter          | Type            | Description |
|--------------------|-----------------|-------------|
| `N`                | `int`           | Size of the generated surface (`N x N`) |
| `rms`              | `float`         | Target **root-mean-square roughness** |
| `skewness`         | `float`         | Target skewness of the height distribution (ignored if `non_Gauss=False`) |
| `kurtosis`         | `float`         | Target kurtosis of the height distribution (ignored if `non_Gauss=False`) |
| `corlength_x`      | `float`         | Correlation length along the **x-axis** |
| `corlength_y`      | `float`         | Correlation length along the **y-axis** |
| `alpha`            | `float`         | **Roughness exponent** (Hurst-like exponent) |
| `non_Gauss`        | `bool`          | `True` (default) generates the full non-Gaussian surface; `False` returns the plain Gaussian scaffold |
| `threshold`        | `float or None` | Height used for the binarized view. Defaults to the surface's own mean if not given |
| `corr`             | `bool`          | Whether to compute/plot the height-height correlation function |
| `invert`           | `bool`          | Flip the surface's sign before shifting non-negative (also flips skewness — off by default) |
| `seed`             | `int or None`   | Random seed, for reproducible surfaces |
| `show`             | `bool`          | Call `plt.show()` or just build and return the figure(s) (useful headless) |
| `clip_percentiles` | `(low, high) or None` | Optionally clip extreme noise values before rescaling (off by default — see note below) |
| `plot_style`       | `'rich'`, `'classic'`, `'none'` | Which visualization to produce (see **New Features**) |

---

## 📈 Output

Returns a 2D NumPy array:

```python
surface = SAimage_fft_2(...)
```

- Height values distributed exactly as requested (mean, rms, skewness, kurtosis), or Gaussian if `non_Gauss=False`
- Non-negative, shifted so the minimum height is 0
- Spatial correlation matching `corlength_x`, `corlength_y`, and `alpha`

---

## ✨ New Features

This version fixes a real correctness bug and adds a full visual/statistical reporting layer on top of the original port.

### 1. Correct Johnson-type dispatch (the bug fix)

The previous implementation always sampled the non-Gaussian noise with `scipy.stats.johnsonsu`, regardless of which Johnson family the moment-fit actually returned. Since the module's own defaults (`skewness=1, kurtosis=3`) fit to an **SB** (bounded) curve, not SU, this silently fed SB parameters through the SU (sinh) transform — producing a sample with skewness ≈ −29 and kurtosis ≈ 2095 instead of the requested 1 and 3. The generator now dispatches on the fitted type (`SN`, `SL`, `SB`, `SU`, or the degenerate `ST` boundary case) and samples from the correct family every time, verified to land within ~1% of the requested moments across all five types.

### 2. Rich visualization gallery

With `plot_style='rich'` (the default), a single combined figure is produced with:

- Side-by-side **3D Gaussian vs. non-Gaussian surfaces**, same colormap for a fair visual comparison
- A **top-down heightmap** with colorbar
- The **binarized mask** at the chosen threshold
- A **histogram of heights** with the fitted Johnson PDF overlaid (or, for the `ST` boundary case, the two-point probability mass shown as stems)
- A **Q-Q plot** against that fit, annotated with the correlation coefficient
- The **height-height correlation function** (log-log)

### 3. Statistical report card

Alongside the gallery (and also printed to the console via `print_report`), a report compares **requested vs. achieved** mean, rms, skewness, and kurtosis, with percent error, plus the Johnson type that was fitted and your correlation-length/roughness settings — so you can see at a glance how faithfully the surface matches what you asked for.

### 4. Configurable, less lossy noise handling

- `clip_percentiles` is now **off by default** (it used to always clip to the 0.1–99.9th percentile). For heavy-tailed (SU, high-kurtosis) requests, those extreme values *are* the kurtosis — clipping them silently biased the achieved kurtosis downward. Pass `clip_percentiles=(0.1, 99.9)` yourself if you specifically need bounded heights and can accept some undershoot.
- `invert` is now **off by default** (it used to always flip the surface's sign before shifting it positive, which silently flipped the sign of the achieved skewness relative to what was requested). Set `invert=True` if you want that mirrored convention back.

### 5. Quality-of-life additions

- `seed=` for reproducible surfaces
- `show=False` to build figures without blocking / for headless & batch runs
- `plot_style='classic'` keeps the original three separate figures if you prefer them; `plot_style='none'` skips plotting entirely
- The autocorrelation matrix `R` is now built with vectorized NumPy instead of a pure-Python double loop — same output, noticeably faster for large `N`

---

## 🌍 Usage

```python
import numpy as np
from non_gaussian_surfaces import SAimage_fft_2

# Full non-Gaussian surface, rich visualization + report:
z = SAimage_fft_2(N=500, rms=3, skewness=1, kurtosis=3,
                   corlength_x=10, corlength_y=10, alpha=0.9,
                   seed=42, plot_style='rich')

# Just the array, no plots at all (e.g. inside a batch job):
z = SAimage_fft_2(N=500, rms=3, skewness=1.5, kurtosis=6,
                   corlength_x=15, corlength_y=8, alpha=0.7,
                   seed=1, show=False, plot_style='none')

# Plain Gaussian surface only:
z_gauss = SAimage_fft_2(N=256, rms=1.0, corlength_x=20, corlength_y=20,
                          alpha=0.8, non_Gauss=False)
```

---

## ⚠️ Notes

- Skewness and kurtosis aren't a free choice of *distribution family* — the Johnson system picks the family for you based on where your requested (skewness, kurtosis) falls. The console output and report card always tell you which type (`SN`/`SL`/`SB`/`SU`/`ST`) was used.
- The report's **mean** row shows `n/a` for error — that's expected: the surface is deliberately shifted so its minimum height is 0, so mean isn't directly comparable to the request, while rms/skewness/kurtosis are shift-invariant and remain meaningful.
- The rank-order mapping assumes a Gaussian-copula-like relationship between the spatial field and the target marginal — it preserves *where* the extremes are, not a rigorous joint model beyond that.
- Spatial correlation is implemented via FFT filtering based on the **Wiener–Khinchin theorem**.

---

## 📜 References

- Yang, F., et al. *Statistical generation of 3D rough surfaces with arbitrary correlation*, CMES, vol. 103, no. 4, 2014.
- Persson, B.N.J. *Theory of rubber friction and contact mechanics*.
- Dave (2021). *Johnson Curve Toolbox*, MATLAB Central File Exchange.

## 👥 Credits

- Python translation by Max Pierini @ [EpiData.it](https://epidata.it)
- Ported from original MATLAB toolbox by Dave (2021)
- Matlab code from dr. Vasilis Costantoudis
- Extended and maintained by Sebastian Korsak
- Johnson-type dispatch bug fix, rich visualization gallery, and statistical report card added with Claude

---

# Results

![image](https://github.com/user-attachments/assets/be9a6f7b-90c5-4240-89bc-78a295d5c9bb)

![image](https://github.com/user-attachments/assets/a328b67c-366d-4d67-ab96-d95e5e2e6f12)