# Olympus 593 — LPC Mean-Line Reverse Engineering

MATLAB study of the seven-stage low-pressure compressor (LPC) of the Rolls-Royce Olympus 593, using published engine-cycle data to reconstruct a stage-by-stage mean-line design. The method follows Cohen, Rogers & Saravanamuttoo, *Gas Turbine Theory*, with radial analysis using the classical vortex laws from Dixon.

> **Status:** the LPC mean-line workflow is implemented and self-checked. The companion scripts preserved under `development/` and `extended_study/` are earlier or extended variants, not additional validated engines. There is **no turbine model** in this repository and **no experimental validation** of any result.

---

## Engineering objective

Starting from inlet conditions and overall pressure ratio, reconstruct:

- inlet static state and annulus geometry from continuity;
- per-stage temperature rise, velocity triangles and reaction;
- stage pressure ratio using polytropic efficiency;
- De Haller number and rotor deflection as diffusion-risk indicators;
- annulus variation through the compressor;
- hub-to-tip variation of the flow using selectable vortex laws;
- self-consistency between the radial solution and the mean-line solution.

The script does not model blade profiles, losses, cooling, secondary flow or off-design operation. It is a mean-line reconstruction, not a CFD or rig test.

---

## Input parameters

Values below are read from the `INPUT` block of `LPC_RUN1.m` (lines 20–50). The corresponding line numbers are given so a reviewer can trace each value to source.

| Parameter | Symbol | Value | Unit | Line |
|---|---|---:|---|---:|
| Ratio of specific heats | `gam` | 1.4 | — | 20 |
| Gas constant | `R` | 287 | J/(kg·K) | 21 |
| Mass flow | `mdot` | 131.2447 | kg/s | 24 |
| LP spool speed | `N_rpm` | 5819 | rpm | 25 |
| Inlet total temperature | `T01` | 389.86 | K | 26 |
| Inlet total pressure | `P01` | 73.66e3 | Pa | 27 |
| Polytropic efficiency | `eta_p` | 0.8782 | — | 28 |
| Inlet hub/tip ratio | `nu_ht` | 0.27 | — | 29 |
| Inlet Mach number (no swirl) | `M1` | 0.549 | — | 30 |
| Overall pressure ratio | `PR_total` | 4.1 | — | 31 |
| Stage count | `N_stages` | 7 | — | 32 |
| Stage 1 temperature rise | `dT0_stage1` | 20.0 | K | 36 |
| Stage 2 temperature rise | `dT0_stage2` | 25.0 | K | 37 |
| Work-done factors (stages 1–7) | `lambda` | 0.98, 0.93, 0.88, 0.83, 0.83, 0.83, 0.83 | — | 39 |
| Reaction (stage 1 is an output) | `Lambda` | NaN, 0.70, 0.50, 0.50, 0.50, 0.50, 0.50 | — | 41 |
| Radial vortex law | `VORTEX_LAW` | `'exponential'` | — | 43 |
| De Haller limit | `DEHALLER_LIMIT` | 0.72 | — | 50 |

Thermodynamic assumptions: perfect gas with constant `gam` and `R`; polytropic (not isentropic) stage pressure-ratio calculation; no inlet guide vanes on stage 1 (absolute inlet swirl is zero, reaction is solved); constant mean radius through the annulus; annulus areas evaluated at **stator exit** (station 3), not rotor exit; ideal radial equilibrium of the selected vortex law. The stage-1 tip is already transonic in the input condition, which is why stages 1 and 2 are prescribed manually rather than equalised.

---

## Mean-line methodology (`LPC_RUN1.m`)

`LPC_RUN1.m` is the primary script (644 lines). The calculation sequence is:

1. **Inlet state (lines 62–94).** Static temperature, density, speed of sound and axial velocity from `M1`. Annulus radii from continuity at the specified hub/tip ratio. Blade speed `U = ω·r_mean`. Flow coefficient `φ = Ca/U` is an **output**, not an input.
2. **Stage loading balance (lines 97–123).** Stages 3–6 share a common temperature rise found by bisection; stage 7 is then recovered so the prescribed overall pressure ratio of 4.1 is met exactly.
3. **Mean-line velocity triangles (lines 126–210).**
   - Stage 1: `α₁ = 0`, `Cw₁ = 0`, `Cw₂ = Cp·ΔT₀/(λ·U)`; reaction is a computed output.
   - Stages 2–7: reaction is prescribed and the `tan β₁ ± tan β₂` relation is solved for the triangle.
   - Stage pressure ratio: `PR_st = (1 + ΔT₀/T₀₁)^(γ·η_p/(γ−1))`.
   - De Haller number: `W₂/W₁ = hypot(Ca, U−Cw₂) / hypot(Ca, U−Cw₁)`.
   - Rotor deflection: `β₁ − β₂`, in degrees.
4. **Annulus through-flow (lines 213–253).** Stator-exit swirl is `Cw₃(i) = Cw₁(i+1)` for stages 1–6 and zero at exit. With constant mean radius, `A₃ = ṁ/(ρ₃·Ca)` and blade height `h = A₃/(2π·r_mean)`.
5. **Radial variation (lines 256–372).** Root, mean and tip sections with:
   - `free`: `Cw·r = const`, `Ca` uniform;
   - `first_power`: `Cw₁ = a·r − b/r`, `Cw₂ = a·r + b/r` (constant reaction), with `Ca²` from radial equilibrium;
   - `exponential` (Dixon zero-power): `Cw₁ = a − b/r`, `Cw₂ = a + b/r` (constant work), with `Ca²` from radial equilibrium.
   Stage 1 uses the simplified assignment `α₁ = 0`, uniform `Ca`, `U = r·ω`.
6. **Self-check (lines 358–372).** Mean-row angles from the radial solution are compared with the mean-line solution; tolerance 0.05°.

---

## Outputs

Running `LPC_RUN1.m` opens **13 output windows**: seven plots and six tabular windows. The saved copies are in `outputs/figures/`.

| # | File | Content |
|---:|---|---|
| 1 | `outputs/figures/fig01_annulus_meridional_section.fig` | Tip/mean/root radii from inlet to stage 7 exit |
| 2 | `outputs/figures/fig02_de_haller_and_deflection.fig` | De Haller number per stage with 0.72 limit, and rotor deflection |
| 3 | `outputs/figures/fig03_velocity_triangles.fig` | Velocity triangles for all seven stages |
| 4 | `outputs/figures/fig04_air_angles_vs_radius.fig` | Absolute and relative angles at root, mean and tip |
| 5 | `outputs/figures/fig05_stage_coefficients.fig` | Work-done factor, reaction, loading and flow coefficients |
| 6 | `outputs/figures/fig06_mean_line_stages_table.fig` | Ten-column mean-line summary, seven stages |
| 7 | `outputs/figures/fig07_stage1_root_mean_tip_table.fig` | Stage 1 root/mean/tip data |
| 8 | `outputs/figures/fig08_exit_annulus_table.fig` | Exit annulus geometry per stage |
| 9 | `outputs/figures/fig09_all_stages_root_mean_tip_table.fig` | All stages at root, mean and tip |
| 10 | `outputs/figures/fig10_inlet_summary_table.fig` | Inlet calculation summary |
| 11 | `outputs/figures/fig11_exponential_vs_first_power_angles.fig` | Blade-angle comparison between vortex laws |
| 12 | `outputs/figures/fig12_de_haller_root_tip_comparison.fig` | De Haller at root and tip for both laws |
| 13 | `outputs/figures/fig13_exponential_vs_first_power_table.fig` | Tabulated vortex-law comparison |

The `.fig` files are MATLAB figure containers written by the scripts; they are outputs, not inputs. `outputs/data/PROTO.mat` is a saved workspace snapshot. `outputs/reports/LPC_Modelling_report.pdf` is a scanned/image-based design report — no selectable text could be extracted from it, so it is not used as a machine-readable source here.

---

## Script classification

| Path | Status | What it actually is |
|---|---|---|
| `LPC_RUN1.m` | **Primary** | Complete seven-stage LPC mean-line workflow, three vortex laws, radial analysis, self-check, 13 outputs. Use this script. |
| `development/TESTLPC1.m` | Historical (235 lines) | Simplified early variant: first-power law only, 3 plots, approximate De Haller (`cos β₁/cos β₂`), no tables or self-check. |
| `development/TESTLPC2.m` | Historical (248 lines) | `development/TESTLPC1.m` with radicand guards against imaginary axial velocity; otherwise the same early method. |
| `extended_study/TESTHPC1.m` | Extended study, **not a verified HPC model** (644 lines) | Nearly identical to `LPC_RUN1.m` but at a different operating point (`N = 8382 rpm`, `T₀₁ = 610.294 K`, `P₀₁ = 289.988 kPa`, `PR = 2.9`). It retains the LPC header, the LPC mass flow of 131.2447 kg/s and the seven-stage LPC architecture, and adds Lieblein/Howell/tip-clearance loss models and a Design-A/B comparison. Treat it as an exploratory redesign study using HPC-like inlet conditions, not as a validated high-pressure compressor. |
| `extended_study/LPC_Experimental.mlx` | Live Script counterpart | MATLAB Live Script containing essentially the same algorithm plus the Design-A/B loss section, with saved outputs. Binary OOXML container; it is a source file, not a generated figure. |

The earlier README's statement that "LPC and HPC compressor analysis [is] complete" is not supported: the HPC-labelled script does not implement an independently validated HPC architecture. No turbine source file exists, so the former "turbine analysis ongoing" note is now recorded as **not started in this repository**.

---

## Running the analysis

Requirements: base MATLAB (the scripts use no additional toolboxes; an earlier README recommended R2020b or later). There is **no Simulink model in this repository** — the previous instruction to open `ElectroMechanics_GIT.slx` referred to a file that does not exist here.

```matlab
% From the repository root:
run('LPC_RUN1.m')
```

Edit the `INPUT` block at the top of `LPC_RUN1.m` to change operating conditions. To exercise a different radial law, set `VORTEX_LAW` to `'free'`, `'first_power'`, or `'exponential'`. The historical scripts in `development/` and the extended study in `extended_study/` are run the same way and are not required for the primary analysis.

The scripts write figures to the MATLAB figure windows; they do not overwrite the saved `.fig` files under `outputs/figures/` automatically.

---

## Assumptions and limitations

1. Mean-line, ideal-gas, polytropic model; no blade-profile or viscous-loss model in the primary script.
2. Constant mean radius; annulus computed at stator exit.
3. Stage 1 radial distribution is simplified, not a full radial-equilibrium solution.
4. Loss models (Lieblein profile, Howell secondary, tip clearance) exist only in the extended study and are not validated here.
5. No comparison against published per-stage Olympus 593 data is included.
6. No turbine model, no off-design points, no CFD, no physical test data.
7. The PDF notes cannot be searched or extracted as text; numeric provenance rests on the MATLAB code and saved figures.
8. Figure counts differ by script: 13 for `LPC_RUN1.m`, 3 for each historical `TESTLPC` script, and 17 for `extended_study/TESTHPC1.m` / the Live Script.

---

## Repository layout

| Path | Contents |
|---|---|
| `LPC_RUN1.m` | Primary LPC mean-line design script |
| `development/` | Historical LPC iterations (`development/TESTLPC1.m`, `development/TESTLPC2.m`) |
| `extended_study/` | HPC-labelled variant and Live Script with loss models |
| `outputs/figures/` | Thirteen saved MATLAB figure containers |
| `outputs/data/` | Saved workspace snapshot (`outputs/data/PROTO.mat`) |
| `outputs/reports/` | Image-based design report PDF |

## References

- Cohen, H., Rogers, G. F. C., Saravanamuttoo, H. I. H., *Gas Turbine Theory*, 5th ed. — mean-line compressor analysis and stage stacking.
- Dixon, S. L., *Fluid Mechanics and Thermodynamics of Turbomachinery* — free-vortex, first-power and exponential vortex laws.
- Lieblein, S. — profile loss correlation (used only in `extended_study/`).
- Howell, A. R. — secondary-flow loss correlation (used only in `extended_study/`).

The engine parameters are published-condition values entered in the script; this repository does not contain a component datasheet or rig measurement.
