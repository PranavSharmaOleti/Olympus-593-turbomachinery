# Reverse-Engineering of the Rolls-Royce Olympus 593 LPC
### Stage-by-Stage Mean-Line Compressor Design | MATLAB

# For Handwritten Notes refer the Mathematical Model Document
---

## What This Is

The Rolls-Royce Olympus 593 was the engine that powered the Concorde — the only supersonic commercial aircraft to enter service. This project reverse-engineers the aerodynamic design of its 7-stage Low Pressure Compressor (LPC) using publicly available engine performance data.

Starting from datasheet parameters — mass flow, shaft speed, inlet conditions, and overall pressure ratio — the code reconstructs the internal stage-by-stage design: how fast the blades spin, what angles the air flows through, how much pressure each stage adds, and how the annulus geometry changes from front to back.

The method follows Cohen, Rogers & Saravanamuttoo's *Gas Turbine Theory*, the standard reference for axial compressor mean-line analysis.

> **Status:** LPC and HPC compressor analysis complete. Turbine analysis ongoing.

---

## Engine Parameters Used

| Parameter | Value |
|---|---|
| Engine | RR Olympus 593 Mk.531 (LPC) |
| Mass flow | 131.24 kg/s |
| Shaft speed | 5,819 RPM |
| Inlet total temperature | 389.86 K |
| Inlet total pressure | 73.66 kPa |
| Inlet Mach number | 0.549 |
| Overall pressure ratio | 4.1 |
| Polytropic efficiency | 0.8782 |
| Number of stages | 7 |

---

## What the Code Does

**Step 1 — Inlet**
Computes the inlet flow conditions from the given Mach number: static temperature, density, speed of sound, and axial velocity. Derives the annulus geometry from continuity — tip radius, hub radius, mean radius — and computes blade speed from the shaft RPM. Flow coefficient is an output, not an input.

**Step 2 — Stage Loading**
Temperature rise for stages 3–7 is solved automatically so that all rear stages carry equal load and the overall pressure ratio is met exactly. Stages 1 and 2 are set manually (lower loading, since the stage-1 blade tip is already near transonic).

**Step 3 — Mean-Line Velocity Triangles**
For each stage, computes the full velocity triangle: absolute and relative air angles at rotor inlet and exit, whirl velocities, stage reaction, blade loading coefficient, De Haller number, and rotor deflection. Uses polytropic — not isentropic — efficiency throughout.

**Step 4 — Annulus Geometry**
Computes the exit annulus area, blade height, hub and tip radii at the stator exit of every stage. Mean radius is held constant (Cohen's constant-mean-radius assumption); hub-to-tip ratio rises as density increases toward the rear — expected behaviour.

**Step 5 — Radial Analysis (Hub to Tip)**
Extends the mean-line solution to root, mean, and tip sections using three radial vortex laws:
- **Free vortex** — whirl velocity × radius = constant, axial velocity uniform
- **First-power (Dixon)** — constant reaction from hub to tip
- **Exponential (Dixon zero-power)** — constant work, milder blade twist than free vortex

All three are cross-validated against each other. Mean-line values are self-checked to ensure consistency.

---

## Outputs

Running the script opens **13 figure windows** in MATLAB:

| # | Output |
|---|---|
| 1 | Annulus profile — tip, mean, and hub radii from inlet to stage 7 exit |
| 2 | De Haller number per stage — with 0.72 stall limit marked |
| 3 | Rotor deflection per stage |
| 4 | Velocity triangles — all 7 stages plotted |
| 5 | Air angles vs. stage number — root, mean, tip |
| 6 | Stage coefficients — loading ψ, reaction Λ, work-done factor λ, flow coefficient φ |
| 7–10 | Tabular summaries — mean-line results, stage 1 root/mean/tip, exit annulus, all stages |
| 11 | Inlet summary table |
| 12–13 | Exponential vs. first-power comparison — angles and De Haller at root and tip |

---

## How to Run

1. Open `ElectroMechanics_GIT.slx` in MATLAB (R2020b or later recommended)
2. Edit the `INPUT` block at the top if needed — all tunable parameters are there
3. Press **Run**
4. 13 figure windows will open; use the Windows menu in MATLAB to navigate between them

No additional toolboxes required beyond base MATLAB.

---

## Key Design Choices

- **One blade speed** — U is derived from shaft RPM and the inlet mean radius. It is not a free parameter.
- **Polytropic efficiency** — used throughout for pressure ratio calculations, not isentropic, which would underestimate work in a real multi-stage machine.
- **Stage 1 has no IGV** — stage reaction is an output for stage 1, not prescribed.
- **Vortex law is user-selectable** — switch between `'free'`, `'first_power'`, or `'exponential'` in the INPUT block.

---

## References

- Cohen, H., Rogers, G.F.C., Saravanamuttoo, H.I.H. — *Gas Turbine Theory*, 5th ed.
- Dixon, S.L. — *Fluid Mechanics and Thermodynamics of Turbomachinery*
- Publicly available Olympus 593 performance datasheet (Concorde technical documentation)

---

## Author

**Oleti Pranava Sharma**  
B.Tech Mechanical Engineering — NIT Durgapur  
[GitHub](https://github.com/PranavSharmaOleti)
