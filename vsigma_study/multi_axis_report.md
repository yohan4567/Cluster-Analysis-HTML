# 多旋轉軸 (Multiple / Nested Rotation Axes) — Findings

Two testable versions of the idea (see METHODS.md §5):
1. **Sign/amplitude structure of V_rot(R)** — does the line-of-sight spin
   component change with radius (e.g. a counter-rotating core)?
2. **3D spin-axis tilt between radial shells** — with Gaia radial velocities,
   fit the full spin vector per shell; the angle between shells is the
   'nested axis' tilt.

## Per-cluster results

| Cluster | rotation detected (PM) | counter-rotation | N_RV | 3D axis incl (deg) | axis PA (deg) | shell twist ± err (deg) | twist significant? |
|---|---|---|---|---|---|---|---|
| NGC5139 | YES | no | 951 | 144 | 201 | 11 ± 9 | no |
| NGC104 | YES | no | 660 | 20 | 301 | 21 ± 11 | no |
| NGC7078 | YES | no | 57 | 168 | 260 | — | — |
| NGC7089 | no | no | 50 | 44 | 230 | — | — |
| NGC5904 | YES | no | 101 | 150 | 119 | — | — |
| NGC6656 | YES | no | 421 | 153 | 89 | 25 ± 16 | no |
| NGC6273 | YES | no | 79 | 141 | 249 | — | — |
| NGC5272 | YES | no | 82 | 17 | 117 | — | — |
| NGC6752 | YES | no | 216 | 20 | 341 | 49 ± 40 | no |
| NGC6809 | YES | no | 55 | 125 | 291 | — | — |
| NGC288 | no | no | — | — | — | — | — |
| NGC362 | no | no | — | — | — | — | — |
| NGC6397 | no | no | 160 | 91 | 133 | 96 ± 48 | no |
| NGC6341 | YES | no | 62 | 55 | 135 | — | — |
| NGC5024 | no | no | — | — | — | — | — |
| NGC2419 | no | no | — | — | — | — | — |
| NGC1851 | no | no | — | — | — | — | — |
| NGC6266 | YES | no | 192 | 148 | 168 | 47 ± 25 | no |
| Terzan5 | no | no | — | — | — | — | — |
| NGC6440 | no | no | — | — | — | — | — |
| NGC5286 | no | no | — | — | — | — | — |
| NGC6121 | YES | no | 89 | 105 | 259 | — | — |
| Pal5 | no | no | — | — | — | — | — |
| NGC6544 | no | no | — | — | — | — | — |

**How to read this:** a large shell twist with small errors means the inner and
outer parts of the cluster do NOT share one rotation axis — the honest version
of the team's nested-axes idea. A 'YES' in counter-rotation means the
line-of-sight spin actually reverses sign with radius.

**Caveats:** RV samples are small (only bright stars have Gaia RVs); the
perspective-rotation correction has been applied (it would otherwise fake a
rotation gradient of up to ~1.6 km/s for nearby fast-moving clusters).
Stars do not literally orbit 'sub-poles' (the moon-earth analogy) — relaxation
destroys such hierarchies — but a radius-dependent spin vector is real physics
and is exactly what these two tests measure.