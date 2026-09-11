NGC 5139 (Omega Cen) — two 3D-plot filtered files
===================================================
Both contain the same 50,865 high-confidence members (membership_prob > 0.90).
They differ ONLY in how the depth (Z) axis is defined.

ngc5139_3d_MODELED.csv
  X_pc, Y_pc      : real sky positions, converted to parsecs at d = 5.2 kpc
  Z_modeled_pc    : MODELED depth (spherical symmetry) — a stand-in, not a
                    measurement. Makes the cloud a true-shaped sphere.
  assumed_dist_pc : 5200 (the single literature distance used for all stars)

ngc5139_3d_RAW.csv
  parallax        : Gaia's raw per-star parallax (the colleague's method)
  dist_pc_raw     : 1000/parallax — WARNING: 28.5% of stars are negative,
                    only 2.4% have S/N >= 5. This column is mostly noise.

Which is usable for V/sigma?
----------------------------
NEITHER depth column is used for V/sigma. V/sigma is computed from the
PROPER MOTIONS (pmra, pmdec — carried in both files), and it is
DISTANCE-INDEPENDENT: converting mas/yr to km/s multiplies both V and
sigma by the same factor (4.74 * distance), which cancels in the ratio.

So for V/sigma:
  - Use either file (the PM columns are identical). Depth is irrelevant.
  - Do NOT feed raw parallax / dist_pc_raw into any velocity — it's noise.
  - If you ever need V or sigma in physical km/s SEPARATELY, use the single
    cluster distance (5.2 kpc, the MODELED file's assumed_dist_pc) for every
    star — never per-star parallax.

Columns for V/sigma (both files):
  dpmra, dpmdec   : proper motion minus the cluster's systemic PM
  pm_tan          : tangential component (the ROTATION signal)
  pm_rad          : radial component
  R_arcmin, theta_deg : position for binning by radius / azimuth

Note: V/sigma is a POPULATION quantity (global, or a profile vs radius) —
not a per-star value. And a correct sigma must be de-convolved for Gaia's
measurement errors; a raw standard deviation of pmra/pmdec overestimates
sigma (it includes the error scatter), which makes an uncorrected V/sigma
come out far too low.
