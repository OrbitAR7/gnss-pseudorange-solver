# gnss-pseudorange-solver

[![CI](https://github.com/OrbitAR7/gnss-pseudorange-solver/actions/workflows/CI.yml/badge.svg)](https://github.com/OrbitAR7/gnss-pseudorange-solver/actions/workflows/CI.yml)

A header-only C++17 GNSS single-point positioning solver. It estimates receiver position and clock bias from pseudoranges by Gauss–Newton least squares (optionally weighted, including by satellite elevation), with Earth-rotation (Sagnac) correction, DOP and WGS-84 coordinate conversion. It is verified by a test suite against known truth (see [Validation](#validation)).

## Quick example

```cpp
#include "gnss_solver.hpp"
using namespace gnss;

int main() {
    // True receiver: Riyadh, with a 1.5 km receiver clock bias (expressed in metres)
    LLA truth_lla; truth_lla.lat = 24.7 * constants::DEG2RAD; truth_lla.lon = 46.7 * constants::DEG2RAD; truth_lla.alt = 600.0;
    ECEF truth = llaToECEF(truth_lla);

    std::vector<SatelliteData> sats = /* satellite ECEF positions, e.g. from broadcast ephemeris */;

    // Simulated pseudoranges with 1 m noise (the simulator does not model Earth rotation)
    auto meas = simulatePseudoranges(truth, 1500.0, sats, 1.0);

    SolverConfig cfg;
    cfg.apply_earth_rotation = false;     // match the simulator; keep true for real measurements
    NavSolution sol = GNSSSolver(cfg).solve(sats, meas);

    printSolution(sol, "Riyadh test");
}
```

`tests/test_gnss_solver.cpp` contains a complete, runnable version that builds a satellite geometry around the receiver.

## Build and test

```bash
cmake -S . -B build && cmake --build build
cd build && ctest --output-on-failure   # validation tests (also runs both examples)
./gps_position_fix                       # 24-satellite constellation, 100 noisy trials
./dop_analysis                           # geometry, satellite count, urban canyon, elevation mask
```

Or copy `include/gnss_solver.hpp` into your project. It has no dependencies.

## Validation

`tests/test_gnss_solver.cpp` uses 8 satellites at GPS orbital radius, spread in azimuth and elevation above a receiver in Riyadh:

| Test | Result |
|---|---|
| Noise-free pseudoranges, no Earth rotation | Position recovered to ~5 nm, clock to ~1 nm (5 iterations) |
| Pseudoranges generated with Earth rotation during signal transit; solver's Sagnac correction on | 2 mm position error |
| Same data with the Sagnac correction turned off | 32 m error (the size of the effect) |
| 1 m white noise, 500 Monte Carlo runs | 3-D RMS error 1.71 m vs. PDOP × σ = 1.71 m |
| Elevation weighting on, noise-free | Position recovered to ~1 nm (6 iterations) |
| Noise σ₀ / sin(el) with σ₀ = 1 m, 1000 runs, unweighted | 4.11 m vs. 4.07 m predicted |
| Same data, elevation-weighted | 3.34 m vs. 3.33 m predicted (19% smaller) |

## How it works

Pseudorange model for satellite *i*:

```
ρᵢ = ‖r_sat,i − r_rx‖ + c·δt_rx − c·δt_sat,i + ε
```

Linearize about the current estimate and iterate:

```
δx = (HᵀWH)⁻¹ HᵀW δρ
```

Here **H** holds the negated unit line-of-sight vectors and a column of ones for the clock term. DOPs come from **Q = (HᵀH)⁻¹**, rotated into local east-north-up axes for HDOP/VDOP.

**Elevation weighting:** low satellites pass through more atmosphere and pick up more multipath, so their pseudoranges are noisier. With `elevation_weighting = true` the solver assumes σᵢ = σ₀ / sin(elᵢ) and weights each satellite by sin²(elᵢ), recomputing elevations from the current estimate on every iteration. If that noise model holds, the weighted error covariance is (HᵀWH)⁻¹. The unweighted solution has the larger (HᵀH)⁻¹ HᵀΣH (HᵀH)⁻¹. The tests check both against Monte Carlo runs.

**Earth rotation:** during the ~70 ms signal transit the Earth turns by ω·τ, so each satellite position is rotated about the z-axis by that angle before the range is computed. τ is approximated by ρ/c, which includes the receiver clock bias. That is accurate to centimetres for typical clock biases.

## API summary

| Item | Description |
|---|---|
| `ECEF`, `LLA` | Coordinates (m; rad, rad, m) with `ecefToLLA` / `llaToECEF` (WGS-84) |
| `SatelliteData` | PRN, ECEF position [m], satellite clock bias [s] |
| `Pseudorange` | PRN, range [m], C/N0, weight |
| `SolverConfig` | Convergence threshold, max iterations, `apply_earth_rotation`, `weighted_least_squares`, `elevation_weighting`, initial guess |
| `GNSSSolver::solve` | Returns `NavSolution`: position, clock bias [m], LLA, DOP, iterations, residuals, validity |
| Helpers | `elevationAngle`, `simulatePseudoranges` (5° elevation mask, Gaussian noise, optionally σ / sin(el)), `printSolution` |

## Limitations and next steps

- No ionospheric or tropospheric delay models, and no broadcast-ephemeris or RINEX parsing: satellite positions are inputs.
- Single constellation, code-phase only (no carrier phase).
- The simulator omits Earth rotation; turn off `apply_earth_rotation` when solving simulated data.
- Natural extensions are integrity monitoring (see [`gnss-integrity-raim`](https://github.com/OrbitAR7/gnss-integrity-raim)), Klobuchar and Saastamoinen corrections, and multi-constellation inter-system biases.

## References

- Misra & Enge, *Global Positioning System: Signals, Measurements, and Performance*, 2nd ed., 2006
- Kaplan & Hegarty (eds.), *Understanding GPS/GNSS: Principles and Applications*, 3rd ed., 2017

## License

MIT. See `LICENSE`.
