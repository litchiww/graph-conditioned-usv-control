[README.md](https://github.com/user-attachments/files/30381623/README.md)
# Anonymized Reproducibility Package

## Associated manuscript

Graph-Conditioned Barrier-Nussbaum Control With Local Direction Handoff for Cooperative USVs

This package contains the MATLAB source, deterministic seeds, final audit outputs, canonical case files, and plotting scripts used for the reported results. It is anonymized for double anonymized peer review.

## Requirements

- MATLAB R2024a or later
- Base MATLAB for simulation, data processing, and plotting
- Parallel Computing Toolbox is optional and only accelerates the exhaustive sign and long-horizon audits
- Sufficient free disk space for regenerated case and figure files

No external MATLAB packages are required.

## Quick verification

Extract the archive, open MATLAB in the extracted package directory, and run:

```matlab
VERIFY_PACKAGE
```

This checks the packaged audit tables, verifies the headline numerical results, and runs a short simulation through the final controller implementation.

Expected summary:

- 64 of 64 fixed sign assignments accepted
- worst fixed-sign internal-envelope ratio below 0.857
- 30 of 30 randomized sign, gain-magnitude, and initial-state trials accepted
- worst randomized internal-envelope ratio below 0.919
- 5 of 5 long-horizon audits accepted
- offline 3-DOF second-half lateral RMS of approximately 0.530488 m
- offline 3-DOF heading-reference RMS of approximately 1.82846 degrees
- 25-degree rudder-equivalent limit reached on 110 of 43203 agentwise
  20-Hz control updates, or approximately 0.254612 percent
- direct 3-DOF yaw-loop lateral RMS of approximately 2.37461 m
- direct 3-DOF yaw-loop heading RMS of approximately 5.92920 degrees
- direct 3-DOF internal-envelope ratio below 0.864 with zero yaw-moment
  saturation in 21603 agentwise 20-Hz control updates

## Reproduction entry points

```matlab
REPRODUCE('verify')
REPRODUCE('figures')
REPRODUCE('canonical')
REPRODUCE('closedloop')
REPRODUCE('signs')
REPRODUCE('random')
REPRODUCE('solver')
REPRODUCE('long')
REPRODUCE('audits')
REPRODUCE('all')
```

`figures` regenerates all manuscript plots from the packaged case and audit
data. `closedloop` reruns only the direct 3-DOF yaw-loop test. `canonical`
reruns the nominal, ablation, topology, sea-state, graph-conditioning, offline
3-DOF replay, and direct 3-DOF yaw-loop studies before plotting. `audits`
reruns the exhaustive fixed-sign, randomized, solver-refinement, and five
600-s checks. `all` executes both groups and can require several hours.

When Parallel Computing Toolbox is unavailable, the fixed-sign audit automatically uses the sequential implementation.

## Directory structure

```text
MAIN.m
REPRODUCE.m
VERIFY_PACKAGE.m
run_all.m
src/
cases/
data/
MANIFEST.txt
COMMENT_REMOVAL_REPORT.txt
```

`src` contains the final controller, simulation, audit, metric, and plotting
functions. `cases` contains plotting-resolution copies of the canonical cases
at a 0.05-s saved interval, with the original run metrics and parameter
structures preserved. These files support immediate figure regeneration
without altering any reported statistic. Full-resolution cases can be
regenerated with `REPRODUCE('canonical')`. `data` contains compact
machine-readable outputs for the analytical certificate, exhaustive and
randomized direction audits, long-horizon checks, solver refinement,
graph-conditioning scan, ablations, the offline 3-DOF replay, and the direct
3-DOF yaw-loop closure.

## Numerical scope

The normalized closed-loop simulations support the theorem-domain and finite
verified-domain claims stated in the manuscript. Two distinct 3-DOF paths are
included. The offline sampled autopilot replay smooths the saved C2 lane
trajectory on the physical time grid and converts it to kinematically
consistent heading and yaw-rate references. It includes actuator limits,
sensor noise, and packet loss, but it neither inherits the normalized `x2`
envelope nor constitutes end-to-end validation of the theorem controller.

The direct yaw-loop path maps heading and yaw-rate errors to the theorem
coordinates and applies the controller output directly to nonidentical signed
physical yaw-moment channels at 20 Hz. A finite, zero-mean commissioning probe
is used only until the local direction handoff locks. That initialization is
an empirical implementation condition and is not presented as an extension of
the theorem. The direct path validates the yaw-loop realization, not field
operation, collision avoidance, or the complete outer-loop inspection system.

The demand-dependent envelope is an internal operational envelope. It is not an externally certified collision-avoidance or physical safe set.

## Anonymization and source comments

Author names, affiliations, email addresses, local absolute paths, development logs, tuning traces, reference papers, and version-control metadata are excluded. MATLAB comments are removed from every packaged `.m` file. Character data containing percent signs, such as format strings and axis labels, are retained because they are executable content rather than comments.

## Integrity

`MANIFEST.txt` provides a SHA-256 digest for every packaged file except the manifest itself. `COMMENT_REMOVAL_REPORT.txt` records the number of MATLAB comment tokens removed from each source file and confirms that no comment tokens remain.
