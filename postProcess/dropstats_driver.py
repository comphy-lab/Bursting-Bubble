#!/usr/bin/env python3
"""dropstats_driver.py — a case's snapshots to a tidy, tracked drop table.

`getDropStats` is memoryless: it reduces one snapshot and knows nothing about
which fragment is which. Everything that needs memory lives here — drop
identity across snapshots, the emission index `n` that Fig-11-style plots are
drawn against, the rule that decides *when* a drop's radius and velocity are
read, and the ascending-drop filter used for the emitted totals.

The measurement rule matters more than it looks. A fragment's volume is not
meaningful in the snapshot where it detaches: it is still exchanging mass with
the ligament, and the connected-component count flickers as the neck thins
through the `f = 0.5` contour. So a track is only measured once its volume has
settled, and never after it re-merges or leaves the domain.

Usage
-----
    dropstats_driver.py --case simulationCases/1101 \\
        --tool postProcess/getDropStats [--jobs 8] [--ztop 4.0]

Writes, inside the case directory:

- `dropstats.csv`      every drop sample, with a stable `track` and its `n`
- `dropstats_main.csv`  the main body and jet tip, one row per snapshot
- `dropstats_summary.json`  first-drop radius and velocity, drop count, and
  the emitted surface / volume / kinetic-energy totals

Radii are in `R_0`, velocities in the inertio-capillary `V_c`. The experimental
viscous-capillary velocity is `V/V_mu = Oh * (V/V_c)`; `Oh` is not known to
this script, so the conversion is left to the plotting step.
"""

import argparse, csv, glob, json, math, os, re, subprocess, sys
from concurrent.futures import ThreadPoolExecutor

SNAP_RE = re.compile(r"snapshot-([0-9.]+)$")

# A track is continued if the candidate lies within MATCH_DZ of where the
# previous sample would have drifted to, and its volume is within MATCH_VOL.
# Both are deliberately loose: a mismatch splits one drop into two tracks,
# which is visible in the output, whereas a too-tight tolerance silently
# renumbers every drop and corrupts the emission index.
MATCH_DZ = 0.25
MATCH_VOL = 3.0
SETTLE_TOL = 1e-2       # |dV/dt| / V, in inverse capillary times
SPHERICITY_WARN = 1.15  # Rs/Rv above this is a ligament, not a drop
MIN_CELLS = 8.0         # cells per drop radius below which the radius is the mesh


def snapshots(case):
    out = []
    for p in glob.glob(os.path.join(case, "intermediate", "snapshot-*")):
        m = SNAP_RE.search(os.path.basename(p))
        if m:
            out.append((float(m.group(1)), p))
    return sorted(out)


def run_tool(tool, path):
    """One snapshot -> (main row dict, [drop row dicts]). Tool writes stderr."""
    r = subprocess.run([tool, path], capture_output=True, text=True)
    main, drops = None, []
    for line in r.stderr.splitlines():
        f = line.split()
        if not f:
            continue
        if f[0] == "MAIN" and len(f) >= 19:
            main = dict(t=float(f[1]), nliq=int(f[2]), ndrop=int(f[3]),
                        V=float(f[4]), Rv=float(f[5]), S=float(f[6]),
                        Rs=float(f[7]), zc=float(f[8]), zmin=float(f[9]),
                        zmax=float(f[10]), rmax=float(f[11]), vz=float(f[12]),
                        vr=float(f[13]), Ek=float(f[14]), ztip=float(f[15]),
                        rtip=float(f[16]), vtip=float(f[17]), Vtot=float(f[18]))
        elif f[0] == "DROP" and len(f) >= 14:
            # `dmin`/`cells` are trailing columns added later; tolerate their
            # absence so tables written by the earlier binary still parse.
            drops.append(dict(t=float(f[1]), id=int(f[2]), V=float(f[3]),
                              Rv=float(f[4]), S=float(f[5]), Rs=float(f[6]),
                              zc=float(f[7]), zmin=float(f[8]), zmax=float(f[9]),
                              rmax=float(f[10]), vz=float(f[11]),
                              vr=float(f[12]), Ek=float(f[13]),
                              dmin=float(f[14]) if len(f) > 14 else float("nan"),
                              cells=float(f[15]) if len(f) > 15 else float("nan")))
    if main is None:
        print(f"WARNING: no MAIN row from {path}", file=sys.stderr)
    return main, drops


def track(frames):
    """Greedy nearest-neighbour linking of drops across consecutive frames.

    Component ids from `tag()` are not stable between snapshots — the same
    physical drop is routinely relabelled — so linking is by predicted
    position and volume, not by id.
    """
    tracks, live = [], {}          # track index -> last sample
    for (t, drops), (tprev, _) in zip(frames, [(None, None)] + frames[:-1]):
        dt = 0.0 if tprev is None else t - tprev
        taken, assigned = set(), {}
        for d in drops:
            best, best_cost = None, None
            for ti, last in live.items():
                if ti in taken:
                    continue
                zpred = last["zc"] + last["vz"] * dt
                dz = abs(d["zc"] - zpred)
                vr = d["V"] / last["V"] if last["V"] > 0 else math.inf
                if dz > MATCH_DZ or vr > MATCH_VOL or vr < 1.0 / MATCH_VOL:
                    continue
                cost = dz + abs(math.log(vr))
                if best_cost is None or cost < best_cost:
                    best, best_cost = ti, cost
            if best is None:                       # a new drop is born here
                best = len(tracks)
                tracks.append([])
            taken.add(best)
            assigned[best] = d
        for ti, d in assigned.items():
            d = dict(d, track=ti)
            tracks[ti].append(d)
            live[ti] = d
        for ti in list(live):                      # merged, or left the domain
            if ti not in assigned:
                del live[ti]
    return tracks


def measure(tr):
    """First sample at which the track's volume has settled.

    Returns None for a track that never settles — a fragment that merges back
    or leaves the domain within a couple of snapshots is not a measured drop.
    """
    for a, b in zip(tr, tr[1:]):
        dt = b["t"] - a["t"]
        if dt <= 0 or b["V"] <= 0:
            continue
        if abs(b["V"] - a["V"]) / (dt * b["V"]) < SETTLE_TOL:
            return b
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--case", required=True, help="case directory (holds intermediate/)")
    ap.add_argument("--tool", required=True, help="compiled getDropStats binary")
    ap.add_argument("--jobs", type=int, default=os.cpu_count() or 4)
    ap.add_argument("--ztop", type=float, default=4.0,
                    help="top boundary; drops within one radius of it have left")
    ap.add_argument("--min-cells", type=float, default=MIN_CELLS,
                    help="cells per radius below which a drop is treated as unresolved")
    a = ap.parse_args()

    snaps = snapshots(a.case)
    if not snaps:
        sys.exit(f"no snapshots under {a.case}/intermediate")
    print(f"{len(snaps)} snapshots, {a.jobs} workers", file=sys.stderr)

    with ThreadPoolExecutor(max_workers=a.jobs) as ex:
        results = list(ex.map(lambda tp: run_tool(a.tool, tp[1]), snaps))

    mains = [m for m, _ in results if m]
    frames = [(m["t"], d) for (m, d) in results if m]
    tracks = track(frames)

    # Emission index: order tracks by the time they first appear.
    order = sorted(range(len(tracks)), key=lambda i: tracks[i][0]["t"] if tracks[i] else math.inf)
    nof = {ti: k + 1 for k, ti in enumerate(order) if tracks[ti]}

    rows = []
    for ti, tr in enumerate(tracks):
        for d in tr:
            rows.append(dict(d, n=nof.get(ti, 0),
                             sphericity=d["Rs"] / d["Rv"] if d["Rv"] else 0.0))
    rows.sort(key=lambda r: (r["t"], r["n"]))

    with open(os.path.join(a.case, "dropstats.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["t", "track", "n", "id", "V", "Rv", "S", "Rs",
                                           "sphericity", "cells", "dmin", "zc", "zmin",
                                           "zmax", "rmax", "vz", "vr", "Ek"])
        w.writeheader()
        w.writerows(rows)

    with open(os.path.join(a.case, "dropstats_main.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(mains[0].keys()))
        w.writeheader()
        w.writerows(mains)

    # Per-drop measured values, and the emitted totals over ASCENDING drops
    # only — the experimental counterpart discards drops that fall back.
    # `n` in the CSV is the raw track order, transients included. The emission
    # index reported here counts only drops that actually settled, so that a
    # ligament which briefly detaches and merges back does not shift every
    # subsequent drop's index in a Fig-11-style plot.
    measured = []
    for ti in order:
        if not tracks[ti]:
            continue
        m = measure(tracks[ti])
        if m is None:
            continue
        cells = m.get("cells", float("nan"))
        measured.append(dict(track=ti, n=len(measured) + 1, n_raw=nof[ti],
                             t=m["t"], Rd=m["Rv"], S=m["S"], V=m["V"],
                             vz=m["vz"], Ek=m["Ek"], cells=cells,
                             resolved=not (cells == cells and cells < a.min_cells),
                             sphericity=m["Rs"] / m["Rv"] if m["Rv"] else 0.0,
                             left_domain=m["zc"] > a.ztop - m["Rv"]))
    rising = [m for m in measured if m["vz"] > 0.0]
    # Aggregates are reported over RESOLVED rising drops. An unresolved speck
    # contributes a mesh-determined volume and area to N, S_t, V_t and E_kt, so
    # including it silently corrupts every Fig-10 observable.
    rising_resolved = [m for m in rising if m["resolved"]]
    Sb, Vb = 4.0 * math.pi, 4.0 * math.pi / 3.0   # R_0 = 1, so E_sb = sigma*Sb = Sb

    # Two distinct "first drops", reported separately rather than conflated.
    # `first_detached` is the first fragment to leave the body, which at low Oh
    # is routinely a sub-grid speck shed by the retracting film. The quantity
    # the experiment reports is the first drop OF THE WORTHINGTON JET, and the
    # nearest defensible proxy is the first *resolved* rising drop. Selecting
    # the chronologically-first detachment silently measures the wrong object.
    first = rising_resolved[0] if rising_resolved else None
    first_detached = rising[0] if rising else None
    summary = dict(
        case=os.path.abspath(a.case),
        n_snapshots=len(snaps),
        n_tracks=len(tracks),
        n_measured=len(measured),
        n_rising=len(rising),
        min_cells=a.min_cells,
        n_rising_resolved=len(rising_resolved),
        first_detached=None if first_detached is None else dict(
            t=first_detached["t"], Rd_over_R0=first_detached["Rd"],
            vz_over_Vc=first_detached["vz"], cells_per_radius=first_detached["cells"],
            resolved=first_detached["resolved"],
            note="first fragment to detach; may be a sub-grid speck, not the jet drop"),
        first_drop=None if first is None else dict(
            t=first["t"], Rd_over_R0=first["Rd"], vz_over_Vc=first["vz"],
            sphericity=first["sphericity"], cells_per_radius=first["cells"],
            resolved=first["resolved"], n_among_all_rising=first["n"],
            ligament_warning=first["sphericity"] > SPHERICITY_WARN,
            definition="first RESOLVED rising drop (>= min_cells per radius)"),
        emitted_totals_rising_resolved=dict(
            N=len(rising_resolved),
            St_over_Sb=sum(m["S"] for m in rising_resolved) / Sb,
            Vt_over_Vb=sum(m["V"] for m in rising_resolved) / Vb,
            Ekt_over_Esb=sum(m["Ek"] for m in rising_resolved) / Sb),
        emitted_totals_rising_all=dict(
            N=len(rising),
            St_over_Sb=sum(m["S"] for m in rising) / Sb,
            Vt_over_Vb=sum(m["V"] for m in rising) / Vb,
            Ekt_over_Esb=sum(m["Ek"] for m in rising) / Sb),
        # Vtot - V is liquid outside the main body, i.e. the detached volume.
        # It grows as drops are emitted and is physics, not error. It was once
        # named mass_closure_max_residual, which invited exactly the wrong
        # reading: an absolute 0.04 on a pool of 1880 looks alarming next to a
        # relative tolerance and is in fact 2e-5.
        detached_volume_max_abs=max(
            abs(m["Vtot"] - m["V"]) for m in mains) if mains else None,
        # The real closure test: drift of the total liquid volume against its
        # initial value. Nothing physical removes liquid from the domain here,
        # so any drift is solver or reduction error.
        mass_closure_max_rel=(
            max(abs(m["Vtot"] - mains[0]["Vtot"]) for m in mains)
            / mains[0]["Vtot"]) if mains and mains[0]["Vtot"] else None,
        n_transient_tracks=len(tracks) - len(measured),
        per_drop=[dict(n=m["n"], track=m["track"], Rd_over_R0=m["Rd"],
                       vz_over_Vc=m["vz"], sphericity=m["sphericity"],
                       cells_per_radius=m["cells"], resolved=m["resolved"],
                       rising=m["vz"] > 0.0) for m in measured],
    )
    with open(os.path.join(a.case, "dropstats_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2)

    print(json.dumps(summary["first_drop"], indent=2), file=sys.stderr)

    print(f"N(rising)={len(rising)}  N(rising,resolved)={len(rising_resolved)}  "
          f"tracks={len(tracks)}  measured={len(measured)}", file=sys.stderr)
    if first is not None and not first["resolved"]:
        print(f"WARNING: first drop is {first['cells']:.1f} cells per radius "
              f"(< {a.min_cells}); its radius is mesh-limited, not physical.",
              file=sys.stderr)


if __name__ == "__main__":
    main()
