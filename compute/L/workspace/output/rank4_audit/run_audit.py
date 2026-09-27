#!/usr/bin/env python3
"""Rebuild the rank-4 audit in a fresh isolated project layout."""
from __future__ import annotations
import argparse
import collections
import datetime
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import time

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
DEPENDENCY_HASH = "eed37d27b0950a5c97bd03f59767e116545daf2cec4baaa234c8c3a84e92dc8e"

def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()

def rows(path):
    return [json.loads(s) for s in Path(path).read_text().splitlines() if s.strip()]

def require(ok, message):
    if not ok:
        raise RuntimeError(message)

def flat(targets, path):
    Path(path).write_text("".join(str(t["id"])+" "+" ".join(str(x) for r in t["H"] for e in r for x in e)+"\n" for t in targets))

def gp_matrix(H):
    return "["+";".join(",".join(f"({a})+({b})*w" for a,b in r) for r in H)+"]"

def tool(name, alternatives):
    value = shutil.which(name)
    if value:
        return value
    for value in alternatives:
        if Path(value).is_file():
            return value
    raise RuntimeError(f"Required executable not found: {name}")

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", type=Path, help="New, non-existing directory; never overwrites a prior run.")
    parser.add_argument("--check-d3-inputs", action="store_true", help="Re-enumerate all 30,423 existing exact balls; does not regenerate the d3 reduction.")
    parser.add_argument("--ubsan", action="store_true", help="Also rebuild and rerun the h=3 abstract sieve with undefined-behavior checks.")
    parser.add_argument("--gp", help="PARI/GP executable override.")
    parser.add_argument("--cxx", help="C++17 compiler override.")
    parser.add_argument("--prepare-only", action="store_true", help="Copy inputs and sources, write the manifest, and stop before calculation.")
    args = parser.parse_args()
    gp = args.gp or tool("gp", ["/opt/homebrew/bin/gp", "/usr/local/bin/gp"])
    cxx = args.cxx or tool("clang++", ["/usr/bin/clang++"])
    if args.run_dir:
        root = args.run_dir.expanduser().resolve()
        root.mkdir(parents=True, exist_ok=False)
    else:
        runs = HERE/"runs"
        runs.mkdir(exist_ok=True)
        root = Path(tempfile.mkdtemp(prefix=datetime.datetime.now().strftime("%Y%m%d-%H%M%S-"), dir=runs))
    audit = root/"output/rank4_audit"
    audit.mkdir(parents=True)
    logs = audit/"driver_logs"
    logs.mkdir()
    print(f"RUN_DIRECTORY {root}", flush=True)
    # Only source files are copied from the present audit. Generated JSON,
    # shells and enumeration output are never seeded into the new calculation.
    for source in HERE.iterdir():
        if source.is_file() and source.suffix in {".py", ".cpp", ".gp"}:
            shutil.copy2(source, audit/source.name)
    deps = [
        "rank88_exact_data.gp",
        "output/d3/benchmark/ambient.gp",
        "output/d3/benchmark/exact_enum.cpp",
        "output/d3/benchmark/all_exact_inputs.txt",
        "output/d3/benchmark/all_exact_results.jsonl",
        "output/d3/benchmark/certificate_summary.json",
        "output/direct_extremality/countermodel_data.gp",
        "output/direct_extremality/basis_search_target_result.gp",
        "output/d4/feasibility/free_counterexample.gp",
    ]
    for rel in deps:
        source, dest = PROJECT/rel, root/rel
        require(source.is_file(), f"Missing dependency: {source}")
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source, dest)
    needed = [
        "verify_ambient_input.gp", "independent_chart_audit.gp", "export_charts_json.gp",
        "independent_abstract.cpp", "independent_h2_abstract.py", "h2_export.gp",
        "h2_binary_check.cpp", "h2_candidate_witnesses.py", "independent_M_shell.gp",
        "ambient_qfauto.gp", "check_M_shell_orbits.cpp", "independent_embedding_sieve.cpp",
        "verify_binary_witnesses.py", "check_binary_absence.cpp", "h3_source_orbits.gp",
        "verify_source_orbits.py", "independent_simultaneous.cpp", "verify_simultaneous.py",
        "exact_affine_enum.cpp", "prepare_lowspan_target.py", "affine_embedding_audit.py",
        "affine_metadata_verify.py", "affine_stage_generator.gp", "test_affine_enum.py",
        "embedding_controls.gp", "check_final_coverage.py",
    ]
    require(all((audit/p).is_file() for p in needed), "Some audit sources are missing; finish source preparation before running.")
    manifest = {str(p.relative_to(root)): sha(p) for p in root.rglob("*") if p.is_file()}
    (audit/"INPUT_MANIFEST.json").write_text(json.dumps(manifest, indent=2)+"\n")
    actual_dep = sha(root/"output/d3/benchmark/all_exact_inputs.txt")
    require(actual_dep == DEPENDENCY_HASH, "The existing d2/d3 exact-ball input differs from the audited input.")
    summary = {
        "status": "prepared", "run_directory": str(root), "gp": gp, "compiler": cxx,
        "python": sys.executable, "assertions_enabled": True,
        "dependency": {"input_sha256": actual_dep, "rechecked": False,
                       "scope": "Existing d2/d3 certificate inputs; no d3 reduction regeneration."},
        "stages": [],
    }
    summary_path = audit/"FULL_RUN_SUMMARY.json"
    def save():
        summary_path.write_text(json.dumps(summary, indent=2)+"\n")
    save()
    if args.prepare_only:
        print(f"PREPARED {summary_path}", flush=True)
        return
    env = os.environ.copy()
    env["AUDIT_GP"] = gp
    env["PYTHONUNBUFFERED"] = "1"
    start = time.perf_counter()
    def run(label, command, stdout=None, gp_check=False, timeout=600):
        print(f"START {label}", flush=True)
        ts = time.perf_counter()
        destination = Path(stdout) if stdout else logs/(label+".stdout")
        error = logs/(label+".stderr")
        with destination.open("w") as out, error.open("w") as err:
            p = subprocess.run([str(v) for v in command], cwd=root, env=env, stdout=out, stderr=err, timeout=timeout)
        item = {"stage": label, "seconds": time.perf_counter()-ts, "returncode": p.returncode,
                "stdout": str(destination.relative_to(root)), "stderr": str(error.relative_to(root))}
        summary["stages"].append(item)
        save()
        require(p.returncode == 0, f"{label} failed. See {error}")
        if gp_check:
            msg = destination.read_text(errors="replace")+"\n"+error.read_text(errors="replace")
            require(not any(line.lstrip().startswith("***") and "Warning:" not in line for line in msg.splitlines()), f"PARI reported an error in {label}: {error}")
        print(f"DONE {label} {item['seconds']:.3f}s", flush=True)
    def pyg(name, *extra, stdout=None):
        run(Path(name).stem, [sys.executable, audit/name, *extra], stdout=stdout)
    def gpr(name, label=None):
        run(label or Path(name).stem, [gp, "-fq", audit/name], gp_check=True)
    # Homebrew arm64 and Intel prefixes are detected, while ordinary Unix
    # installations can use the compiler's default include/library paths.
    flags = ["-std=c++17", "-O3"]
    gmpflags = []
    for prefix in ["/opt/homebrew", "/usr/local"]:
        if (Path(prefix)/"include/gmpxx.h").is_file():
            gmpflags += ["-I"+prefix+"/include", "-L"+prefix+"/lib"]
            break
    gmpflags += ["-lgmpxx", "-lgmp"]
    binaries = {
        "independent_abstract": False, "check_M_shell_orbits": False,
        "check_binary_absence": False, "h2_binary_check": True,
        "independent_embedding_sieve": True, "independent_simultaneous": True,
        "exact_affine_enum": True,
    }
    try:
        for name, uses_gmp in binaries.items():
            run("compile_"+name, [cxx, *flags, audit/(name+".cpp"), "-o", audit/name, *(gmpflags if uses_gmp else [])])
        exact = root/"output/d3/benchmark/exact_enum"
        run("compile_exact_enum", [cxx, *flags, root/"output/d3/benchmark/exact_enum.cpp", "-o", exact, *gmpflags])
        run("verify_ambient_input", [gp, "-fq", audit/"verify_ambient_input.gp"], stdout=audit/"ambient_input_verification.log", gp_check=True)
        # Every used M generator must preserve both the trace form and O action.
        (audit/"driver_generator_audit.gp").write_text(
            'read("output/d3/benchmark/ambient.gp");\nread("output/rank4_audit/ambient_qfauto.gp");\n'
            'for(i=1,#AUT[2],A=AUT[2][i];if(A~*G*A!=G || A*W!=W*A || abs(matdet(A))!=1,error("invalid M generator")));\n'
            'print("M_GENERATORS_VERIFIED");\nquit;\n')
        gpr("driver_generator_audit.gp")
        gpr("independent_chart_audit.gp")
        gpr("export_charts_json.gp")
        pyg("independent_h2_abstract.py")
        for h in (2,3):
            run(f"abstract_h{h}", [audit/"independent_abstract", str(h)], stdout=audit/f"independent_h{h}_cpp.log")
        py_h2 = json.loads((audit/"independent_h2_abstract.json").read_text())["candidates"]
        cpp_h2 = rows(audit/"independent_h2_HY.jsonl")
        require({json.dumps(r["HY"]) for r in py_h2} == {json.dumps(r["H"]) for r in cpp_h2}, "Python/C++ h2 abstract candidate sets differ.")
        targets2, targets3 = rows(audit/"independent_h2_targets.jsonl"), rows(audit/"independent_h3_targets.jsonl")
        require(len(targets2)==23 and len(targets3)==13272, "Unexpected abstract counts.")
        if args.ubsan:
            ub = audit/"independent_abstract_ubsan"
            run("compile_abstract_ubsan", [cxx, "-std=c++17", "-O2", "-fsanitize=undefined", "-fno-sanitize-recover=all", audit/"independent_abstract.cpp", "-o", ub])
            run("abstract_h3_ubsan", [ub, "3", "output/rank4_audit/ubsan_h3"])
            require((audit/"ubsan_h3_targets.jsonl").read_bytes()==(audit/"independent_h3_targets.jsonl").read_bytes(), "UBSan output differs.")
        gpr("h2_export.gp")
        for label in ("M", "N"):
            run("h2_binary_"+label, [audit/"h2_binary_check", audit/f"h2_{label}_input.txt"], stdout=audit/f"h2_{label}_result.json")
        mcheck = json.loads((audit/"h2_M_result.json").read_text())
        ncheck = json.loads((audit/"h2_N_result.json").read_text())
        require(mcheck["complete"] and mcheck["obstruction_representations"]==[0,0,0,0], "M binary absence check failed.")
        require(ncheck["complete"] and any(ncheck["obstruction_representations"]), "N positive control failed.")
        pyg("h2_candidate_witnesses.py")
        gpr("independent_M_shell.gp")
        # PARI write(vector) uses brackets and commas; convert only that syntax.
        for name in ("M_integer_data.txt", "M_shell12.txt", "M_shell12_count_input.txt"):
            path = audit/name
            path.write_text(re.sub(r"[\[\],;]", " ", path.read_text()))
        run("exact_M_shell_count", [exact, audit/"M_shell12_count_input.txt", "55"], stdout=audit/"M_shell12_exact_count.json")
        shellcount = json.loads((audit/"M_shell12_exact_count.json").read_text())
        require(shellcount["complete"] and shellcount["pairs"]==283613, "Incomplete M sphere enumeration.")
        run("M_shell_orbits", [audit/"check_M_shell_orbits"], stdout=audit/"M_shell12_orbit_audit.json")
        orb = json.loads((audit/"M_shell12_orbit_audit.json").read_text())
        require(orb["unique_valid"] and orb["generator_closed"] and sum(s["orbits"] for s in orb["shells"])==78, "M shell orbit audit failed.")
        for h, data, name in [(2,targets2,"independent_h2_embeddings"),(3,targets3,"independent_h3_binary_sieve")]:
            f = audit/f"independent_h{h}_targets.flat"
            flat(data,f)
            run(f"binary_sieve_h{h}", [audit/"independent_embedding_sieve", f], stdout=audit/(name+".jsonl"))
        pyg("verify_binary_witnesses.py")
        run("independent_binary_absence", [audit/"check_binary_absence"], stdout=audit/"binary_absence_independent.json")
        require(json.loads((audit/"binary_absence_independent.json").read_text())["hits"]==0, "A forbidden M binary type was found.")
        b3 = rows(audit/"independent_h3_binary_sieve.jsonl")
        byid = {t["id"]:t for t in targets3}
        survivors = [byid[r["id"]] for r in b3 if r["status"]=="survives"]
        require(len(survivors)==1377, "Unexpected binary survivors.")
        (audit/"h3_source_inputs.gp").write_text("ids="+str([t["id"] for t in survivors])+";\nHS=["+",".join(gp_matrix(t["H"]) for t in survivors)+"];\n")
        gpr("h3_source_orbits.gp")
        pyg("verify_source_orbits.py")
        reps = rows(audit/"independent_h3_representatives.jsonl")
        rank = {r["id"]:r.get("short_rank") for r in b3}
        high = [t for t in reps if rank[t["id"]]==4]
        low = [t for t in reps if rank[t["id"]]!=4]
        require(len(reps)==363 and len(high)==362 and [t["id"] for t in low]==[11622], "Unexpected final representative partition.")
        highflat = audit/"independent_h3_representatives_highspan.flat"
        flat(high, highflat)
        run("simultaneous_baseline", [audit/"independent_simultaneous", highflat], stdout=audit/"independent_h3_simultaneous_baseline.jsonl")
        pyg("verify_simultaneous.py", "--targets", audit/"independent_h3_targets.jsonl",
            "--results", audit/"independent_h3_simultaneous_baseline.jsonl",
            "--quadruples", str(highflat)+".baseline_quads",
            "--output", audit/"independent_h3_simultaneous_verification.json")
        pyg("test_affine_enum.py")
        pyg("prepare_lowspan_target.py")
        pyg("affine_embedding_audit.py", audit/"h3_lowspan_rebased.json", "lowspan_affine")
        # Controls contain a known actual M sublattice and the abstract N.
        gpr("embedding_controls.gp")
        run("embedding_controls", [audit/"independent_embedding_sieve", audit/"embedding_controls.flat"], stdout=audit/"embedding_controls_result.jsonl")
        controls = rows(audit/"embedding_controls_result.jsonl")
        require([r["status"] for r in controls]==["survives","binary_absent","binary_absent","det_lt16"], "Embedding controls failed.")
        t64flat = audit/"T64_control.flat"
        t64flat.write_text((audit/"embedding_controls.flat").read_text().splitlines()[0]+"\n")
        run("T64_simultaneous_positive", [audit/"independent_simultaneous", t64flat], stdout=audit/"T64_simultaneous.jsonl")
        require(rows(audit/"T64_simultaneous.jsonl")[0]["integral"]>0, "Known integral embedding was incorrectly rejected.")
        if args.check_d3_inputs:
            run("d3_existing_balls_exact_recheck", [exact, "--batch", root/"output/d3/benchmark/all_exact_inputs.txt", "55"], stdout=audit/"d3_dependency_recheck.jsonl")
            dep = rows(audit/"d3_dependency_recheck.jsonl")
            old = rows(root/"output/d3/benchmark/all_exact_results.jsonl")
            require(len(dep)==len(old)==30423 and all(r["complete"] for r in dep), "Incomplete d3 dependency sphere check.")
            require(all((a["dimension"],a["pairs"])==(b["dimension"],b["pairs"]) for a,b in zip(dep,old)), "d3 exact-ball counts differ from the existing certificate.")
            require(all(r["pairs"]==0 for r in dep if r["dimension"]==18), "A forbidden rank18 quotient vector was found.")
            summary["dependency"].update(rechecked=True,forms=30423,dimensions=dict(collections.Counter(r["dimension"] for r in dep)))
        # The coverage decision uses regenerated data, never saved conclusions.
        h2result = json.loads((audit/"h2_candidate_witnesses_summary.json").read_text())
        sim = rows(audit/"independent_h3_simultaneous_baseline.jsonl")
        inversecheck = json.loads((audit/"independent_h3_simultaneous_verification.json").read_text())
        affine = json.loads((audit/"lowspan_affine/summary.json").read_text())
        smap = json.loads((audit/"h3_source_orbit_map.json").read_text())
        require(h2result["all_excluded"] and h2result["candidates"]==23, "Incomplete h2 coverage.")
        require({r["id"] for r in b3}==set(byid) and len(b3)==len(byid), "Incomplete binary coverage.")
        require({r[0] for r in smap}=={t["id"] for t in survivors} and len(smap)==len(survivors), "Incomplete source-orbit coverage.")
        require({r["id"] for r in sim}=={t["id"] for t in high} and len(sim)==len(high), "Incomplete simultaneous coverage.")
        require(all(r["integral"]==0 for r in sim) and inversecheck["complete"] and inversecheck["integral_embeddings"]==0, "Simultaneous integral exclusions failed.")
        require(affine["completed"] and affine["embeddings_found"]==0 and [r["output_states"] for r in affine["stages"]]==[49,9,0], "Lowspan exclusion failed.")
        excluded = {r["id"] for r in b3 if r["status"]!="survives"}
        eliminated_reps = {r["id"] for r in sim}|{11622}
        excluded |= {t for t,rep,_ in smap if rep in eliminated_reps}
        require(excluded==set(byid), "Some h3 target lacks a completed exclusion.")
        save()
        pyg("check_final_coverage.py")
        require(json.loads((audit/"FINAL_AUDIT_SUMMARY.json").read_text())["verified"], "Final coverage aggregator did not certify the result.")
        summary.update(status="complete",rank4_h_le3_excluded=True,h1="23h>=24 excludes h=1",
                       h2_candidates=23,h3_candidates=13272,h3_all_targets_excluded=len(excluded),
                       source_representatives=363,simultaneous_representatives=362,
                       lowspan_affine_output_counts=[49,9,0],wall_seconds=time.perf_counter()-start)
        summary["output_hashes"] = {name:sha(audit/name) for name in
            ["independent_h2_targets.jsonl","independent_h3_targets.jsonl","independent_h3_binary_sieve.jsonl",
             "h3_source_orbit_map.json","independent_h3_simultaneous_baseline.jsonl","lowspan_affine/summary.json"]}
        save()
        print(f"AUDIT_COMPLETE {summary_path}", flush=True)
    except BaseException as exc:
        summary.update(status="failed",error=str(exc),wall_seconds=time.perf_counter()-start)
        save()
        print(f"AUDIT_FAILED {summary_path}: {exc}", file=sys.stderr, flush=True)
        raise

if __name__ == "__main__":
    main()
