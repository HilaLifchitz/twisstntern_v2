#!/usr/bin/env python3
"""
COMPREHENSIVE CLI TEST SUITE FOR TWISSTNTERN & TWISSTNTERN_SIMULATE

Covers (at least once) every flag shown in --help for both CLIs:
- twisstntern:
  INPUT/--input, -o/--output, --verbose, --granularity, --downsample, --axis, --colormap,
  tree inputs with --taxon-names, --outgroup, --topology-mapping
- twisstntern_simulate:
  --help, --get-config, -c/--config, -o/--output, --verbose, --quiet, --log-file, --seed,
  --granularity, --colormap, --density-colormap, --override, --downsample, --downsampleKB,
  --topology-mapping

Usage:
  conda activate twisstntern_env
  cd /home/hlifchit/projects/twissting_baby
  python comprehensive_cli_test.py

Outputs:
  - All run outputs: CLI_TEST_RUNS_YYYYMMDD_HHMMSS/
  - JSON report: CLI_TEST_RUNS_.../cli_test_report.json
  - Non-zero exit code if any critical (expected-success) test fails
"""

import sys
import time
import json
import subprocess
from pathlib import Path

try:
    import yaml
except Exception:
    yaml = None  # simulate phases that require writing configs will be skipped if yaml missing


class CLITestSuite:
    def __init__(self):
        self.project_root = Path(__file__).resolve().parent
        self.py = sys.executable  # use active environment (twisstntern_env)
        ts = time.strftime("%Y%m%d_%H%M%S")
        self.run_root = self.project_root / f"CLI_TEST_RUNS_{ts}"
        self.run_root.mkdir(parents=True, exist_ok=True)
        self.results = []
        self.t0 = time.time()

        self.data_dir = self.project_root / "Examples" / "data_files"
        self.csv_variants = [
            ("A_m0.csv", "numeric headers CSV"),
            ("migration_topology_weights.csv", "T1/T2/T3 headers CSV"),
            ("neurospora_weights.csv", "topo1/topo2/topo3 headers CSV"),
        ]

    def outdir(self, name: str) -> Path:
        p = self.run_root / name
        p.mkdir(parents=True, exist_ok=True)
        return p

    def run(self, argv, desc, expect_success=True, timeout=240, cwd=None):
        if cwd is None:
            cwd = str(self.project_root)
        if argv and argv[0] == "-m":
            cmd = [self.py] + argv
        elif argv and argv[0] == self.py:
            cmd = argv
        else:
            cmd = [self.py] + argv

        print(f"\n🔍 {desc}")
        print(f"   Command: {' '.join(cmd)}")
        try:
            proc = subprocess.run(
                cmd, cwd=cwd, capture_output=True, text=True, timeout=timeout
            )
            ok = (proc.returncode == 0) if expect_success else (proc.returncode != 0)
            status = "✅ PASS" if ok else "❌ FAIL"
            print(f"   {status} (rc={proc.returncode})")
            if proc.stdout.strip():
                print(f"   📤 stdout: {proc.stdout.strip()[:240]}{'...' if len(proc.stdout) > 240 else ''}")
            if proc.stderr.strip():
                print(f"   🚨 stderr: {proc.stderr.strip()[:240]}{'...' if len(proc.stderr) > 240 else ''}")
            self.results.append(
                {
                    "desc": desc,
                    "cmd": cmd,
                    "ok": ok,
                    "expected_success": expect_success,
                    "rc": proc.returncode,
                    "stdout": proc.stdout,
                    "stderr": proc.stderr,
                    "t": round(time.time() - self.t0, 1),
                }
            )
            return ok, proc
        except subprocess.TimeoutExpired:
            print(f"   ⏰ TIMEOUT after {timeout}s")
            self.results.append(
                {
                    "desc": desc,
                    "cmd": cmd,
                    "ok": False,
                    "expected_success": expect_success,
                    "rc": -1,
                    "stdout": "",
                    "stderr": f"Timeout after {timeout}s",
                    "t": round(time.time() - self.t0, 1),
                }
            )
            return False, None
        except Exception as e:
            print(f"   💥 EXCEPTION: {e}")
            self.results.append(
                {
                    "desc": desc,
                    "cmd": cmd,
                    "ok": False,
                    "expected_success": expect_success,
                    "rc": -2,
                    "stdout": "",
                    "stderr": str(e),
                    "t": round(time.time() - self.t0, 1),
                }
            )
            return False, None

    # PHASE 1: Basics / help
    def phase_help(self):
        print("\n" + "=" * 80)
        print("📖 PHASE 1: HELP / USAGE")
        print("=" * 80)
        tests = [
            (["-m", "twisstntern", "--help"], "twisstntern --help"),
            (["-m", "twisstntern", "-h"], "twisstntern -h"),
            (["-m", "twisstntern_simulate", "--help"], "twisstntern_simulate --help"),
            (["-m", "twisstntern_simulate", "-h"], "twisstntern_simulate -h"),
            (["-m", "twisstntern_simulate", "--get-config"], "twisstntern_simulate --get-config"),
        ]
        for argv, desc in tests:
            self.run(argv, desc, expect_success=True, timeout=40)

    # PHASE 2: Invalid args (expected failures)
    def phase_invalid(self):
        print("\n" + "=" * 80)
        print("🚨 PHASE 2: INVALID ARGUMENTS (EXPECTED FAIL)")
        print("=" * 80)
        tests = [
            (["-m", "twisstntern"], "twisstntern (no args) — expect fail"),
            (["-m", "twisstntern_simulate"], "twisstntern_simulate (no args) — expect fail"),
            (["-m", "twisstntern", "DOES_NOT_EXIST.csv"], "twisstntern nonexistent input — expect fail"),
            (["-m", "twisstntern_simulate", "-c", "DOES_NOT_EXIST.yaml"], "simulate nonexistent config — expect fail"),
            (["-m", "twisstntern", str(self.data_dir / "A_m0.csv"), "--downsample", "invalid"], "invalid downsample format — expect fail"),
            (["-m", "twisstntern", str(self.data_dir / "A_m0.csv"), "--downsample", "5+10"], "invalid downsample constraint i>=N — expect fail"),
            (["-m", "twisstntern", str(self.data_dir / "A_m0.csv"), "--granularity", "nope"], "invalid granularity value — expect fail"),
        ]
        for argv, desc in tests:
            self.run(argv, desc, expect_success=False, timeout=40)

    # PHASE 3: CSV analysis basics (different header variants)
    def phase_csv_basic(self):
        print("\n" + "=" * 80)
        print("📊 PHASE 3: CSV ANALYSIS (BASIC)")
        print("=" * 80)
        for fname, label in self.csv_variants:
            out = self.outdir(f"csv_basic_{Path(fname).stem}")
            argv = [
                "-m", "twisstntern",
                str(self.data_dir / fname),
                "-o", str(out),
                "--verbose",
            ]
            self.run(argv, f"CSV basic: {label}", expect_success=True, timeout=120)

    # PHASE 4: Advanced twisstntern flags
    def phase_advanced(self):
        print("\n" + "=" * 80)
        print("🔬 PHASE 4: ADVANCED TWISSTNTERN OPTIONS")
        print("=" * 80)
        test_file = str(self.data_dir / "A_m0.csv")
        cases = [
            (["--granularity", "superfine"], "granularity superfine"),
            (["--granularity", "coarse"], "granularity coarse"),
            (["--granularity", "0.05"], "granularity 0.05"),
            (["--downsample", "10"], "downsample 10"),
            (["--downsample", "10+3"], "downsample 10+3"),
            (["--axis", "T2", "T1", "T3"], "axis order T2 T1 T3"),
            (["--colormap", "viridis"], "colormap viridis"),
            (["--colormap", "plasma"], "colormap plasma"),
            (["--colormap", "inferno"], "colormap inferno"),
        ]
        for extra, label in cases:
            out = self.outdir(f"advanced_{label.replace(' ', '_')}")
            argv = ["-m", "twisstntern", test_file, "-o", str(out)] + extra
            self.run(argv, f"twisstntern advanced: {label}", expect_success=True, timeout=120)

    # PHASE 5: Tree inputs & flags (twisstntern)
    def phase_tree_inputs(self):
        print("\n" + "=" * 80)
        print("🌳 PHASE 5: TREE INPUTS & TOPOLOGY FLAGS (TWISSTNTERN)")
        print("=" * 80)
        out = self.outdir("tree_inputs")
        newick_file = out / "tiny.newick"
        newick_file.write_text("(O,(P1,(P2,P3)));")

        # A) positional + taxa/outgroup
        self.run(
            [
                "-m", "twisstntern",
                str(newick_file),
                "--taxon-names", "O", "P1", "P2", "P3",
                "--outgroup", "O",
                "-o", str(out / "positional"),
                "--verbose",
            ],
            "newick positional + --taxon-names/--outgroup",
            expect_success=True,
            timeout=180,
        )

        # B) --input alias
        self.run(
            [
                "-m", "twisstntern",
                "--input", str(newick_file),
                "--taxon-names", "O", "P1", "P2", "P3",
                "--outgroup", "O",
                "-o", str(out / "input_flag"),
            ],
            "newick via --input alias + taxa/outgroup",
            expect_success=True,
            timeout=180,
        )

        # C) topology mapping
        self.run(
            [
                "-m", "twisstntern",
                str(newick_file),
                "--taxon-names", "O", "P1", "P2", "P3",
                "--outgroup", "O",
                "--topology-mapping",
                'T1="(O,(P3,(P1,P2)))"; T2="(O,(P1,(P2,P3)))"; T3="(O,(P2,(P1,P3)))";',
                "-o", str(out / "topology_mapping"),
            ],
            "newick + --topology-mapping",
            expect_success=True,
            timeout=180,
        )

    # Helper: write a minimal simulate config (if yaml available)
    def write_minimal_sim_config(self, path: Path):
        if yaml is None:
            return False
        cfg = {
            "population_labels": {"p1": "A", "p2": "B", "p3": "C", "p23": "AB_anc", "p123": "ABC_anc", "O": "Outgroup", "ANC": "Root"},
            "splits": [
                {"time": 50, "derived_pop1": "p2", "derived_pop2": "p3", "ancestral_pop": "p23"},
                {"time": 100, "derived_pop1": "p23", "derived_pop2": "p1", "ancestral_pop": "p123"},
                {"time": 150, "derived_pop1": "p123", "derived_pop2": "O", "ancestral_pop": "ANC"},
            ],
            "populations": [
                {"name": "O", "Ne": 1000, "sample_size": 5},
                {"name": "p1", "Ne": 1000, "sample_size": 5},
                {"name": "p2", "Ne": 1000, "sample_size": 5},
                {"name": "p3", "Ne": 1000, "sample_size": 5},
                {"name": "p23", "Ne": 1000},
                {"name": "p123", "Ne": 1000},
                {"name": "ANC", "Ne": 1000},
            ],
            "migration": {f"{a}>{b}": 0.0 for a in ["p1", "p2", "p3", "O"] for b in ["p1", "p2", "p3", "O"] if a != b},
            "simulation_mode": "locus",
            "ploidy": 1,
            "seed": 12345,
            "n_loci": 20,
            "locus_length": 1,
            "rec_rate": 1e-8,
            "chromosome_length": 1_000_000,
            "mutation_rate": 0,
        }
        path.write_text(yaml.safe_dump(cfg))
        return True

    # PHASE 6: Simulate basics + overrides
    def phase_simulate(self):
        print("\n" + "=" * 80)
        print("🧬 PHASE 6: SIMULATION (BASICS & OVERRIDES)")
        print("=" * 80)

        # --get-config to file
        out = self.outdir("simulate_get_config")
        cfg_dl = out / "downloaded_config.yaml"
        self.run(["-m", "twisstntern_simulate", "--get-config", str(cfg_dl)], "simulate --get-config → file", expect_success=True, timeout=60)

        # pick config: project template or minimal
        project_cfg = self.project_root / "config_template.yaml"
        if project_cfg.exists():
            use_cfg = project_cfg
        else:
            synth_cfg = out / "minimal_config.yaml"
            if not self.write_minimal_sim_config(synth_cfg):
                print("   ⚠️  YAML missing, skipping simulate phases")
                return
            use_cfg = synth_cfg

        # basic simulate (verbose)
        sim_out = self.outdir("simulate_basic")
        self.run(
            ["-m", "twisstntern_simulate", "-c", str(use_cfg), "-o", str(sim_out), "--verbose"],
            "simulate basic (verbose)",
            expect_success=True,
            timeout=360,
        )

        # simulate with overrides + colormap + granularity
        sim_out2 = self.outdir("simulate_overrides")
        self.run(
            [
                "-m", "twisstntern_simulate",
                "-c", str(use_cfg),
                "-o", str(sim_out2),
                "--override", "n_loci=10",
                "--override", "populations.p1.Ne=2000",
                "--granularity", "0.1",
                "--colormap", "plasma",
            ],
            "simulate with overrides, granularity, colormap",
            expect_success=True,
            timeout=360,
        )

    # PHASE 7: Simulate extra flags (quiet/log/seed/density/downsample/downsampleKB/topology-mapping)
    def phase_simulate_extras(self):
        print("\n" + "=" * 80)
        print("🧪 PHASE 7: SIMULATE EXTRA FLAGS")
        print("=" * 80)
        out = self.outdir("simulate_extra")
        project_cfg = self.project_root / "config_template.yaml"

        # base config for locus mode
        if project_cfg.exists():
            base_cfg = project_cfg
        else:
            base_cfg = out / "minimal_config.yaml"
            if not self.write_minimal_sim_config(base_cfg):
                print("   ⚠️  YAML missing, skipping simulate extra flags")
                return

        # --quiet + --log-file + --seed
        self.run(
            [
                "-m", "twisstntern_simulate",
                "-c", str(base_cfg),
                "-o", str(out / "quiet"),
                "--quiet",
                "--log-file", str(out / "quiet.log"),
                "--seed", "4242",
            ],
            "simulate --quiet --log-file --seed",
            expect_success=True,
            timeout=360,
        )

        # --density-colormap
        self.run(
            [
                "-m", "twisstntern_simulate",
                "-c", str(base_cfg),
                "-o", str(out / "density_cmap"),
                "--density-colormap", "coolwarm",
            ],
            "simulate --density-colormap",
            expect_success=True,
            timeout=360,
        )

        # --downsample (locus mode)
        self.run(
            [
                "-m", "twisstntern_simulate",
                "-c", str(base_cfg),
                "-o", str(out / "downsample_locus"),
                "--downsample", "5",
            ],
            "simulate --downsample (locus mode)",
            expect_success=True,
            timeout=360,
        )

        # chromosome mode + --downsampleKB
        if yaml is not None:
            try:
                base = yaml.safe_load(Path(base_cfg).read_text())
                base["simulation_mode"] = "chromosome"
                base["chromosome_length"] = 2_000_000
                base["rec_rate"] = 1e-8
                base["mutation_rate"] = 0
                chrom_cfg = out / "chromosome_config.yaml"
                chrom_cfg.write_text(yaml.safe_dump(base))
                self.run(
                    [
                        "-m", "twisstntern_simulate",
                        "-c", str(chrom_cfg),
                        "-o", str(out / "downsample_kb"),
                        "--downsampleKB", "100kb+50kb",
                    ],
                    "simulate chromosome --downsampleKB 100kb+50kb",
                    expect_success=True,
                    timeout=420,
                )
            except Exception as e:
                print(f"   ⚠️  Skipping chromosome downsampleKB test: {e}")

        # --topology-mapping string
        self.run(
            [
                "-m", "twisstntern_simulate",
                "-c", str(base_cfg),
                "-o", str(out / "topology_mapping"),
                "--topology-mapping",
                'T1="(0,(3,(1,2)))"; T2="(0,(1,(2,3)))"; T3="(0,(2,(1,3)))";',
            ],
            "simulate --topology-mapping",
            expect_success=True,
            timeout=360,
        )

    # PHASE 8: Edge cases (csv)
    def phase_edge(self):
        print("\n" + "=" * 80)
        print("⚠️  PHASE 8: EDGE CASES (CSV)")
        print("=" * 80)
        out = self.outdir("edge_cases")

        # minimal valid CSV
        minimal = out / "minimal.csv"
        minimal.write_text("T1,T2,T3\n1,0,0\n0,1,0\n0,0,1\n")
        self.run(["-m", "twisstntern", str(minimal), "-o", str(out / "minimal_out")], "minimal valid CSV", expect_success=True, timeout=60)

        # empty CSV — expect fail
        empty = out / "empty.csv"
        empty.write_text("")
        self.run(["-m", "twisstntern", str(empty), "-o", str(out / "empty_out")], "empty CSV (expect fail)", expect_success=False, timeout=40)

        # headers-only — expect fail
        headers_only = out / "headers_only.csv"
        headers_only.write_text("T1,T2,T3\n")
        self.run(["-m", "twisstntern", str(headers_only), "-o", str(out / "headers_only_out")], "headers-only CSV (expect fail)", expect_success=False, timeout=40)

        # invalid granularity (constraint)
        self.run(
            ["-m", "twisstntern", str(self.data_dir / "A_m0.csv"), "--granularity", "0.3", "-o", str(out / "bad_granularity")],
            "invalid granularity 0.3 (expect fail)",
            expect_success=False,
            timeout=60,
        )

    # REPORT
    def report(self):
        print("\n" + "=" * 80)
        print("📊 COMPREHENSIVE REPORT")
        print("=" * 80)
        total = len(self.results)
        passed = sum(1 for r in self.results if r["ok"])
        failed = total - passed
        print(f"\nTotal tests: {total}")
        print(f"✅ Passed: {passed}")
        print(f"❌ Failed: {failed}")
        print(f"⏱️  Duration: {round(time.time() - self.t0, 1)}s")

        critical = [r for r in self.results if not r["ok"] and r["expected_success"]]
        if critical:
            print(f"\n🚨 Critical failures ({len(critical)}):")
            for r in critical:
                print(f" - {r['desc']}")
                if r["stderr"]:
                    print(f"   Error: {r['stderr'].strip()[:200]}{'...' if len(r['stderr']) > 200 else ''}")

        report_path = self.run_root / "cli_test_report.json"
        report_path.write_text(json.dumps(self.results, indent=2))
        print(f"\n📝 JSON report: {report_path}")
        print(f"📁 Outputs dir: {self.run_root}")
        sys.exit(1 if critical else 0)


def main():
    print("🧪 Comprehensive CLI Test Suite (twisstntern + twisstntern_simulate)")
    print(f"📁 Project: {Path(__file__).resolve().parent}")
    print(f"🐍 Python: {sys.executable}")
    suite = CLITestSuite()
    suite.phase_help()
    suite.phase_invalid()
    suite.phase_csv_basic()
    suite.phase_advanced()
    suite.phase_tree_inputs()
    suite.phase_simulate()
    suite.phase_simulate_extras()
    suite.phase_edge()
    suite.report()


if __name__ == "__main__":
    main()
