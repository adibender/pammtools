# One-time setup on LRZ CoolMUC-4

Run these from your **local** machine with the multiplexed `ssh lrz` connection
active (open `ssh lrz` in a separate terminal first if needed).

```bash
# 1. clone the repo and check out the study branch (adjust remote URL if needed)
ssh lrz "git clone https://github.com/adibender/pammtools.git ~/pammtools 2>/dev/null; \
         cd ~/pammtools && git fetch origin claude/confident-cori-dng5at && \
         git checkout claude/confident-cori-dng5at && git pull"

# 2. install pammtools + all dependencies into the user library (login node is
#    fine for this; pak resolves CRAN deps and installs the local package)
ssh lrz 'module load slurm_setup && module load r/4.3.3-gcc13-mkl && \
  Rscript -e "if (!requireNamespace(\"pak\", quietly = TRUE)) install.packages(\"pak\", repos = \"https://cloud.r-project.org\")" && \
  Rscript -e "pak::pkg_install(\"local::~/pammtools\", lib = \"~/R/x86_64-pc-linux-gnu-library/4.3\", ask = FALSE)"'

# 3. quick smoke test (small, login node OK)
ssh lrz 'module load slurm_setup && module load r/4.3.3-gcc13-mkl && \
  Rscript -e "library(pammtools); packageVersion(\"pammtools\")"'

# 4. create the logs dir and submit the full run
ssh lrz "mkdir -p ~/pammtools/attic/simulations/ci-surv-prob/logs ~/pammtools/attic/simulations/ci-surv-prob/results/raw"
ssh lrz "cd ~/pammtools/attic/simulations/ci-surv-prob && sbatch lrz/job-full.slurm"
# note the job id, then monitor (sparingly; minutes-scale intervals):
#   ssh lrz "squeue --clusters=serial -u \$USER"
#   ssh lrz "sacct --clusters=serial -j <JOBID> --format=JobID,State,ExitCode,Elapsed,MaxRSS"
#   ssh lrz "cat ~/pammtools/attic/simulations/ci-surv-prob/logs/pamm-ci-surv.<JOBID>.out"

# 5. fetch results when done
scp -r lrz:~/pammtools/attic/simulations/ci-surv-prob/output ./output-full
# (or commit output/ on LRZ and push, since the clone is a git repo)
```

Notes

- `run-sim.R` checkpoints per scenario (`results/raw/scen-XX.rds`) and only
  computes missing replications, so resubmitting the same job after a timeout
  or failure resumes the run. RNG streams are tied to (scenario, rep), not to
  scheduling, so resumed runs reproduce exactly what a single run would give.
- To split the work across several jobs instead, submit multiple copies with
  disjoint `--scenarios` ranges, e.g.
  `sed 's/--reps 500/--reps 500 --scenarios 1-12/' lrz/job-full.slurm | sbatch`
  (each scenario file is written by exactly one job, so ranges must not
  overlap).
- For more than 500 replications in contested cells, rerun with e.g.
  `--reps 1000 --scenarios <ids>`; streams support up to 1000 reps per
  scenario without changing existing results.
