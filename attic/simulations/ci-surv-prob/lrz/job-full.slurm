#!/bin/bash
#SBATCH -J pamm-ci-surv
#SBATCH -o ./logs/%x.%j.out
#SBATCH -e ./logs/%x.%j.err
#SBATCH -D /dss/dsshome1/lxc04/ri68suz2/pammtools/attic/simulations/ci-surv-prob
#SBATCH --get-user-env
#SBATCH --clusters=serial
#SBATCH --partition=serial_std
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16000mb
#SBATCH --time=08:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=fabian.scheipl@stat.uni-muenchen.de
#SBATCH --export=NONE

# Full run of the CI coverage study (pammtools issue #285):
# all scenarios x 500 replications. Estimated ~1-2 h on 16 cores; results are
# checkpointed per scenario, so a resubmission resumes where it stopped.

module load slurm_setup
module load r/4.3.3-gcc13-mkl

# --export=NONE strips env vars; ensure user R library is discoverable
export R_LIBS_USER="${HOME}/R/x86_64-pc-linux-gnu-library/4.3"

Rscript run-sim.R \
  --reps 500 \
  --cores "${SLURM_CPUS_PER_TASK}" \
  --seed 285 \
  --out results/raw

Rscript eval-sim.R --raw results/raw --out output --tag -full
