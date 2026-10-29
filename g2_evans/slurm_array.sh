#!/bin/bash
#SBATCH --job-name=g2_evans
#SBATCH --array=0-49
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=24:00:00
#SBATCH --output=logs/g2_%A_%a.out
# Varredura G2 de K' para (1;1;0,5).  Cada tarefa processa as linhas task, task+NTASKS, ... da lista $BOXES
# (tasks_A.txt, tasks_B.txt, tasks_C.txt ou todo.txt; ver README), com subdivisao adaptativa.
# Certificados em certs/<lista>_task_XXXX.jsonl.  RETOMAVEL: resubmeta o mesmo array apos queda/limite de tempo.
#   sbatch --array=0-49 --export=ALL,BOXES=tasks_A.txt,NTASKS=50 slurm_array.sh
BOXES=${BOXES:-tasks_A.txt}
NTASKS=${NTASKS:-50}
cd "$SLURM_SUBMIT_DIR"
mkdir -p certs logs
python3 sweep.py run --modo "${MODO:-contorno}" --boxes "$BOXES" --lamstar lamstar.txt --task "$SLURM_ARRAY_TASK_ID" --ntasks "$NTASKS" \
        --out certs/$(basename "$BOXES" .txt)_task_$(printf %04d "$SLURM_ARRAY_TASK_ID").jsonl --timeout 21600 --maxdepth 4
