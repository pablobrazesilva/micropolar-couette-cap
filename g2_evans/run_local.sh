#!/bin/bash
# Sem SLURM: P processos em paralelo na maquina local, uma faixa por vez (ex.: ./run_local.sh 16 tasks_A.txt)
P=${1:-$(nproc)}; BOXES=${2:-tasks_A.txt}
mkdir -p certs logs
for i in $(seq 0 $((P-1))); do
  python3 sweep.py run --modo "${MODO:-contorno}" --boxes "$BOXES" --lamstar lamstar.txt --task $i --ntasks $P \
     --out certs/$(basename "$BOXES" .txt)_task_$(printf %04d $i).jsonl > logs/local_$(basename "$BOXES" .txt)_$i.log 2>&1 &
done
wait
python3 sweep.py verify --boxes boxes.txt --certs 'certs/*.jsonl'
