#!/bin/zsh
# Le lanceur de la nuit du 2026-10-17. Il relit queue.txt avant chaque tache,
# prend la premiere dont l'identifiant n'est pas dans done.txt, verifie
# irrisor et la place disque, la lance, et consigne le resultat. Une tache qui
# echoue n'arrete pas la nuit. Rien n'est jamais efface.
#
# queue.txt, une tache par ligne:  ID|GO_IRRISOR|GO_LOCAL|commande
# (les Go libres qu'il faut sur irrisor et sur le disque interne, marge comprise)

REPO=/Users/eartigau/GitHubProjects/pca2d_preclean
HERE=$REPO/outputs/reports/nuit17
PY=/opt/anaconda3/envs/pca2d-preclean/bin/python
QUEUE=$HERE/queue.txt
DONE=$HERE/done.txt
LOG=$HERE/runner.log
IRRISOR=/Volumes/irrisor/pca2d_preclean/science

ts() { $PY -c "import datetime as d; n=d.datetime.now(); print(n.strftime('%y%m%d %H:%M:%S.')+'%02d' % (n.microsecond//10000))" }
say() {   # say niveau message: info (vert), value (bleu), warn (orange), stop (rouge)
  local c; case $1 in info) c='\e[32m';; value) c='\e[34m';; warn) c='\e[38;5;208m';; stop) c='\e[31m';; esac
  print -r -- "$(ts) | $2" >> $LOG
  print -- "${c}$(ts) | $2\e[0m"
}
free_gb() { df -g "$1" | awk 'NR==2{print $4}' }

cd $REPO
# UNE VERSION POUR TOUTE LA NUIT: refus de demarrer sur du code non commite
dirty=$(git status --porcelain --untracked-files=all pca2d tools config.yaml variants tests)
if [[ -n "$dirty" ]]; then
  say stop "uncommitted code, the night does not start:"; print -r -- "$dirty"
  exit 1
fi
ldirty=$(git -C /Users/eartigau/lbl-rapide status --porcelain lbl)
if [[ -n "$ldirty" ]]; then
  say stop "uncommitted LBL code, the night does not start"; exit 1
fi
print -r -- "pca2d $(git rev-parse HEAD)  lbl-rapide $(git -C /Users/eartigau/lbl-rapide rev-parse HEAD)  started $(date '+%Y-%m-%d %H:%M:%S')" >> $HERE/VERSIONS.txt
say value "frozen: pca2d $(git rev-parse --short HEAD), lbl-rapide $(git -C /Users/eartigau/lbl-rapide rev-parse --short HEAD)"
touch $DONE
caffeinate -dimsu -w $$ &
say info "runner started, pid $$, queue $QUEUE"

while true; do
  next=""
  while IFS= read -r line; do
    [[ -z "$line" || "$line" == \#* ]] && continue
    id=${line%%|*}
    grep -q "^$id|" $DONE && continue
    next=$line; break
  done < $QUEUE
  if [[ -z "$next" ]]; then
    say info "queue empty: NUIT TERMINEE"
    break
  fi
  id=${next%%|*}; rest=${next#*|}
  need_i=${rest%%|*}; rest=${rest#*|}
  need_l=${rest%%|*}; cmd=${rest#*|}

  # irrisor doit etre la; on attend jusqu'a 30 minutes, puis on arrete
  waited=0
  while [[ ! -d $IRRISOR ]]; do
    [[ $waited -eq 0 ]] && say warn "$id: irrisor absent, waiting"
    sleep 30; waited=$((waited + 30))
    if [[ $waited -ge 1800 ]]; then
      say stop "$id: irrisor absent for 30 minutes, the night stops here"
      exit 1
    fi
  done

  fi_=$(free_gb /Volumes/irrisor); fl=$(free_gb /)
  if [[ $fi_ -lt $need_i || $fl -lt $need_l ]]; then
    say warn "$id: SKIPPED, not enough disk (irrisor $fi_ Go free for $need_i, internal $fl for $need_l)"
    print -r -- "$id|skipped-disk|$(ts)|$(ts)|-|$(git rev-parse --short HEAD)" >> $DONE
    continue
  fi

  say info "$id: START  (irrisor $fi_ Go free, internal $fl)  $cmd"
  start=$(ts)
  zsh -c "$cmd" > $HERE/$id.log 2>&1
  code=$?
  end=$(ts)
  if [[ $code -eq 0 ]]; then
    say info "$id: DONE   code 0"
  else
    say warn "$id: FAILED code $code, see $HERE/$id.log"
  fi
  print -r -- "$id|code=$code|$start|$end|$cmd|$(git rev-parse --short HEAD)" >> $DONE
done
