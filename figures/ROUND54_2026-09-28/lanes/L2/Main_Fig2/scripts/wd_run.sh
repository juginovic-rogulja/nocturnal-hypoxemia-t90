#!/bin/bash
# Locations come from the repository's paths.py: run   eval "$(python3 paths.py)"   from the repository root first (it defines the T90_* variables).
# wd_run.sh NAME MAX_S RSS_KB -- CMD ARGS...
# Round 38 L5 SHEETS watchdog (2026-09-16, after the 12:26 kernel panic). Runs ONE command and polls the resident set of the
# process AND its descendants every 10 s (ps -o rss=). Kills the whole tree with SIGKILL when the tree's summed RSS exceeds
# RSS_KB or the wall time exceeds MAX_S. Logs: $WD_LOGDIR/NAME.log (one poll line per 10 s), NAME.out (stdout+stderr of the
# command), NAME.status (rc wall_s peak_rss_kb peak_single_kb reason), and one line appended to $WD_LOGDIR/RSS_LOG.tsv.
# Refuses to start while another wd_run.sh command is still running (one sheet, one build, one verifier at a time).
NAME=$1; LIMIT=$2; RSS_LIMIT=$3; shift 3
[ "$1" = "--" ] && shift
W=${WD_LOGDIR:-${T90_FIGURE_ROOT}/ROUND38_2026-09-16/figures/logs/watchdog}
mkdir -p "$W"
LOG=$W/$NAME.log; ST=$W/$NAME.status; OUT=$W/$NAME.out; LOCK=$W/.running
rm -f "$ST"
if [ -f "$LOCK" ] && kill -0 "$(cat "$LOCK" 2>/dev/null)" 2>/dev/null; then
  echo "rc=-1 wall_s=0 peak_rss_kb=0 peak_single_kb=0 reason=refused_another_wd_run_active($(cat "$LOCK"))" > "$ST"; cat "$ST"; exit 2
fi
tree_pids() { local p=$1; echo "$p"; local c; for c in $(pgrep -P "$p" 2>/dev/null); do tree_pids "$c"; done; }
start=$(date +%s)
"$@" > "$OUT" 2>&1 &
pid=$!
echo $pid > "$LOCK"
echo "$(date +%T) START pid=$pid limit=${LIMIT}s rss_limit=${RSS_LIMIT}KB cmd: $*" >> "$LOG"
peak=0; peak1=0; reason=finished; npoll=0
while kill -0 "$pid" 2>/dev/null; do
  npoll=$((npoll + 1)); if [ "$npoll" -le 20 ]; then sleep 1; else sleep 3; fi   # coordinator 18:14: 3-second polls
  pids=$(tree_pids "$pid" | tr '\n' ' ')
  total=0; single=0
  for p in $pids; do
    r=$(ps -o rss= -p "$p" 2>/dev/null | tr -d ' '); [ -z "$r" ] && continue
    total=$((total + r)); [ "$r" -gt "$single" ] && single=$r
  done
  now=$(date +%s); el=$((now - start))
  [ "$total" -gt "$peak" ] && peak=$total
  [ "$single" -gt "$peak1" ] && peak1=$single
  echo "$(date +%T) t=${el}s rss_tree=${total}KB rss_max_proc=${single}KB peak=${peak}KB pids=$(echo $pids | wc -w | tr -d ' ')" >> "$LOG"
  if [ "$total" -gt "$RSS_LIMIT" ]; then reason=killed_rss; for p in $pids; do kill -9 "$p" 2>/dev/null; done; break; fi
  if [ "$el" -gt "$LIMIT" ]; then reason=killed_timeout; for p in $pids; do kill -9 "$p" 2>/dev/null; done; break; fi
done
wait "$pid" 2>/dev/null; rc=$?
end=$(date +%s); wall=$((end - start))
rm -f "$LOCK"
echo "$(date +%T) END rc=$rc wall=${wall}s peak_rss_tree=${peak}KB peak_single=${peak1}KB reason=$reason" >> "$LOG"
echo "rc=$rc wall_s=$wall peak_rss_kb=$peak peak_single_kb=$peak1 reason=$reason" > "$ST"
[ -f "$W/RSS_LOG.tsv" ] || printf "name\tstart\trc\twall_s\tpeak_rss_tree_kb\tpeak_single_kb\treason\tcmd\n" > "$W/RSS_LOG.tsv"
printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" "$NAME" "$(date -r $start +%FT%T)" "$rc" "$wall" "$peak" "$peak1" "$reason" "$*" >> "$W/RSS_LOG.tsv"
cat "$ST"
exit $rc
