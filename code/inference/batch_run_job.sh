# std data + MCMC + DIVA-clado
tmux new-session -d -s "anolis_nr9_2000" ./single_run_job.sh anolis_nr9_ns383 2000 true false 4 5;sleep 3
tmux new-session -d -s "anolis_nr9_2001" ./single_run_job.sh anolis_nr9_ns383 2001 true false 4 5;sleep 3
tmux new-session -d -s "anolis_nr9_2002" ./single_run_job.sh anolis_nr9_ns383 2002 true false 4 5;sleep 3

# std data + RJMCMC + DIVA-clado
tmux new-session -d -s "anolis_nr9_2100" ./single_run_job.sh anolis_nr9_ns383 2100 true true 4 5;sleep 3
tmux new-session -d -s "anolis_nr9_2101" ./single_run_job.sh anolis_nr9_ns383 2101 true true 4 5;sleep 3
tmux new-session -d -s "anolis_nr9_2102" ./single_run_job.sh anolis_nr9_ns383 2102 true true 4 5;sleep 3

# thin data + MCMC + DIVA-clado
tmux new-session -d -s "anolis_nr9_2200" ./single_run_job.sh anolis_nr9_ns378_major_edits 2200 true false 4 6;sleep 3
tmux new-session -d -s "anolis_nr9_2201" ./single_run_job.sh anolis_nr9_ns378_major_edits 2201 true false 4 6;sleep 3
tmux new-session -d -s "anolis_nr9_2202" ./single_run_job.sh anolis_nr9_ns379_major_edits 2202 true false 4 6;sleep 3

# thin data + RJMCMC + DIVA-clado
tmux new-session -d -s "anolis_nr9_2300" ./single_run_job.sh anolis_nr9_ns378_major_edits 2300 true true 4 6;sleep 3
tmux new-session -d -s "anolis_nr9_2301" ./single_run_job.sh anolis_nr9_ns378_major_edits 2301 true true 4 6;sleep 3
tmux new-session -d -s "anolis_nr9_2302" ./single_run_job.sh anolis_nr9_ns378_major_edits 2302 true true 4 6;sleep 3

