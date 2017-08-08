# Scripts for .bashrc on torquelord
#
# To use, run the following:
# echo "source nanopore-scripts/torquelord_scripts.sh" >> ~/.bashrc
#
# Commands:
# q: views all torque jobs under your usernam
# qsub: modified to send all output to a folder ~/torque_logs/
# host: ssh into the host of a job, or list all nodes if there are more than one
# stdout: print stdout of an ongoing job
# stderr: print stderr of an ongoing job
# qalter: change the priority of a queued job (must be in queue: node_waiting and not yet be building. Will only affect priority between jobs you own)
alias q="qstat -u $USER"
mkdir -p ~/torque_logs/
alias qsub="qsub -o ~/torque_logs/ -e ~/torque_logs/"
host () {
  PROC=$1
  if [ $(q -rnt1 $PROC | tail -n +6 | sed "s/.* //" | tr + "\n" | wc -l) -eq 1 ]
    then
    ssh `q -rnt1 $PROC | tail -n1 | sed "s/.* //"`
  else
    q -rnt1 $PROC | tail -n +6 | sed "s/.* //" | tr + "\n"
  fi
}; export -f host
# % host
# torque_big-95435-n1.hpc.wehi.edu.au
# torque_big-95447-n1.hpc.wehi.edu.au
# torque_big-95453-n1.hpc.wehi.edu.au
# % host 95345
# runs `ssh torque_big-95435-n1.hpc.wehi.edu.au`
stdout () {
  PROC=$(echo "$1" | sed 's/\[\([0-9]*\)\]/-\1/g')
  cat /stornext/HPCScratch/.torque/spool/$PROC.$HOSTNAME.OU
}; export -f stdout
stderr () {
  PROC=$(echo "$1" | sed 's/\[\([0-9]*\)\]/-\1/g')
  cat /stornext/HPCScratch/.torque/spool/$PROC.$HOSTNAME.ER
}; export -f stderr
# stderr 95345 | less
# stdout 95345 > temp.log
qalter () {
  OPTIND=1
  PRIORITY=1023
  ATTEMPTS=5
  WAIT_TIME=10
  HELP="qalter (modified for Milton)
Scott Gigante, 2017-06-15
Usage: qalter [-h] [-p PRIORITY] [-t WAIT_TIME] [-n NUM_ATTEMPTS] PROC_ID"
  local arg
  while getopts 'hp:n:t:' arg
  do
    case ${arg} in
      h) echo "$HELP"; return 0;;
      p) PRIORITY=${OPTARG};;
      n) ATTEMPTS=${OPTARG};;
      t) WAIT_TIME=${OPTARG};;
      *) echo "$HELP" >&2; return 1 # illegal option
      esac
  done
  shift $(($OPTIND - 1))
  for PROC in $@; do
    PBS_ID=$(curl -ks "https://hpc1-pbs.wehi.edu.au/idmapping/get_pbs_id.php?headnode=$HOSTNAME&jobid=$PROC")
    if [ $( echo $PBS_ID | wc -c ) -eq 1 ]; then
      i=1
      while [ $( echo $PBS_ID | wc -c ) -eq 1 -a $i -lt $ATTEMPTS ]; do
        echo "PBS node ID request failed... (attempt $i of $ATTEMPTS) Trying again in 10 seconds."
        sleep $WAIT_TIME
        PBS_ID=$(curl -ks "https://hpc1-pbs.wehi.edu.au/idmapping/get_pbs_id.php?headnode=$HOSTNAME&jobid=$PROC")
        i=$(($i+1))
      done
      if [ $( echo $PBS_ID | wc -c ) -eq 1 ]; then
        echo "PBS ID request failed for job $PROC. Skipping."
        continue 2
      fi
    fi
    ssh hpc1-pbs qalter -p $PRIORITY $PBS_ID
    echo "Job $PROC (PBS ID $PBS_ID) set to priority $PRIORITY"
  done
}; export -f qalter
# % qalter 42330
# Job 42330 (PBS ID 101048) set to priority -1023
# % qalter -p 100 42330
# Job 42330 (PBS ID 101048) set to priority 100
