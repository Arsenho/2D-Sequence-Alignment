#!/bin/bash
#SBATCH --job-name=SEQUENTIAL-150
#SBATCH --partition=normal            # submission queue (normal or bigmem)
#SBATCH --time=02-00:00:00            # 1-1 means one day and one hour
#SBATCH --mail-type=ALL    # Type can be BEGIN, END, FAIL, ALL(any state change).
#SBATCH --mail-user=arsenenanfackvadele@yahoo.fr # e-mail notification
#SBATCH --output=job_seq-%j.out        # if --error is absent, includes also the errors
#SBATCH --nodes=1                   # Number of nodes
#SBATCH --ntasks=1                # Number of MPI ranks
#SBATCH --ntasks-per-node=1 #Number of MPI ranks per node
#SBATCH --cpus-per-task=1            # Number of cores per MPI rank 
#SBATCH --distribution=cyclic:cyclic # Distribute tasks cyclically on nodes and sockets
#SBATCH --mem-per-cpu=4G          # Memory per processor

echo "----"
echo "hostname                    =   $(hostname)"
echo "SLURM_JOB_NAME              =   $SLURM_JOB_NAME"
echo "SLURM_SUBMIT_DIR            =   $SLURM_SUBMIT_DIR"
echo "SLURM_JOBID                 =   $SLURM_JOBID"
echo "SLURM_JOB_ID                =   $SLURM_JOB_ID"
echo "SLURM_NODELIST              =   $SLURM_NODELIST"
echo "SLURM_JOB_NODELIST          =   $SLURM_JOB_NODELIST"
echo "SLURM_TASKS_PER_NODE        =   $SLURM_TASKS_PER_NODE"
echo "SLURM_JOB_CPUS_PER_NODE     =   $SLURM_JOB_CPUS_PER_NODE"
echo "SLURM_TOPOLOGY_ADDR_PATTERN =   $SLURM_TOPOLOGY_ADDR_PATTERN"
echo "SLURM_TOPOLOGY_ADDR         =   $SLURM_TOPOLOGY_ADDR"
echo "SLURM_CPUS_ON_NODE          =   $SLURM_CPUS_ON_NODE"
echo "SLURM_NNODES                =   $SLURM_NNODES"
echo "SLURM_JOB_NUM_NODES         =   $SLURM_JOB_NUM_NODES"
echo "SLURMD_NODENAME             =   $SLURMD_NODENAME"
echo "SLURM_NTASKS                =   $SLURM_NTASKS"
echo "SLURM_NPROCS                =   $SLURM_NPROCS"
echo "SLURM_MEM_PER_NODE          =   $SLURM_MEM_PER_NODE"
echo "SLURM_PRIO_PROCESS          =   $SLURM_PRIO_PROCESS"
echo "----"

# special commands for openmpi/intel
module load openmpi/intel-opa/gcc/64/1.10.4-hfi

export CC=/cm/local/apps/gcc/8.2.0/bin/gcc && export CXX=/cm/local/apps/gcc/8.2.0/bin/g++

#make clean all

#pip install mpi4py

python3 serial.py --rows 150 --columns 150