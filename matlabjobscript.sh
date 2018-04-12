#!/bin/bash
#SBATCH -t 24:00:00
#SBATCH -n 10
#SBATCH -N 1
#SBATCH -A physics-hi
#SBATCH --mail-user=paraj@umd.edu
#SBATCH --mail-type=END
#SBATCH --mem-per-cpu=1024
#SBATCH --share
#SBATCH --array=1-3
. ~/.profile
module load matlab

echo "%%% real space Hamiltonian of graphene with gauge field A(sin(wt),cos(wt)).
%see notability note: Floquet graphene circularly polarized.
clear
tic
A=1.26
%M=0.385
Marr1=[1.4,1.62,1.7]
%[0.29,0.295,0.3,0.305,0.31,0.315,0.32,0.325,0.33,0.335,0.34,0.345,0.35,0.355,0.36,0.365,0.37,0.375,0.38];%,0.385,0.39,0.395,0.4,0.405,0.41,0.415,0.42]
Marr=[Marr1(str2num(getenv('SLURM_ARRAY_TASK_ID')))]
w=(2*pi)
T=2*pi/w
tnn=1.8 %nearest neighbour coupling
tnnn=0.13*tnn; %next nearest neighbour coupling, next nearest neighbor coupling shifts dirac point in energy by 3t
Tdiv=100;
dt=T/Tdiv;
Lx=60% even
Ly=30%even
PBCx=0;
PBCy=1; % for PBC=1 for OBC =0 along y direction (zigzag)
disavmax=100;
%Vrand=0.1;
seedvalue=14; %15*(str2num(getenv('SLURM_ARRAY_TASK_ID')));
rng(seedvalue);
Vrandarr=[1.2];
tnn_disarr=0.0*Vrandarr; %same length as Vrandarr
tnnn_disarr=0.0*0.25*Vrandarr; %samelength  as Vrandarr.
fixedbound=-3.013%-pi%-2.906 % sets a bound at E=-12 where the gap is trivial
movingboundarr=[0]; % use moving bound to scan through the energy BZ.
energywidthtolerance=0.025; %sets energy window in which to average the IPR for each realization
JobID='$SLURM_JOB_ID';
">>floqgraphene-params-$SLURM_JOB_ID.m

##cat floqgraphene-params-$SLURM_JOB_ID.m cluster-floqgrapheneunitary.m > cluster-floqgrapheneunitary-$SLURM_JOB_ID.m
##cat floqgraphene-params-$SLURM_JOB_ID.m cluster-floqgrapheneunitary_disordered_hopping.m > cluster-floqgrapheneunitary_disordered_hopping-$SLURM_JOB_ID.m

cat floqgraphene-params-$SLURM_JOB_ID.m cluster_time_evolution_floqgrapheneunitary_disordered_hopping.m > cluster_time_evolution_floqgrapheneunitary_disordered_hopping-$SLURM_JOB_ID.m

##matlab -nosplash -nodesktop -nojvm < cluster-floqgrapheneunitary_disordered_hopping-$SLURM_JOB_ID.m > output-$SLURM_JOB_ID.txt

matlab -nosplash -nodesktop -nojvm < cluster_time_evolution_floqgrapheneunitary_disordered_hopping-$SLURM_JOB_ID.m > output-$SLURM_JOB_ID.txt
rm *$SLURM_JOB_ID.m

hostname
date

