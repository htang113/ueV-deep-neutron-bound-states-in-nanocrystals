#!/bin/bash 
#SBATCH -J vasp          
#SBATCH -o vasp.%j.out     
#SBATCH -e vasp.%j.err 
#SBATCH -n 112         
#SBATCH -N 2 
#SBATCH -p small      
#SBATCH -t 4:00:00        
#SBATCH -A DMR20023

module load vasp/5.4.4
ibrun vasp_std > vasp_test.out
