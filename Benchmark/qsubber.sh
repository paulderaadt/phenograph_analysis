#!/bin/bash 
#$ -S /bin/bash  
#WELKE SHELL
#$ -q all.q
#WELKE SERVER
#$ -N FullBM 
#NAAM JOB
#$ -cwd 
#OUTPUT FILES NAAR CWD
#$ -j Y
# JOIN STDERROR EN STDOUT
#$ -V 
#GEBRUIK ENV
#$ -m e 
#MAIL OP EINDE b VOOR BEGIN
#VALIDE EINDMAIL
#$ -l h_vmem=8G
#1 GIGABYTE MAX GEHEUGEN


echo Start time : `date` 
echo 'Starting job...' 

echo python main.py "$@"
echo 'end of script.' 
echo End time : `date`
