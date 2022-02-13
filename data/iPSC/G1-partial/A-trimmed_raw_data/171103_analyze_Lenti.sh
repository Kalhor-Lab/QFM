echo >.tracking ""
for i in `ls ../0-raw_data_Lenti/*R1_001.fastq`; do
	prename="$(cut -d'/' -f3 <<<$i)"; name="$(cut -d'-' -f1 <<<$prename)"
        date=`date`
        counter=`squeue -u kreza1@jhu.edu | grep "kreza1@jhu.edu" | wc -l`
        size=`wc -l $i | awk '{print $1}'`
        runtimeinseconds=`expr $size / 1000 + 2`    # On average, last I checked, each 5000 FASTQ lines took about one second to process. To be safe, I gave 1 sec to 1000 lines
        echo submitting $i on $date. Total requested walltime is $runtimeinseconds seconds.  Total submitted counter is: $counter
        echo "#!/bin/bash
#SBATCH --job-name=$name
#SBATCH -p shared
#SBATCH -t 0-0:0:$runtimeinseconds
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=200mb
#SBATCH -o $name-%j.out
#SBATCH -e $name-%j.err
#SBATCH --mail-type=FAIL
#SBATCH --mail-use=kalhor@jhu.edu
perl -w 170530_genotype1_Lenti.pl $i" >.tempscript_Lenti
	sbatch .tempscript_Lenti
        sleep 2
        while [ `squeue -u kreza1@jhu.edu | grep "PD" | wc -l` -ge 30 ] 
        do
                time=`date`
                counter=`squeue -u kreza1@jhu.edu | grep "PD" | wc -l`
                echo >>.tracking "Pending counter is $counter, waiting... $time"; sleep 10
        done
done
