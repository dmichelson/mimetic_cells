module load gcc/6.2.0
module load python/2.7.12
module load macs2/2.1.1.20160309
# module load python/3.6.0
# module load idr/2.0.2

/n/groups/cbdm-db/dam41/mec_scatac/peaks/2020-10-01/macs2

ls *.narrowPeak > /n/groups/cbdm-db/dam41/mec_scatac/peaks/2020-10-01/macs2/narrowPeak.idx

while read narrowPeak
do
    echo $narrowPeak
    sort -k8,8nr unsorted/$narrowPeak > sorted/$narrowPeak
done < narrowPeak.idx

for n in 1 2 3 4 5 6 7 8 9 10 11 12 13
do
    echo $n >> logs/run_log_2020_06_03.txt
    idr --samples sorted/WT1_KO1_merged_peaks${n}_peaks.narrowPeak sorted/WT2_KO2_merged_peaks${n}_peaks.narrowPeak \
        --input-file-type narrowPeak \
        --rank p.value \
        --output-file idr/${n}_idr_peaks \
        --plot \
        --log-output-file logs/${n}_idr.log
    echo "there are this many signficant peaks for sample " $n "..." >> logs/run_log_2020_10_11.txt
    awk '{if($5 >= 540) print $0}' idr/${n}_idr_peaks | wc -l >> logs/run_log_2020_10_11.txt
done