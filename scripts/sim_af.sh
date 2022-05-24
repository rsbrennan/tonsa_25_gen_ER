
for pop in aa ah ha hh; do

    cd ~/tonsa_genomics/analysis/cmh_simulation/${pop}

    ls *.sync |  parallel -j 12 ' perl ~/bin/popoolation2/cmh-test.pl --input {} --output {}.cmh --min-count 2 \
        --min-coverage 2 --max-coverage 50000 --remove-temp \
        --population 1-2,3-4,5-6,7-8'

# ls *f25.*.sync |  parallel -j 10 ' perl ~/bin/popoolation2/cmh-test.pl --input {} --output {}.cmh --min-count 2 \
#    --min-coverage 40 --max-coverage 50000 --remove-temp \
#    --population 1-2,3-4,5-6,7-8'

# parse these files down so we just have the last column, and first 2 columns

    for i in $(ls *.sync.cmh | cut -f 1-2 -d '.'); do

        awk '{OFS="\t"} {print $1":"$2,$12}' ${i}.sync.cmh > ${i}.pval

    rm ${i}.sync.cmh

    done

done
#for i in $(ls *f25.*.sync | cut -f 1-3 -d '.'); do

#  cut -f 12 ${i}.sync.cmh > ${i}.pval

#done

#rm simulated.f25.*.sync.cmh
#rm *sync.cmh
