cd ~/tonsa_genomics/analysis/

# have already done the filtering, so settings here are very liberal.
perl ~/bin/popoolation2/cmh-test.pl --input variants.sync --output AA.cmh --min-count 1 \
    --min-coverage 2 --max-coverage 5000 --remove-temp \
    --population 1-5,2-6,3-7,4-8

perl ~/bin/popoolation2/cmh-test.pl --input variants.sync --output AH.cmh --min-count 1 \
    --min-coverage 2 --max-coverage 5000 --remove-temp \
    --population 1-9,2-10,3-11,4-12

perl ~/bin/popoolation2/cmh-test.pl --input variants.sync --output HA.cmh --min-count 1 \
    --min-coverage 2 --max-coverage 5000 --remove-temp \
    --population 1-13,2-14,3-15,4-16

perl ~/bin/popoolation2/cmh-test.pl --input variants.sync --output HH.f3.cmh --min-count 1 \
    --min-coverage 2 --max-coverage 5000 --remove-temp \
    --population 1-21,2-22,3-23,4-24

perl ~/bin/popoolation2/cmh-test.pl --input variants.sync --output HH.cmh --min-count 1 \
    --min-coverage 2 --max-coverage 5000 --remove-temp \
    --population 1-25,2-26,3-27,4-28

# run cmh between f25 samples

perl ~/bin/popoolation2/cmh-test.pl --input variants.sync --output AH.f25.cmh --min-count 1 \
    --min-coverage 2 --max-coverage 5000 --remove-temp \
    --population 5-9,6-10,7-11,8-12

perl ~/bin/popoolation2/cmh-test.pl --input variants.sync --output HA.f25.cmh --min-count 1 \
    --min-coverage 2 --max-coverage 5000 --remove-temp \
    --population 5-13,6-14,7-15,8-16

perl ~/bin/popoolation2/cmh-test.pl --input variants.sync --output HH.f25.cmh --min-count 1 \
    --min-coverage 2 --max-coverage 5000 --remove-temp \
    --population 5-25,6-26,7-27,8-28
