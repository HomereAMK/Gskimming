## % abundance file KRANK
```bash
FILE_NAME=SRR26417692_2.fastq.gz
CLASSFILE_NAME_PATH=/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/skmer1/echarvel_oedulis_skimpip_results_17.08_nonrandlib_tresh0.2_0.2/krank_output/krank_reports


grep 'kingdom' ${CLASSFILE_NAME_PATH}/classification_info-${FILE_NAME} | column -t
```



## correct_classification.sh
```bash
set -x
grep 'Homo sapiens' ${1} | awk -v var=${2} -F"\t" '$6>50' | cut -f1 > ${1}-corrected

grep -v 'U' ${1} | grep -v 'Homo sapiens' | grep -v 'SEQ_ID' |  awk -v var=${3} -F"\t" '$6>var' | cut -f1 >> ${1}-corrected
```

## run correct_classification.sh
```bash
#var
FILE_NAME=SRR26417692_2.fastq.gz
NEW_BACTERIAL_THRESHOLD=0.2
NEW_BACTERIAL_THRESHOLD=0.2
CLASSFILE_NAME_PATH=/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/skmer1/echarvel_oedulis_skimpip_results_17.08_nonrandlib_tresh0.2_0.2/krank_output/krank_reports

#run correct_classification.sh
bash /projects/mjolnir1/people/sjr729/Skmer_ms/scripts/correct_classification.sh \
${CLASSFILE_NAME_PATH}/classification_info-${FILE_NAME} \
${NEW_BACTERIAL_THRESHOLD} \
${NEW_HUMAN_THESHOLD}
```

## 
```bash
RAW_FASTQ_FOLDER=/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/oedulis/fastq/oedulis_RawFastq_14.08.24
FILE_NAME=SRR26417692_2.fastq.gz
CLASSFILE_NAME_PATH=/projects/mjolnir1/people/sjr729/Skmer_ms/Oedulis/skmer1/echarvel_oedulis_skimpip_results_17.08_nonrandlib_tresh0.2_0.2/krank_output/krank_reports

/projects/mjolnir1/people/sjr729/skimming_scripts-echarvel/bbmap/filterbyname.sh \
in=${RAW_FASTQ_FOLDER}/${FILE_NAME} \
names=${CLASSFILE_NAME_PATH}/classification_info-${FILE_NAME}-corrected \
out=./new_${FILE_NAME} \
-include=true #change this for include=false