
## installation of Skmer2
```bash
~/skmer-2-test
conda install skmer==3.2.1
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge


pwd = /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea

skmer --debug reference /input_directory/ -l ./library/ -r /path/to/reference.fna -p threads -o output
GENOME="/projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/genomeClupea/ncbi_dataset/data/GCA_900700415.2/GCA_900700415.2_Ch_v2.0.2_genomic.fna"
DATE=$(date +%d.%m)

#
skmer --debug reference skims_processing_pipeline_jan24/kraken/ -r $GENOME -p 5 -o skmer2_kakenfq_$DATE



#cd /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/

conda activate skmer_2_test
conda install jellyfish seqtk mash
python /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/Skmer-2/skmer/__main_TESTING.py --debug reference /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/skims_processing_pipeline_jan24/kraken/ -r $GENOME -p 15 -o skmer2_kakenfq_$DATE -l /projects/mjolnir1/people/sjr729/tutorial/skimming_scripts/testClupea/library_skmer2_$DATE


python Skmer-2/skmer/__main_TESTING.py --debug reference ../skimming_scripts/testClupea/skims_processing_pipeline_jan24/consult/ -r $GENOME -p 5 -o skmer2_kakenfq_$DATE -l ../skimming_scripts/testClupea/library_skmer2

``` 