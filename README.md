# Gskimming
 Analyses of Nuclear Reads Obtained using Genome skimming

# Module list in Mjolnir 
```  1) skmer/3.3.0   2) apples/20231016   3) fasttree/2.1.11   4) mash/2.3   5) fastme/2.1.6.1  /projects/mjolnir1/people/sjr729/softwares/guppy-64 ```

#if CONSULT used, it needs a lot of memory
```
#salloc --qos=debug --nodes=1 -c 4 --mem-per-cpu 18000 -t 1400
srun --cpus-per-task=4 --mem=180G --time=2:00:00 --pty bash

srun --pty -n 1 -c 1 bash -i
 1) sratoolkit/3.0.0   2) sra-tools/3.0.3
```

#Phacochoerus dataset
```
accession_list=(
        SRR19174497
        SRR19174498
        SRR19174499
        SRR19174500
        SRR19174501
        SRR19174502
        SRR19174503
        SRR19174504
        SRR19174505
        SRR19174506
        SRR19174507
        SRR19174508
        SRR19174509
        SRR19174510
        SRR19174511
        SRR19174512
        SRR19174513
        SRR19174514
        SRR19174515
        SRR19174516
        SRR19174517
        SRR19174518
        SRR19174519
        SRR19174520
        SRR19174521
        SRR19174522
        SRR19174523
        SRR19174524
        SRR19174525
        SRR19174526
        SRR19174527
        SRR19174528
        SRR19174529
        SRR19174530
        SRR19174531
        SRR19174532
        SRR19174533
        SRR19174534
        SRR19174535
        SRR19174536
        SRR19174537
        SRR19174538
        SRR19174539
        SRR19174540
        SRR19174541
        SRR19174542
        SRR19174543
        SRR19174544
        SRR19174545
        SRR19174546
        SRR19174547
        SRR19174548
        SRR19174549
        SRR19174550
        SRR19174551
        SRR19174552
        SRR19174553
        SRR19174554
)
```
