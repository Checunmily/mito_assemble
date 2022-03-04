# mito_assemble
1.To use this script, first you need to prepare the pair-end reads files of each taxon.
like danio_rerio_R1.fq , danio_rerio_R2.fq ......


2.Then you should write a list of all the species' name into a .txt file and each of them must in a single line.
like this:
danio_rerio
name2
......



3.You also need to prepare a .fas file for BLASTN. This file is usually the 11 sites of sequences which you have chosen as a qeury.
like this:
>site1
AAAAATTTTCCCGGG
>site2
GGCCTTCCGTAG
......



4.After all the files above have been prepared, you can use the script to assemble the mitogenome of each taxon.
like this:

perl assemble_mito_novo.pl -- query taxon.txt --reference mitoref.fas



5.Don't forget to put data.pm module into the same dir.

