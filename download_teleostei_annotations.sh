#!/bin/bash

outdir=teleostei_annotations_ncbi
zipfile=${outdir}.zip
infofile=info_${outdir}.tab

datasets download genome taxon teleostei --reference --annotated --include gff3 --dehydrated --filename $zipfile

dataformat tsv genome --package $zipfile --fields --fields accession,assminfo-level,assmstats-contig-n50,assmstats-scaffold-n50,organism-name,assminfo-submission-date,assminfo-status,annotinfo status,assminfo-type > $infofile

unzip -d $outdir $zipfile

datasets rehydrate --directory $outdir