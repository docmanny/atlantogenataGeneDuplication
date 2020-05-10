#!/bin/bash

luigi --module cat RunCat \
--hal=../s3-multipart/halAnc105root.hal \
--ref-genome=Homo_sapiens --workers=4 --config=basic_human_chimp.config \
--work-dir humanchimp4 --out-dir humanchimp_out --local-scheduler \
--target-genomes='("Pan_troglodytes",)' --maxCores=4 > log3.out 2> log3.err
