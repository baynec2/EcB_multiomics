#!/usr/bin/env bash

# Copying final figures to directory 
cp 02_figures/**/Figure*.pdf 07_final_manuscript/

# Copying supplemental 
## Figures
 cp 03_supplemental/**/S[[:digit:]].pdf 07_final_manuscript/supplemental/

## Tables
cp 03_supplemental/tables/*.xlsx 07_final_manuscript/supplemental/
