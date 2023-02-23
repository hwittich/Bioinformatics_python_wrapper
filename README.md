# Bioinformatics Python Wrapper
This project was created for COMP 383 - Bioinformatics, taught by Dr. Wheeler at Loyola University Chicago. The goal of the assignment was to learn how to construct a wrapper for running data through a bioinformatics pipeline. The bioinformatics tasks involved in this pipeline include quantifying reads from RNA Seq data using kallisto and sleuth, alignment with bowtie, assembly with SPAdes, and querying with BLAST.

---
## Languages and Packages
---
**Python**
```
import os
import argparse
import sys
from Bio.Seq import Seq
from Bio import Entrez
from Bio import SeqIO
```
**R**
```
library(sleuth)
library(dplyr)
```

## Software ##
---
- [Kallisto](https://pachterlab.github.io/kallisto/manual)
- [Sleuth](https://pachterlab.github.io/sleuth/manual)
- [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
- [SPAdes](https://cab.spbu.ru/files/release3.13.0/manual.html)
- [BLAST+](https://www.ncbi.nlm.nih.gov/books/NBK279691/)

## Running the Pipeline ##
* Clone this repository to your machine
	```
	git clone https://github.com/hwittich/COMP383_mini_project.git
	```

* Move into the COMP383_mini_project directory
	```
	cd COMP383_mini_project
	```

* Run the python wrapper, optional argument -t or --test. When this flag is specified, the program will run on the test data included in the data folder, writing outputs to the test_outputs directory. Otherwise, the program will download larger data files to run and write outputs to the miniProject_Henry_Wittich directory.   
	**Test run option 1:**
	```
	python3 pipeline.py -t
	```
	**Test run option 2:**
	```
	python3 pipeline.py --test
	```
	**Running with full dataset:**
	```
	python3 pipeline.py
	```
