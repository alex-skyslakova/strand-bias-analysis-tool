# Strand Bias Analysis Tool

## Overview
SBAT is a Python command-line tool for detection of strand bias. 
Strand bias is a situation when information from one strand of DNA is overrepresented compared to the information from 
the other strand. It is one of the types of bias that occur in next-generation sequencing data. Strand bias might lead 
to incorrect evaluation of results gained from sequencing data, if the bias is high. 
This tool offers a way of validating quality of the data against strand bias. More about strand bias and development 
of this tool can be found  [here](path to bachelor thesis once it is public). 

The tool uses [Jellyfish](https://github.com/gmarcais/Jellyfish) k-mer counting tool for counting k-mers in the NGS data and compares frequencies of k-mers 
and their complements, creating both statistics and visual analysis of the results of strand bias. 


## Installation 
#### First, <b>Jellyfish</b> must be installed. 
On Debian and Ubuntu with `apt`:
```Shell
sudo apt update
sudo apt install jellyfish
```

On MacOS with `brew`:
```Shell
brew install jellyfish
```

On Arch, it is available from [AUR](https://aur.archlinux.org/packages/jellyfish/).

On Windows, the best option is to use WSL. For other OS or installation from source code, see [here](https://github.com/gmarcais/Jellyfish)

#### After Jellyfish is installed, proceed with SBAT itself:
Installation from pip
```Shell
pip install sbat
```

To install from source code, download the code and run following in the root of the source tree:
```Shell
python3 -m pip install --upgrade build
python3 -m build
pip install -e .
```

## Usage

In order to perform analysis on one or multiple files, use command ```sbat``` followed by your files. By default, k-mers for k in 
range 5 to 10 are analysed. The figures and statistics of the tool are saved into ```sbat_out``` directory, which will be created if does not exist already. By default all the partial results of the SBAT (Jellyfish-generated files, statistics for each size of k) are deleted at the end of the run. To prevent that,
use argument ```-c```
```Shell
sbat my_file.fasta my_file2.fasta my_file3.fastq
```

Following command additionally specifies output directory with ```-o``` and keeps partial results of computations 
using parameter ```-c```. To speed up SBAT run time, use parameter ```-t T``` with specified number of threads
you wish to pass to the application. To specify size of k-mers for which you want to run analyses, use parameter ```-m START END```. 
If one argument is passed to it, SBAT runs only for this size of k. If two arguments are passed, application analyses 
k-mers in range [START, END]
```Shell
sbat my_file.fasta my_file2.fasta my_file3.fastq -o output_dir -c -t 10 -m 5 8
```

If you want to analyse Nanopore dataset, add ```-n``` in order to run more specific, time-based analysis. As part of this
analysis, dataset is divided into one-hour long bins. Each of them is then analysed on its own. The time duration of 
one bin can be set by ```-i H``` parameter followed by number of hours. If you wish to subsample your data, you can use 
parameters ```-r N``` or ```-b N``` to take only first N reads or bases of each bin. 
```Shell
sbat my_nanopore.fastq -o output_dir -b 500M -i 4 -n
```

To see all possible options, run:
```Shell
sbat -h
```

From version 0.0.9, ```-p``` parameter enables creation of interactive plots as well as .jpg results. After analysis 
finishes, SBAT creates Bokeh server on http://localhost:5006/
