# NextVpower
Next V-Power, adapted from "Virus Phylogenetic Resolver of Wastwater-based Epidemiology (V-Power).
=======
Next V-Power is a tool for SARS-CoV-2 lineage demixing via amplicon sequencing data.

Requirements:
---------------
 - python3
 - numpy
 - pandas
 - scipy

This tool was developed on Windows, and tested on Linux.

Installation:
---------------
1.Clone this repository to your local directory. 
2.Install requiments in a python environment. Skip if you have got these requirments.
```sh
pip install numpy pandas scipy
```
3.Change directory to where you cloned this repository to.
```sh
cd /path/of/NextVPower
```
4.Extract `usher_barcodes.zip` and `var_anno.zip` under current folder.
```sh
unzip usher_barcodes.zip
```
```sh
unzip var_anno.tsv
```
5.Run python command to check installation and see help.
```sh
python NextVpower.py -h
```

Usage: 
---------------
NextVpower.py -i INPUT [-o OUTPUT] [-b BARCODE] [-l LINEAGES] [-r MINRATE] [-m] [-v] [-d MINDEPTH]
              [-a ANNFILE] [--barfilter2] [--k_lineages K_LINEAGES] [--n_sites N_SITES]
              [--ann_outpath ANN_OUTPATH] [--vcsample VCSAMPLE] [--fsample FSAMPLE] [--fbarcode FBARCODE]
              [--potentials POTENTIALS] [-h] [--version]

#options:
  -h, --help            show this help message and exit
  
  --version             show program's version number and exit

#Demixing arguments for the demix solver:

  -i INPUT, --input INPUT	[File/Dir] path of input sample table file or vcfs folder
  
  -o OUTPUT, --output OUTPUT	[File] path of output table file (default: ./demix_result.tsv)
  
  -b BARCODE, --barcode BARCODE	[File] specify a usher_barcodes.csv as input barcode matrix (default: ./usher_barcodes.csv)
  
  -l LINEAGES, --lineages LINEAGES	[Int] maximum number of demixing lineages (default: 100)

#Sample Processing arguments for the [--input] sample handler:

  -r MINRATE, --minrate MINRATE	[Int] filter mutation sites with mutation rate lower than setting threshold in sample vectors (default: 0)
  
  -m, --merge           [Flag] merge lineages with completely identical mutation sites in the barcode matrix
  
  -v, --vcfs            [Flag] parse *.vcf files under input folder
  
  -d MINDEPTH, --mindepth MINDEPTH	[Int] filter mutation sites with depth lower than setting threshold in *.vcf files (default: 0)
  
  -a ANN_FILE, --ann_file ANN_FILE	[File] specify a var_anno.tsv as input variation annotation table (default: ./var_anno.tsv)

#Barcode Processing arguments for the barcode filter:

  --barfilter2          [Flag] use another barcode filter (authored by Kun Yang) to handle barcode matrix
  
  --k_lineages K_LINEAGES	[Int] filter lineages with fewer than [Int] mutation sites (default: 200)
  
  --n_sites N_SITES     [Int] retain "key" mutation sites present in more than [Int] lineages (default: 20)

#Middle file output arguments for middle processes:

  --ann_outpath ANN_OUTPATH	[Dir] if not None, add save the annotated *.vcf table files under a folder (optional)
  
  --vcsample VCSAMPLE   [File] save the sample table file converted from *.vcf files (optional)
  
  --fsample FSAMPLE     [File] save the filtered sample table file (optional)
  
  --fbarcode FBARCODE   [File] save the filtered barcode matrix file (optional)
  
  --potentials POTENTIALS	[File] save potential sites not recorded in barcode but present in samples (optional)

Example:
---------------
1. Demix from input sample table file, and save result in result.tsv:
```sh
python NextVpower.py -i PP_raw_example.tsv -o demix_result_example.tsv
```
2. Demix from input *.vcf files under a folder, retain top 200 lineages, and save result in result.tsv:
```sh
python NextVpower.py -i vcf_example -v -l example -o demix_result_vcf_example_top200.tsv
```
3. Merge lineages with completely identical mutation sites when Demixing:
```sh
python NextVpower.py -i PP_raw_example.tsv -m -o demix_result_example_merged.tsv
```
4. Add annotation to *.vcf files according to variation annotation table and demix:
```sh
python NextVpower.py -i vcf_example -v -o demix_result_vcf_example.tsv --ann_outpath ann_tab_example
```
5. Demix from input *.vcf files under a folder, save result in result.tsv, and save middle data to files:
```sh
python NextVpower.py -i vcf_example -v -o demix_result_vcf_example.tsv --vcsample PP_raw_example.tsv --fbarcode MMFF_example.tsv --fsample PPFF_example.tsv
```

Publications
------------
This project was not published yet —— but you can still have a try on your SARS-CoV-2 amplicon data.

License
-------
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

