# MDPS

Author: Amlan Talukder

Date: July 16,2019

MDPS is used to filter out false positive targets from a given list of miRNA targets. 
It was developed by the computational System biology group at UCF.


INSTALLATION
-------------------------------------------------------------------------------------------
   (1) Install Python 2.7

EXECUTION 
-------------------------------------------------------------------------------------------

   (1) cd to the directory where mdps.py is located.
   (2) You can run the software by running the script "mdps.py" with the following command:
   
   ----------------------------------------------------------------------------------------
   
	usage: mdps.py [-h] -p PREDICTIONS -t TARGET_SEQUENCES -m MIRNA_SEQUENCES
               [-o OUTPUT]

	Filter predicted miRNA targets

	optional arguments:
	  -h, --help           show this help message and exit
	  -o OUTPUT            Path for MDPS outputs

	required arguments:
	  -p PREDICTIONS       Path for miRNA prediction file in a tsv format
	  -t TARGET_SEQUENCES  Path for target sequences in fasta file format
	  -m MIRNA_SEQUENCES   Path for miRNA sequences in fasta file format

	Example: python mdps.py -p examples/test_predictions.txt -t
	examples/test_target_seq.fa -m examples/test_mirna_seq.fa -o
	examples/test_output.txt

Required inputs
---------------------------------------------------------------------------------------------
The tool takes the miRNA prediction file path, corresponding target sequences and miRNA 
sequences as required input.

Input file format
---------------------------------------------------------------------------------------------
The prediction file must be in a tab delimited (tsv) format. Each line of this file must have
at least the following information in the first four columns.
	1. miRNA id
	2. Target id
	3. Target start (based on the provided target sequence)
	4. Target end (based on the provided target sequence)

Example of a prediction file information:

miR-744	ENSG00000117707_ENST00000366958_PROX1	331	351
miR-744	ENSG00000117707_ENST00000366958_PROX1	1175	1195

The target sequence and miRNA sequence files will be in fasta format. The target ids and the 
miRNA ids must be the same ones provided in the prediction file.

Example of a target sequence file information:

>ENSG00000117707_ENST00000366958_PROX1
AAGTAAATCTTGTTGTGGAGCGGAGCCCTCAGCTGAGGGAGCGCTCTGAAATAATACACCATTGCAGCCGGGGAAAGCAGAGCGGCGCAAAAG
AGCTCTCGCCGGGTCCGCCTGCTCCCTCTCCGCTTCGCTCCTCTTCTCTTCTTTACCCTTCTCCTCTCTCCTCCTCTGCTGCTCTCTCCTCTC
CTCCCGCTCTTCTCTCTCCTCCTCTCCTGCTCTCTCCTCTTCCCTTAGCTCCTCTTCTTTTCTTCTCCTCTTCTTCCCTCTCCTCGCCTCTCC
CCTGCTCCTCTTCTCTCGTCTCCCCTCCCCTCCCGCCTCTCTCTCCCCTCTCCCTCTCCCACTCGCCCCGCTCGCTCGCTCGCTGTCGCACAG
ACTCACCGTCCCTTGTCCAATTATCATATTCATCACCCGCAAGATATCACCGTGTGTGCACTCGCGTGTTTTCCTCTCTCTGCCGGGGGAAAA
AAAAGAGAGAGAGAGAGATAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGGCTCGGTCCCACTGCTCCCTGCACCGCGGTCCCGGGATTCTTGA
GCTGTGCCCAGCTGACGAGCTTTTGAAGATGGCACAATAACCGTCCAGTGATGCCTGACCATGACAGCACAGCCCTCTTAAGCCGGCAAACCA
AGAGGAGAAGAGTTGACATTGGAGTGAAAAGGACGGTAGGGACAGCATCTGCATTTTTTGCTAAGGCAAGAGCAACGTTTTTTAGTGCCATGA
ATCCCCAAGGTTCTGAGCAGGATGTTGAGTATTCAGTGGTGCAGCATGCAGATGGGGAAAAGTCAAATGTACTCCGCAAGCTGCTGAAGAGGG

Example of a miRNA sequence file information:

>miR-744
TGCGGGGCTAGGGCTAACAGCA

Optional inputs
---------------------------------------------------------------------------------------------
The output file path is an optional parameter.


MODEL
---------------------------------------------------------------------------------------------
The cost_threshold.txt and mirna_position_wise_knowledge.txt are used as model files to run
mdps.py.


RESULTS
---------------------------------------------------------------------------------------------
If the optional output file parameter is not provided the result file will be stored in the 
same directory where mdps.py is located by the name 'predictions_filtered.txt'. The result 
files are created in the same format as the input prediction file.

LICENSE & CREDITS
-------------------------------------------------------------------------------------------------
The software is a freely available for academic use.
plase contact xiaoman shawn li (xiaoman@mail.ucf.edu) for further information. 


CONTACT INFO
-------------------------------------------------------------------------------------------------
If you are encountering any problem regarding to MDPS, please refer the manual first.
If problem still can not be solved, please feel free to contact us:
Amlan Talukder (amlan@knights.ucf.edu)
Xiaoman Shawn Li (xiaoman@mail.ucf.edu)
Nancy Haiyan Hu (haihu@cs.ucf.edu)
