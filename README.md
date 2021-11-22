# RNASeq-Pipeline
Simplified Bioinformatics Pipeline for Scientists Analyzing RNASeq Data

## Table of Contents
- [Setting Up Your Work Environment](#environment)
  - [Configuring Your Github Account in Terminal](#configuring-your-github-account-in-terminal)
  - [Using this Repository](#using-this-repository)
- [Working in a Cloud Environment](#cloud-computing)
  - [Vanderbilt ACCRE Cluster Cloud Computing](#vanderbilt-accre)
- [Transcript Quantification with Salmon](#salmon)
- [Quality Checking with FastQC Analysis](#fastqc)
- [Differential Gene Expression Analysis](#differential-gene-expression-analysis)
  - [DESeq2](#deseq2)
  - [Plots](#plots)
    - [Principle Component Analysis (PCA)](#principle-component-analysis)
    - [Heat Map](#heat-map--component-analysis--)
  - [Normalized Counts Data](#normalized-counts)
    - [MA Plots](#ma-plots)
    - [Heat Maps](#heat-maps)
  - [Pairwise Analysis](#pairwise-analysis)
    - [Volcano Plots](#volcano-plots)
  - [Gene Set Enrichment Analysis (GSEA)](gene-set-enrichment-analysis)

## Environment
> How to set up the basic work environment for the RNASeq Pipeline,

### Configuring Your Github Account in Terminal

1. Make a free Github Account at github.com
2. Set global configuration values for your computer
option 1: type these commands in terminal
```shell
git config --global user.name "YOUR NAME"                  
git config --global user.email "YOUR-EMAIL@DOMAIN.COM"      

git config --list                                           # check to see all values have updated correctly
```
option 2: edit the `~/.gitconfig` file directly
```shell
[user]
    name = Your Name
    email = youremail@yourdomain.com
```

### Using this Repository
1. Clone this repository
Open a terminal window and navigate to the place you would like to download this package using the `cd` command.
Then, clone the repository using one of 2 commands with the link in the green code button at the top left of this repo.
```bash
git clone https://github.com/parkdj1/RNASeq-Pipeline.git    # HTTPS version (login using email/password)
git clone git@github.com:parkdj1/RNASeq-Pipeline.git        # SSH   version (login using an ssh key pair)
```
Alternatively, download the zip file where you want the folder to be and open it to unzip.

2. Python
3. R


## Cloud Computing
> Working in a cloud computing environment 101

### Basic Linux
Some computing resources have a GUI to easily access the cloud environment. Learning some basic linux commands can help you to really control what you're doing and opens up the possibilities for more functionality.  

Here are some basic command-line tools to help you get started in terminal.
> Tip: Use the `TAB` key to fill in the rest of your file path or see options for doing so  

```shell
# This is a comment. Everything following the `#` symbol is ignored by the computer

# Connecting to a Remote Host (Cloud Computer)
ssh username@login.source                                   # Connect / Log In : Your provider will likely provide these credentials for you.
                                                            # Doing so basically allows you to connect to someone else's computer from your own.
hostname (-I)                                               # Display hostname. Using the '-I' flag displays your IP address

# Transferring Files to/from a Remote Computer
scp fileName.ext username@login.source:/path/to/dest        # Upload files to a remote machine (use from your computer, not the remote host)
scp username@login.source:/path/to/file.ext /path/to/dest   # Download files from a remote machine using SCP

# Navigating your files and directories
pwd                                                         # 'print working directory' print which folder you are currently in
cd /path/to/folder                                          # 'change directory' move into the specified folder
cd ..                                                       # Move to parent folder (i.e. /home/parent/current -> /home/parent)

ls folder                                                   # List the contents of a folder (-l more details  -a list hiddent content too)
ls .                                                        # List the contents of the current folder

# Manipulating files and directories
mkdir folderName                                            # Create a folder with the specified name
touch fileName                                              # Create a file with the specified name

rm fileName                                                 # Delete file
rm prefix*                                                  # Delete all files that start with the prefix
rm -r folderName                                            # Delete folder, including all content within

mv fileName /path/to/dest                                   # Move a file or folder the desired location
mv fileName newFileName                                     # Rename file or folder
cp fileName /path/to/dest                                   # Copy file to another directory

wget link                                                   # Download files from the specified internet link

cat fileName                                                # Print the contents of the file to terminal
head fileName                                               # View the first 10 lines (change the number of lines using the `-n #` flag)
grep something file > newFile                               # Search for occurrences of 'something' in the file and store the output in a new file
diff file1 file2                                            # Output the lines that do not match across the specified files

# Helpful Tips
man command                                                 # Pull up the manual pages for the specified command
command -h or --help                                        # Show how to use command
clear                                                       # Clear up your terminal window. You can still scroll up to see your command history!

# Permissions and Control
sudo command                                                # "SuperUser Do" Use administrative or root privileges for a command
chmod options fileName                                      # Change the permissions for a file and make it executable. (options: +/- r/w/x OR ###)
chown owner[:group] fileName                                # Change the user and/or group ownership of a file
```

### Bash Scripts
> In order to submit work for the computers to run in the cloud, you need to submit a Bash Script.  
1. Create a script
```shell
touch scriptname.bash
```
3. File Header

A shebang (#!) is necessary to specify the type of file / language used
Specify different options with #SBATCH

For example :  
```shell
#!/bin/bash
#SBATCH --mail-user = email-username@vanderbilt.edu
#SBATCH --mail-type=ALL
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=16G
#SBATCH --time=300
```
3. Send the script `sbatch scriptname.bash`
4. Check the Status of Your Scripts `squeue -u username`


### Vanderbilt ACCRE
> Vanderbilt University's Cluster Computing Resource - ACCRE

ACCRE offers 2 types of clusters: a traditional cluster and the Jupyter cluster.  
For the rest of this section, I will mostly be referring to the traditional compute cluster.  

There are 2 ways to access ACCRE.

1. [The ACCRE Visualization Portal](http://portal.accre.vanderbilt.edu/)  
This is a graphical interface, meaning you can access the resource through your web browser and click buttons to navigate.
Find more information on [navigating the ACCRE visualization Portal](https://www.vanderbilt.edu/accre/portal/) on the Vanderbilt ACCRE Website.  

2. Shell access
The other way to access ACCRE is by connecting to the resource through your terminal.

```shell
ssh vunetid@login.accre.vanderbilt.edu                 # ACCRE Login
```

The 1st time you login, the system may prompt you asking if you trust the system. Type "yes" and hit return.  

You will also need to reset your password. You should receive a temporary password in an email from ACCRE.  
Once you log in the first time, enter the following command and follow the instructions to set a new password.  
```shell
accre_password change
```
The password must pass the [ZXCVBN test](https://www.bennish.net/password-strength-checker/) with a perfect score.

Find more commands on the [ACCRE Cheat Sheet](https://cdn.vanderbilt.edu/vu-wp0/wp-content/uploads/sites/157/2018/11/08171621/ACCRE-Cheat-Sheet-November-2018.pdf)

## Salmon
> [Salmon](https://salmon.readthedocs.io/en/latest/index.html) is a transcript quantification tool for RNA-seq data.  
> Because the RNASeq files are very large, we utilized ACCRE for increased computing power and decreased computing time.  

**Input**   : FASTA/FASTQ read files

**Output**  : quantification log file, trimmed read files

### Installation

There are a few different ways to do this in the documentation, but I found this the simplest way.
```shell
cd DESIRED_LOCATION
git clone https://github.com/COMBINE-lab/salmon.git
cd salmon
make
make install
```

### Transcriptome Index Set Up
1. Find a transcriptome index file and download it  
We used [this mus musculus index file from ensembl](ftp://ftp.ensembl.org/pub/release-
99/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz)
```shell
wget TRANSCRIPTOME_INDEX_FILE_URL [-o FILENAME]     # use the `-o` flag to rename the file when it downloads. Do not include the brackets if you decide to include this part!
```
2. Unzip the file
```shell
gunzip filename.gz
```

This file will be used for the `-i` flag in the main salmon command
```shell
salmon index -t filename.fa -i indexname
```

### Post-Run : Verification & Exports  
Export the `quant.sf` file and `log/*.log` files for further analysis

To check Salmon ran properly:  
1. Compare the number of reads in the FQC file to the number of reads in the trimmed file
2. Check the mapping rate (Our mapping rates were in the low to mid 80s)
3. Check read counts and ratios of genes of interest

## FastQC
Download FASTQC software (I cloned to ACCRE)

```
/path/to/FastQC/fastqc /path/to/file_1.fq.gz /path/to/file_2.fq.gz
```

## Differential Gene Expression Analysis

### DeSeq2
Set Up Files for DeSeq2:
1. In Python:
  - extract desired columns
  - write data in each quant file to a more useful format
  - convert all counts to int data type
  - increment all counts by 1 to get rid of any zeros
2. In R:
  - combine  modified quant files into one large table with each sample as a column
  - filter data as necessary/desired (i.e. average counts > 10)
3. Write `colData` table to import to R
```
Sample  Condition Type
  #      Control   #
  ...
  #      Exp       #
  ...
```

### Plots

#### Principle Component Analysis
Determines which component contributes the most variance in the samples used
> i.e. pc1 contributes the most variance, pc2 contributes the next second most amount of variance, etc.   

There are other components than the two shown but the variance contributed from these is less abundant. The changes are likely pretty subtle and noise can contribute a lot to the clustering.  
To minimize skew from noise, look at control data metrics and structure the dataset/analysis

```R
rld <- rlog(dds, blind = TRUE)
plotPCA(rld, intgroup = c("condition", "Type"))   # columns of colData table
```

#### Heat Map (Component Analysis)
Shows similarities between the samples (measure of variance).  
Adjacent groups (especially those connected by hierarchy) are very close.  
Vertical distance is also proportional to actual 'distance' between samples.

```R
rld_sampledist <- dist(t(assay(rld)))

library("RColorBrewer")
library(pheatmap)

rld_sampledistmatrix <- as.matrix(rld_sampledist)
rownames(rld_sampledistmatrix) <-
  paste(rld$Condition1, rld$Type, sep = "-")
colnames(rld_sampledistmatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

pheatmap(
  rld_sampledistmatrix,
  clustering_distance_rows = rld_sampledist,
  clustering_distance_cols = rld_sampledist,
  col = colors
)
```

### Normalized Counts
The normalized counts data can be used for further analysis with additional heat maps

```R
dds <- estimateSizeFactors(dds)
counts <- counts(dds, normalized = TRUE)
write.csv(counts, file = "norm_counts.csv")
```

#### MA Plots
Look at data along the x-axis, where the data on the right are highly expressed.  
Red dots indicate significantly differentially expressed genes.

```R
dds <- DESeq(dds)
plotMA(dds, ylim = c(-10,10))
```

#### Heat Maps
> 2 heatmaps per analysis: (1) upregulated genes (2) downregulated genes

1. Standardize normalized data counts by reads per million
  - sum each column
  - divide each individual read by the column total
  - multiply by one million
  - log transform (log2 or log10) as necessary
2. Set cutoffs as necessary. For example,
  - abs(log2foldchange) > 1
  - basemean > 100
  - padj < 0.01)
3. Plot with gene symbol
4. Use `pheatmap` library in R

```R
pheatmap(up_norm,
  scale = "none",
  cellwidth = 20,
  cellheight = 15,
  cluster_rows = TRUE,
  cluster_cols = FALSE)
```

### Pairwise Analysis
> Post-DESeq modifications: Match ensembl id's to gene names and other details

#### Volcano Plots
`EnhancedVolcano` library in R : simple and quick way to construct volcano plots

> Restrict data to include only genes with baseMean > 100

```
EnhancedVolcano(data,           # name of dataset analyzed
  lab = rownames(data),         # labels for data (i.e. gene names)
  x = 'log2FoldChange',         # data column to be used as x-axis
  y = 'padj',                   # data column to be used as y-axis
  xlim = c(-8, 8),
  ylim = c(-1,22),
  title = 'MY_TITLE',
  pCutoff = 0.01,               # p-value cutoff (y-axis)
  FCcutoff = 1,                 # fold change cutoff (x-axis)
  pointSize = 3.0,
  labSize = 0,
  legendPosition = 'right',
  legendLabSize = 12,
  legendIconSize = 4.0)
```

#### Gene Set Enrichment Analysis
This analysis takes a set of gene expression data and systematically looks for signatures of genes (up/down)-regulated in certain conditions.  
FDR: false discovery rate (comparable to padj because it takes into account extra factors)  

1. Download the program from [the GSEA website](https://www.gsea-msigdb.org/gsea/downloads.jsp)
>The site will require you to register with an email and an organization name

2. Make a `.cls file`
```
[numSamples] [numGroups] 1
# CTR EXP
ctr ctr ctr ctr exp exp exp exp
```

3. Modify the normalized count values file  
  - Only include protein_coding genes  
  - Subset by average reads > 100  
  - Include gene name and description columns, with `ctr_{#}` and `exp_{#}` as column names for the sample counts  
  - Write table as a `.txt` file  
  - Open file in Excel and save
    > Not sure why this works but it does. Otherwise, GSEA complains of about the file type and won't upload it to the software)  

4. Run GSEA with correct parameters  
  - **Expression Dataset**: your modified normalized count file
  - **Gene Sets Database**: Choose based on desired analysis  
  - **Number of permutations**: default (1000)  
  - **Phenotype labels**: CTR_versus_EXP
  - **Collapse/Remap**: default
  - **Permutation type**: gene_set
  - **Chip Platform**: Mouse_Gene_Symbol_Remapping_MSigDB.v7.0.chip

5. Analyze Results  
  - Sort `gsea_report_for_CTR/EXP_#####.xls` files by FDR and filter for genes with FDR < 0.01
  - Look at 'ranked_gene_list.####.xls' file for combined information
