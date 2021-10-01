# RNASeq-Pipeline
Simplified Bioinformatics Pipeline for Scientists Analyzing RNASeq Data

## Table of Contents
- [Setting Up Your Work Environment](#environment)
- [Working in a Cloud Environment](#cloud-computing)
  - [Vanderbilt ACCRE Cluster Cloud Computing](#vanderbilt-accre)
- [Transcript Quantification with Salmon](#salmon)
- [Quality Checking with FastQC Analysis](#fastqc)
- [Differential Gene Expression Analysis](#differential-gene-expression-analysis)
  - [DESeq2](#deseq2)
  - [Plots](#plots)
    - [Principle Component Analysis (PCA)](#principle-component-analysis)
    - [Heat Maps](#heat-maps)
  - [Normalized Counts Data](#normalized-counts)
    - [MA Plots](#ma-plots)
    - [Heat Maps](#heat-maps)
  - [Pairwise Analysis](#pairwise-analysis)
    - [Volcano Plots](#volcano-plots)
  - [Gene Set Enrichment Analysis (GSEA)](gene-set-enrichment-analysis)
- 

## Environment
> How to set up the basic work environment for the RNASeq Pipeline,

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

```bash
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

cat fileName                                                # Print the contents of the file to terminal
head fileName                                               # View the first 10 lines (change the number of lines using the `-n #` flag)
grep something file > newFile                               # Search for occurrences of 'something' in the file and store the output in a new file


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

### Vanderbilt ACCRE
> Vanderbilt University's Cluster Computing Resource - ACCRE

## Salmon

## FastQC

## Differential Gene Expression Analysis

### DeSeq2

### Plots

#### Principle Component Analysis

#### Heat Maps

### Normalized Counts

#### MA Plots

#### Heat Maps

### Pairwise Analysis

#### Volcano Plots

#### Gene Set Enrichment Analysis
