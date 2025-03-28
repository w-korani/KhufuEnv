<div align="center">
  <center><h1>KhufuEnv</h1></center>
</div>

![Image](https://github.com/user-attachments/assets/ba498cd7-e1e3-404a-a326-a62ebb6fa23d)

KhufuEnv is an open-source, Unix-based auxiliary environment consisting of flexible tools for manipulating and analyzing genomic data and other datasets. The KhufuEnv provides a buildable platform for constructing pipelines for several genomic analyses across different species. 

Every tool in KhufuEnv is designed to meet the following criteria:
- Single-thread process.
- Does not require significant computational resources.
- Does a simple job. 
- To be run using bash script.
- Does not require special dependencies or other third-party computational tools.

<!-- You may use this if you like? :)
## Table of Contents
1. [Requirements](#requirements)
2. [Installation](#installation)
3. [Uninstallation](#uninstallation)
4. [Getting Help](#getting-help)
5. [Testing Data](#testing-data)
6. [Citation](#citation)
7. [Available Tools](#available-tools)
-->

## Requirements
### R Packages
- **Basic:** [R (latest version recommended)](https://www.r-project.org/)

- **Additional packages:**
   - data.table
   - tidyr
   - plyr

### Command-Line Tools
- gawk ([see the official GNU Awk website](https://www.gnu.org/software/gawk/manual/gawk.html#Installation))

## Installation

1. Download the package.
   ```
   git clone https://github.com/w-korani/KhufuEnv
   ```
2. Go to the package folder.
   ```
   cd KhufuEnv_main
   ```
3. Run the installer.
   ```
   sudo bash ./installer.sh
   ```
4. Add the source for the Bash Shell Environment.
   ```
   echo "source /etc/KhufuEnv/call.sh"  >>  ~/.bashrc
   ```
5. Refresh the Bash Shell Environment.
   ```
   . ~/.bashrc
   ```


## Uninstallation
1. Go to the package folder.
   ```
   cd KhufuEnv_main
   ```
2. Run the uninstaller.
   ```
   sudo bash ./uninstaller.sh
   ```
3. Remove the source for the Bash Shell Environment.
   ```
   sed -i "/^source \/etc\/KhufuEnv\/call.sh$/d"  ~/.bashrc
   ```
4. Refresh the Bash Shell Environment.
   ```
   . ~/.bashrc  
   ```
   or
   ```
   exec bash
   ```


## Getting Help
- To list all tools:
```
   KhufuEnvHelp
```
- To show the documentation of a specific tool:
```
   KhufuEnvHelp <tool-name>
```
## Testing Data
In order to help users test the functionality of these tools, a directory containing test input files for each section is provided. Each tool is associated with its own documentation, which includes a specific example of how to use it. 

### How to Test a Tool
- Navigate to the main test directory.
```
cd KhufuEnv_main/TestingData
```
- Each section has their own set of test inputs. You will need to navigate to the appropriate subdirectory for the section's tools you want to test, replacing `<section-inputs>` with the name of the section's corresponding directory.
```
cd <section-inputs> 
```
- Display the documentation for the tool using `KhufuEnvHelp <tool-name>`.

- In the tool's documentation, you will find an example command as shown below:
```
#####################
KhufuEnv.Ver1.0.0: hapmapFilterMissingVariant
Usage: hapmapFilterMissingVariant hapmap missing

Description:
 The hapmapFilterMissingVariant function filters the input hapmap to remove sites where the percentage of variants missing is greater than the specified threshold. 

Parameters:
 hapmap ex: test1.hapmap Hapmap to be processed. 
 missing ex: 0.75 Missing percentage to be used. Must be in decimal format.

Example:
 hapmapFilterMissingVariant test1.hapmap 0.75
#####################
 ``` 
- For instance, to test the output of the `hapmapFilterMissingVariant` tool, type or copy and paste `hapmapFilterMissingVariant test1.hapmap 0.75` into the terminal and execute it. 

- The output will be displayed as shown below:
```
$ hapmapFilterMissingVariant test1.hapmap 0.75 | column -t
chr         pos        SM001  SM002  SM003  SM004  SM005  SM006  SM007  SM008  SM009  SM010
TRv2Chr.01  296910     C      -      -      A,C    A,C    -      -      -      A,C    A,C
TRv2Chr.01  9470991    C      -      -      C      -      -      C      -      -      -
TRv2Chr.01  13207367   -      -      T      C      C      -      C      C,T    C      -
TRv2Chr.01  30479088   -      -      C      C      -      -      C      -      C      -
TRv2Chr.01  66933056   G      -      G      A,G    G      -      A,G    A,G    A,G    A,G
TRv2Chr.01  66944054   T      T      -      T      T      T      -      T      T      T
TRv2Chr.01  73021455   C      -      -      -      -      -      -      C      C      -
TRv2Chr.01  95590686   A,G    G      -      A,G    G      -      -      -      -      G
TRv2Chr.01  95593242   -      -      G      C      -      C      -      -      -      -
TRv2Chr.01  95608110   -      -      -      A      -      A,G    -      -      A      -
TRv2Chr.01  112412989  A,G    A,G    A      A,G    A      -      A,G    A,G    A,G    -
TRv2Chr.02  71100      C,T    C,T    C,T    C,T    -      -      C      C,T    C,T    -
TRv2Chr.02  74281      A      A      A      A      -      A      A      A      A      -
TRv2Chr.02  80175      C,T    C,T    C,T    C,T    C,T    C      C,T    C,T    C,T    -
TRv2Chr.02  107994     -      G      A      A      -      -      A,G    A,G    G      -
TRv2Chr.02  119526     -      G      G      G,T    -      -      G,T    G      G,T    -
TRv2Chr.02  27794987   -      -      -      G      -      -      -      G      A,G    -
TRv2Chr.02  67801286   C,T    C      T      C,T    C,T    -      -      -      C,T    C,T
TRv2Chr.02  84715892   -      A      -      A      -      -      -      A      -      -
TRv2Chr.03  31478009   T      T      -      C,T    T      T      C,T    T      C,T    -
TRv2Chr.03  106731401  -      A      -      A      -      -      -      -      A      -
```
**Note:** the command `column -t` has been piped for better display. It is not required.


## Citation
bioxiv

## Available tools
![Image](https://github.com/user-attachments/assets/a52acdee-73b0-4a68-87ee-69dc2a1b65c9)

