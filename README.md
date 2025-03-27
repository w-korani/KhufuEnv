<div align="center">
  <center><h1>KhufuEnv</h1></center>
</div>

![Image](https://github.com/user-attachments/assets/ba498cd7-e1e3-404a-a326-a62ebb6fa23d)

KhufuEnv, an open-source, flexible, auxiliary environment for manipulating and analyzing genomic data, among other datasets, in the Unix environment. The KhufuEnv provides a buildable platform for constructing pipelines for several genomic analyses across different species. 

Every tool in KhufuEnv should meet the following criteria
- single-thread process.
- does not need big computational resources.
- does a simple job.
- to be run using bash script
- does not need special dependencies or other third-party computational tools

## Requirements
R package (basic, in addition to data.table, tidyr & plyr packages)

## Installation

1. download the package
   ```
   git clone https://github.com/w-korani/KhufuEnv
   ```
2. go to the package folder
   ```
   cd KhufuEnv_main
   ```
3. run the installer
   ```
   sudo bash ./installer.sh
   ```
4. add the source for the Bash Shell Environment
   ```
   echo "source /etc/KhufuEnv/call.sh"  >>  ~/.bashrc
   ```
5. refresh the Bash Shell Environment
   ```
   . ~/.bashrc
   ```


## Uninstallation
1. go to the package folder
   ```
   cd KhufuEnv_main
   ```
2. run the uninstaller
   ```
   sudo bash ./uninstaller.sh
   ```
3. remove the source for the Bash Shell Environment
   ```
   sed -i "/^source \/etc\/KhufuEnv\/call.sh$/d"  ~/.bashrc
   ```
4. refresh the Bash Shell Environment
   ```
   . ~/.bashrc  
   ```
   or
   ```
   exec bash
   ```


## Getting Help
- to list all tools
```
   KhufuEnvHelp
```
- to show the documentation of a specific tool
```
   KhufuEnvHelp tool_name
```


## Citation
bioxiv

## Available tools
![Image](https://github.com/user-attachments/assets/a52acdee-73b0-4a68-87ee-69dc2a1b65c9)

