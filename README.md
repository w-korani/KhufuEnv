# KhufuEnv

- single-thread process.
- does not need big compuataional resources.
- does a simple job.

## Requirments
R
bedtools

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
   sudo ./installer.sh
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
   sudo ./uninstaller.sh
```
3. remove the source for the Bash Shell Environment
```
   sed -i "/^source \/etc\/KhufuEnv\/call.sh$/d"  ~/.bashrc
```
4. refresh the Bash Shell Environment
```
   . ~/.bashrc
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
![Image](https://github.com/user-attachments/assets/b0c86fbb-970e-4606-912c-57444a0c41f7)

![Image](https://github.com/user-attachments/assets/384de179-8138-45b5-85cd-42099890530b)

## Citation


