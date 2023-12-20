# ROAST 3.0 MultiPrior Upgrade

To use MultiPrior, call `roast(subj, 'multiprior', 'on')`. This will not work if a previous segmentation was done before the upgrade because the options and option files have changed.

`'multiprior'` -- advanced options of ROAST, for controlling multiprior segmentation
`'on' | 'off' (default)`  
The most recent ROAST update will have a new feature installed: MultiPrior. The purpose of MultiPrior is to provide a more accurate MRI segmentation when dealing with abnormalities such as a lesion or an enlarged brain stem that the original ROAST does not consider while segmenting. The MultiPrior AI associates tissue type with neighboring tissue to make predictions which will successfully segment important tissue types as accurately as possible. The goal is to simply add multiprior as simply as possible as it runs using Python. Scripts will be added to roast-3.0 which allow MultiPrior to run seamlessly. A conda enviornment will have to be created using the system specific .yml file.
(see [Example](#example)).  

#### Example 

    roast('example/subject1.nii',[],'multiprior','on')

Run simulation on subject1 with default recipe, but use multiprior segmentations instead of spm.

## Installing MultiPrior

To add MultiPrior to ROAST you will need to download the MultiPrior zip file and make a few changes to the `roast.m` file and other files.
The MultiPrior zip file contains the following:

- `MultiPrior.m`: Move this to the `roast-3.0` folder.
- `renameFiles.M`: Move this to the `roast-3.0` folder.

- `Defs.m`: Move this to the `spm12` folder in the `lib` folder in `roast-3.0`.

- `MultiPriors_WEB` folder, which includes the following files:

  - `WARP_indiTPM.m`
  - `SEGMENT.m`
  - `SEGMENT.py`
  - `Segmentation_config.py`
  - `best_model_tf2_singleGPUepoch100.h5`
  - `scripts` folder 
  -  Create MultiPriors Enviornment using the batch file specific to your operating system:

```
WINDOWS

Double click and the 'setupWindows.bat' file and it will automatically download the enviornment.
It will use the 'MultiPriors_env.yml' file 
*Anconda must already be installed for this to work.
```

```
LINUX (Ubuntu 22.04)

Run the 'setupLinux.sh' file in bash and it will download the enviornment.
It will use the 'MultiPriors_env_Linux.yml' file 
*Anconda must already be installed for this to work.
```

```
macOS

Run the 'setupMac.sh' file in bash and it will download the enviornment.
It will use the 'MultiPriors_env_Mac.yml' file 
*Anconda must already be installed for this to work.
```

```
Alternative

Create your own conda enviornment named 'MultiPriors_env' with the proper OS tag (None, Linux, or Mac) 
Place the enviornment into your MultiPriors_WEB Folder
Conda Install:
-tensorflow
-scikit-image
-nibabel
```
