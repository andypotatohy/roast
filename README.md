# ROAST: Realistic vOlumetric-Approach-based Simulator for Transcranial electric stimulation

1. [Getting started](#getting-started)
2. [How to use `roast`](#how-to-use-roast)
   - [Synopsis of `roast`](#synopsis-of-roast)
   - [Examples on `roast`](#examples-on-roast)
3. [How to use `roast_target`](#how-to-use-roast_target)
   - [Synopsis of `roast_target`](#synopsis-of-roast_target)
   - [Examples on `roast_target`](#examples-on-roast_target)
4. [More notes on the `capInfo.xls` file](#more-notes-on-the-capInfoxls-file)
5. [Outputs of ROAST software](#outputs-of-roast-software)
   - [Outputs of `roast`](#outputs-of-roast)
   - [Outputs of `roast_target`](#outputs-of-roast_target)
6. [Review of simulation data](#review-of-simulation-data)
7. [How to ask questions](#how-to-ask-questions)
8. [Acknowledgements](#acknowledgements)
9. [Notes](#notes)
10. [License](#license)

## Getting started

After you download the zip file, unzip it, launch your Matlab, make sure you are under the root directory of ROAST (i.e., you can see `example/`, `lib/`, and all other files), and then enter:

`roast`

This will demo a modeling process on the MNI152 head. Specifically, it will use the T1 image of the [6th gen MNI-152 head](http://nist.mni.mcgill.ca/?p=858) to build a TES model with anode on Fp1 (1 mA) and cathode on P4 (-1 mA).

There are 3 main functions that you can call: `roast()`, `roast_target()` and `reviewRes()`. The following sections will cover how to use them.

## How to use `roast`

### Synopsis of `roast`

`roast(subj,recipe,varargin)`

`subj`: file name (can be either relative or absolute path) of the MRI of
the subject that you want to run simulation on. The MRI can be either T1
or T2. If you have both T1 and T2, then put T2 file in the option `'T2'`
(see below options for details). If you do not have any MRI but just want
to run ROAST for a general result, you can use the default subject the
[MNI152 averaged head](http://nist.mni.mcgill.ca/?p=858) (see [Example 1](#example-1)) or the [New York head](https://www.parralab.org/nyhead/) (see [Example 2](#example-2)).

`recipe`: how you want to ROAST the subject you specified above. Default
recipe is anode on Fp1 (1 mA) and cathode on P4 (-1 mA). You can specify
any recipe you want in the format of `electrodeName-injectedCurrent` pair
(see [Example 3](#example-3)). You can pick any electrode from the 10/20, 10/10, 10/05 or BioSemi-256
EEG system (see the Microsoft Excel file `capInfo.xls` under the root directory
of ROAST). The unit of the injected current is in milliampere (mA). Make sure
they sum up to 0. You can also place electrodes at customized locations
on the scalp. See [Example 5](#example-5) for details.

`varargin`: Options for ROAST can be entered as `Name-Value` pairs in the 3rd argument 
(available from ROAST v2.0). The syntax follows the Matlab convention (see `plot()` for example).

*If you do not want to read the detailed info on the options below, you can go to [Example 23](#example-23) for quick reference.*

`'capType'` -- the EEG system that you want to pick any electrode from.  
`'1020' | '1010' (default) | '1005' | 'BioSemi'`  
You can also use customized electrode locations you defined. Just provide
the text file that contains the electrode coordinates. See below [Example 5](#example-5) for details.

`'elecType'` -- the shape of electrode.  
`'disc' (default) | 'pad' | 'ring'`  
Note you can specify different shapes to different electrodes, i.e., you
can place different types of electrodes at the same time. See below
[Example 6](#example-6) for details.

`'elecSize'` -- the size of electrode.  
All sizes are in the unit of millimeter (mm). For disc electrodes, sizes
follow the format of `[radius height]`, and default size is `[6mm 2mm]`; 
for pad electrodes, sizes follow the format of `[length width height]`, and
 default size is `[50mm 30mm 3mm]`; for ring electrodes, sizes follow the 
format of `[innerRadius outterRadius height]`, and default size is `[4mm 6mm 2mm]`.

If you're placing only one type of electrode (e.g., either disc, or pad,
or ring), you can use a *one-row* vector to customize the size, see [Example
7](#example-7); if you want to control the size for each electrode separately
(provided you're placing only one type of electrode), you need to specify
the size for each electrode correspondingly in a *N-row* matrix, where *N* is
the number of electrodes to be placed, see [Example 8](#example-8); if you're placing
more than one type of electrodes and also want to customize the sizes,
you need to put the size of each electrode in a *1-by-N* cell (put `[]` for
any electrode that you want to use the default size), where *N* is the number
of electrodes to be placed, see [Example 9](#example-9).

`'elecOri'` -- the orientation of pad electrode (ONLY applies to pad electrodes).  
`'lr' (default) | 'ap' | 'si' | direction vector of the long axis`  
For pad electrodes, you can define their orientation by giving the
direction of the long axis. You can simply use the three pre-defined keywords:  
- `'lr'`--long axis going left (l) and right (r);
- `'ap'`--long axis pointing front (anterior) and back (posterior);
- `'si'`--long axis going up (superior) and down (inferior).  
For other orientations you can also specify the direction precisely by giving the direction vector of the long axis.

If you're placing pad electrodes only, use the pre-defined keywords
([Example 10](#example-10)) or the direction vector of the long axis ([Example 11](#example-11)) to
customize the orientations; if you want to control the
orientation for each pad electrode separately, you need to specify
the orientation for each pad correspondingly using the pre-defined
keywords in a *1-by-N* cell ([Example 12](#example-12)) or the direction vectors of the
long axis in a *N-by-3* matrix ([Example 13](#example-13)), where *N* is the number of pad 
electrodes to be placed; if you're placing more than one type of electrodes
and also want to customize the pad orientations, you need to put the
orientations into a *N-by-3* matrix ([Example 14](#example-14); or just a *1-by-3* vector or a
single pre-defined keyword if same orientation for all the pads) where *N*
is the number of pad electrodes, or into a *1-by-N* cell ([Example 15](#example-15)), where *N*
is the number of all electrodes to be placed (put `[]` for non-pad electrodes).

`'T2'` -- use a T2-weighted MRI to help segmentation.  
`[] (default) | file path to the T2 MRI`  
If you have a T2 MRI aside of T1, you can put the T2 file in this option,
see [Example 16](#example-16), note you should put the T1 and T2 files in the same
directory.
If you ONLY have a T2 MRI, put the T2 file in the first argument `'subj'`
when you call roast, just like what you would do when you only have a T1.

`'meshOptions'` -- advanced options of ROAST, for controlling mesh parameters
(see [Example 17](#example-17)).  
5 sub-options are available:
- `meshOpt.radbound`: maximal surface element size, default 5;
- `meshOpt.angbound`: mimimal angle of a surface triangle, default 30;
- `meshOpt.distbound`: maximal distance between the center of the surface bounding circle
and center of the element bounding sphere, default 0.3;
- `meshOpt.reratio`: maximal radius-edge ratio, default 3;
- `meshOpt.maxvol`: target maximal tetrahedral element volume, default 10.  
See iso2mesh documentation for more details on these options.

`'simulationTag'` -- a unique tag that identifies each simulation.  
`dateTime string (default) | user-provided string`  
This tag is used by ROAST for managing simulation data. ROAST can
identify if a certain simulation has been already run. If yes, it will
just load the results to save time. You can leave this option empty so 
that ROAST will just use the date and time as the unique tag for the
simulation. Or you can provide your preferred tag for a specific
simulation ([Example 18](#example-18)), then you can find it more easily later. Also all the
simulation history with options info for each simulation are saved in the
log file (named as `"subjName_roastLog"`), parsed by the simulation tags.

`'resampling'` -- re-sample the input MRI to 1mm isotropic resolution.  
`'on' | 'off' (default)`  
Sometimes the input MRI has a resolution of not being 1 mm, e.g., 0.6 mm.
While higher resolution can give more accurate models, the computation
will be more expensive and thus slower. If you want a faster simulation,
you can ask ROAST to resample the MRI into 1 mm isotropic resolution by
turning on this option ([Example 19](#example-19)). Also it is recommended to turn on
this option if your input MRI has anisotropic resolution (e.g., 1 mm by
1.2 mm by 1.2 mm), as the electrode size will not be exact if the model
is built from an MRI with anisotropic resolution.

`'zeroPadding'` -- extend the input MRI by some amount, to avoid
complications when electrodes are placed by the image boundaries. Default
is not padding any slices to the MRI. You can ask ROAST to pad *N* empty
slices in all the six directions of the input MRI (left, right, front,
back, up and down), where *N* is a positive integer. This is very useful
when placing big electrodes on the locations close to image boundaries
([Example 20](#example-20)). This is also useful for MRIs that are cut off at the nose. 
If you specify a zeropadding of, say, 60 slices, ROAST can automatically get the segmentation 
of the lower part of the head, see ([Example 20](#example-20)).

`'conductivities'` -- advanced options of ROAST, the values are stored as a 
structure, with the following field names:
- `white` (default 0.126 S/m);
- `gray` (default 0.276 S/m);
- `csf` (default 1.65 S/m);
- `bone` (default 0.01 S/m);
- `skin` (default 0.465 S/m);
- `air` (default 2.5e-14 S/m);
- `gel` (default 0.3 S/m);
- `electrode` (default 5.9e7 S/m).  
You can use this option to customize the electrical conductivity for each tissue, each electrode, as well as the
conducting medium under each electrode. You can even assign different conductivity
values to different electrodes and their conducting media (e.g., `'gel'`). See
[Example 21](#example-21) and [Example 22](#example-22) for details.

### Examples on `roast`

#### Example 1

    roast

Default call of ROAST, will demo a modeling process on the MNI152 head.
Specifically, this will use the MRI of the MNI152 head to build a model
of transcranial electric stimulation (TES) with anode on Fp1 (1 mA) and cathode
on P4 (-1 mA). Electrodes are modeled by default as small disc electrodes.
See options below for details.

#### Example 2

    roast('nyhead')

ROAST New York head. Again this will run a simulation with anode on Fp1 (1 mA)
and cathode on P4 (-1 mA), but on the 0.5-mm resolution New York head. A decent
machine of 50GB memory and above is recommended for running New York
head. Again electrodes are modeled by default as small disc electrodes.
See options below for details.

#### Example 3

    roast('example/subject1.nii',{'F1',0.3,'P2',0.7,'C5',-0.6,'O2',-0.4})

Build the TES model on any subject with your own "recipe". Here we inject
0.3 mA at electrode F1, 0.7 mA at P2, and we ask 0.6 mA coming out of C5,
and 0.4 mA flowing out of O2. You can define any stimulation montage you want
in the 2nd argument, with `electrodeName-injectedCurrent` pair. Electrodes are
modeled by default as small disc electrodes. You can pick any electrode
from the 10/20, 10/10, 10/05 or BioSemi-256 EEG system. You can find all
the info on electrodes (names, locations, coordinates) in the Microsoft
Excel file `capInfo.xls` under the root directory of ROAST. Note the unit of
the injected current is milliampere (mA). Make sure they sum up to 0.

#### Example 4

    roast('example/subject1.nii',{'G12',1,'J7',-1},'captype','biosemi')

Run simulation on subject1 with anode on G12 (1 mA) and cathode on J7 (-1
mA) from the extended BioSemi-256 system (see `capInfo.xls` under the root
directory of ROAST).
 
#### Example 5

    roast('example/subject1.nii',{'G12',0.25,'J7',-0.25,'Nk1',0.5,'Nk3',-0.5,'custom1',0.25,'custom3',-0.25},'captype','biosemi')

Run simulation on subject1 with recipe that includes: BioSemi electrodes
G12 and J7; neck electrodes Nk1 and Nk3 (see `capInfo.xls`); and
user-provided electrodes custom1 and custom3. You can use a free program
called [MRIcro](http://www.mccauslandcenter.sc.edu/crnl/mricro) to load
the MRI first (note do NOT use MRIcron for this as MRIcron will not give you
the true voxel coordinates) and click the locations on the scalp surface where you want
to place the electrodes, record the voxel coordinates returned by MRIcro
into a text file, and save the text file to the MRI data directory with name
`"subjName_customLocations"` (e.g., here for subject1 it's saved as
`"subject1_customLocations"`). ROAST will load the text file and place the
electrodes you specified. You need to name each customized electrode in
the text file starting with `"custom"` (e.g., for this example they're
named as `custom1`, `custom2`, etc. You can of course do
`"custom_MyPreferredElectrodeName"`).

#### Example 6

    roast([],{'Fp1',1,'FC4',1,'POz',-2},'electype',{'disc','pad','ring'})

Run simulation on the MNI152 averaged head with the specified recipe. A disc
electrode will be placed at location Fp1, a pad electrode will be placed at FC4,
and a ring electrode will be placed at POz. The sizes and orientations will be set
as default.

#### Example 7

    roast('nyhead',[],'electype','ring','elecsize',[7 10 3])

Run simulation on the New York head with default recipe. Ring electrodes
will be placed at Fp1 and P4. The size of each ring is 7mm inner radius,
10mm outter radius and 3mm height.

#### Example 8

    roast('nyhead',{'Fp1',1,'FC4',1,'POz',-2},'electype','ring','elecsize',[7 10 3;6 8 3;4 6 2])

Run simulation on the New York head with the specified recipe. Ring electrode
placed at Fp1 will have size `[7mm 10mm 3mm]`; ring at FC4 will have size
`[6mm 8mm 3mm]`; and ring at POz will have size `[4mm 6mm 2mm]`.

#### Example 9

    roast([],{'Fp1',1,'FC4',1,'POz',-2},'electype',{'disc','pad','ring'},'elecsize',{[8 2],[45 25 4],[5 8 2]})

Run simulation on the MNI152 averaged head with the specified recipe. A disc
electrode will be placed at location Fp1 with size `[8mm 2mm]`, a pad electrode
will be placed at FC4 with size `[45mm 25mm 4mm]`, and a ring electrode will be
placed at POz with size `[5mm 8mm 2mm]`.

#### Example 10

    roast([],[],'electype','pad','elecori','ap')

Run simulation on the MNI152 averaged head with default recipe. Pad
electrodes will be placed at Fp1 and P4, with default size of `[50mm 30mm
3mm]` and the long axis will be oriented in the direction of front to back.

#### Example 11

    roast([],[],'electype','pad','elecori',[0.71 0.71 0])

Run simulation on the MNI152 averaged head with default recipe. Pad
electrodes will be placed at Fp1 and P4, with default size of `[50mm 30mm
3mm]` and the long axis will be oriented in the direction specified by the
vector `[0.71 0.71 0]`.

#### Example 12

    roast('example/subject1.nii',{'Fp1',1,'FC4',1,'POz',-2},'electype','pad','elecori',{'ap','lr','si'})

Run simulation on subject1 with specified recipe. Pad electrodes will be 
placed at Fp1, FC4 and POz, with default size of `[50mm 30mm 3mm]`. The long
axis will be oriented in the direction of front to back for the 1st pad,
left to right for the 2nd pad, and up to down for the 3rd pad.

#### Example 13

    roast('example/subject1.nii',{'Fp1',1,'FC4',1,'POz',-2},'electype','pad','elecori',[0.71 0.71 0;-0.71 0.71 0;0 0.71 0.71])

Run simulation on subject1 with specified recipe. Pad electrodes will be 
placed at Fp1, FC4 and POz, with default size of `[50mm 30mm 3mm]`. The long
axis will be oriented in the direction of `[0.71 0.71 0]` for the 1st pad,
`[-0.71 0.71 0]` for the 2nd pad, and `[0 0.71 0.71]` for the 3rd pad.

#### Example 14

    roast([],{'Fp1',1,'FC4',1,'POz',-2},'electype',{'pad','disc','pad'},'elecori',[0.71 0.71 0;0 0.71 0.71])

Run simulation on the MNI152 averaged head with specified recipe. A disc
electrode will be placed at FC4. Two pad electrodes will be placed at Fp1
and POz, with long axis oriented in the direction of `[0.71 0.71 0]` and 
`[0 0.71 0.71]`, respectively.

#### Example 15

    roast([],{'Fp1',1,'FC4',1,'POz',-2},'electype',{'pad','disc','pad'},'elecori',{'ap',[],[0 0.71 0.71]})

Run simulation on the MNI152 averaged head with specified recipe. A disc
electrode will be placed at FC4. Two pad electrodes will be placed at Fp1
and POz, with long axis oriented in the direction of front-back and `[0 0.71 0.71]`, respectively.

#### Example 16

    roast('example/subject1.nii',[],'T2','example/subject1_T2.nii')

Run simulation on subject1 with default recipe. The T2 image will be used
for segmentation as well.

#### Example 17

    roast([],[],'meshoptions',struct('radbound',4,'maxvol',8))

Run simulation on the MNI152 averaged head with default recipe. Two of
the mesh options are customized.

#### Example 18

    roast([],[],'simulationTag','roastDemo')

Give the default run of ROAST a tag as `'roastDemo'`.

#### Example 19

    roast('example/subject1.nii',[],'resampling','on')

Run simulaiton on subject1 with default recipe, but resample the MRI of
subject1 to 1mm isotropic resolution first (the original MRI of subject1
has resolution of 1mm by 0.99mm by 0.99mm).

#### Example 20

    roast([],{'Exx19',1,'C4',-1},'zeropadding',60,'simulationTag','paddingExample')

Run simulation on the MNI152 averaged head, but add 60 empty slices on
each of the six directions to the MRI first, to allow placement of
electrode Exx19, which is outside of the MRI (i.e., several centimeters
below the most bottom slice of the MRI). This zeropadding also will generate the 
segmentation of the lower part of the head, thanks to the extended TPM coming along 
with ROAST. You can visually check this by 

    reviewRes([],'paddingExample','all')

If you run this without zero-padding first, you'll get strange results.
Note it is always a good practice to add empty slices to the MRI if you 
want to place electrodes close to, or even out of, the image boundary.
ROAST can detect if part or all of your electrode goes out of image boundary,
but sometimes it cannot tell (it's not that smart yet :-). So do a `'zeroPadding'`
of 10 to start with, and if you're not happy with the results, just increase
the amount of zero-padding. But the best solution is to get an MRI that covers
the area where you want to place the electrodes.

#### Example 21

    roast([],{'Fp1',1,'FC4',1,'POz',-2},'conductivities',struct('csf',0.6,'electrode',0.1))

Run simulation on the MNI152 averaged head with specified recipe. The
conductivity values of CSF and electrodes are customized. Conductivities
of other tissues will use the literature values.

#### Example 22

    roast([],{'Fp1',1,'FC4',1,'POz',-2},'electype',{'pad','disc','pad'},'conductivities',struct('gel',[1 0.3 1],'electrode',[0.1 5.9e7 0.1]))

Run simulation on the MNI152 averaged head with specified recipe.
Different conductivities are assigned to pad and disc electrodes. For pad
electrodes, `'gel'` is given 1 S/m and `'electrode'` is 0.1 S/m; for the disc
electrode, `'gel'` is given 0.3 S/m and `'electrode'` is 59000000 S/m. When
you control the conductivity values for each electrode individually, keep
in mind that the values you put in the vector in `'gel'` and `'electrode'`
field in `'conductivities'` option should follow the order of electrodes
you put in the `'recipe'` argument.

#### Example 23

All the options can be combined to meet your specific simulation needs.

    roast('path/to/your/subject.nii',{'Fp1',0.3,'F8',0.2,'POz',-0.4,'Nk1',0.5,'custom1',-0.6},...
        'electype',{'disc','ring','pad','ring','pad'},...
        'elecsize',{[],[7 9 3],[40 20 4],[],[]},...
        'elecori','ap','T2','path/to/your/t2.nii',...
        'meshoptions',struct('radbound',4,'maxvol',8),...
        'conductivities',struct('csf',0.6,'skin',1.0),...
        'resampling','on','zeropadding',30,...
        'simulationTag','awesomeSimulation')

Now you should know what this will do.

## How to use `roast_target`

### Synopsis of `roast_target`

`roast_target(subj,simTag,targetCoord,varargin)`

TO BE ADDED...

### Examples on `roast_target`

TO BE ADDED...

## More notes on the `capInfo.xls` file

A lot of info are hidden in the fancy `capInfo.xls` file under the ROAST root directory. There you can find the comprehensive layouts of both [10-05](https://www.sciencedirect.com/science/article/pii/S1053811906009724?via%3Dihub#fig6) and [BioSemi](https://www.biosemi.com/pics/cap_256_layout_medium.jpg) systems (with my personal drawings), and also visually-striking 3D renderings of the New York head with these electrodes placed on. So make sure to check it out.


## Outputs of ROAST software

### Outputs of `roast`

#### Figure outputs

ROAST outputs 7 or 8 figures for quick visualization of the simulation
results. These figures include the slice view of the MRI (T1 and/or T2) and the segmentation; 3D rendering of the computed voltage and electric field distribution; and the slice view of the voltage and electric field. Note the slice view is always in the MRI voxel space, and the 3D rendering displays the data in the world space.

#### Outputs in Matlab format

ROAST saves the results as `"subjName_simulationTag_roastResult.mat"`, where 3 variables are available:

`vol_all`: the voltage at each pixel in the MRI voxel space, unit in mV.

`ef_all`: the electric field vector at each pixel in the MRI voxel space, unit in V/m. This variable includes 3 volumes, representing the x-, y-, and z-component of the electric field.

`ef_mag`: the magnitude of the electric field at each pixel in the MRI voxel space, unit in V/m.

#### Outputs in [NIfTI](https://nifti.nimh.nih.gov/) format

Voltage: `"subjName_simulationTag_v.nii"`, unit in mV.

E-field: `"subjName_simulationTag_e.nii"`, unit in V/m.

E-field magnitude: `"subjName_simulationTag_emag.nii"`, unit in V/m.

#### Outputs in text files

Voltage: `"subjName_simulationTag_v.pos"`, unit in mV.

E-field: `"subjName_simulationTag_e.pos"`, unit in V/m.

Note in these text files, voltage and electric field are defined at each mesh node, whose location can be found in the mesh file `"subjName_simulationTag.msh"` or `"subjName_simulationTag.mat"`. Also note that in these two mesh files the node coordinates are in the voxel space but with the scaling factors in the MRI header applied, i.e., the unit of the mesh coordinates is millimeter (mm).

### Outputs of `roast_target`

TO BE ADDED...

## Review of simulation data

You can also use the other function `reviewRes()` to review/visualize the
simulations that you already run before. `reviewRes()` has a simpler
interface than `roast()` so that you do not have to enter all the
simulation parameters again as you would have to do in `roast()`. For example, 
to review the results generated from running [Example 23](#example-23), you can
simply type

    reviewRes('path/to/your/subject.nii','awesomeSimulation')

where `'awesomeSimulation'` is the corresponding `'simulationTag'`. For more info 
on how to use this function, type `help reviewRes`.

TO BE ADDED...

## How to ask questions

Please read the [Getting started](#getting-started) and the [Synopsis](#synopsis) sections first. It'll only cost you 10 minutes. If you're lazy or tired, just go to [Example 23](#example-23) for quick reference. Do not forget to check out the fancy [`capInfo.xls`](#more-notes-on-the-capInfoxls-file) file.

If ROAST crashes, first check if there are any warning messages output in the command window. If there are, follow the suggestions in the warning messages. This will usually fix the problems and let you run through to get your model.

If you are still stuck, please subscribe to the ROAST user mailing list by clicking [here](https://groups.google.com/group/roast-users/subscribe). You can ask any questions and share your experiences there.

## Acknowledgements

If you use ROAST in your research, please cite these:

Huang, Y., Datta, A., Bikson, M., Parra, L.C., [Realistic vOlumetric-Approach to Simulate Transcranial Electric Stimulation -- ROAST -- a fully automated open-source pipeline](https://iopscience.iop.org/article/10.1088/1741-2552/ab208d), Journal of Neural Engineering, Vol. 16, No. 5, 2019. (prefered reference)

Huang, Y., Datta, A., Bikson, M., Parra, L.C., [ROAST: an open-source,
fully-automated, Realistic vOlumetric-Approach-based Simulator for TES](https://www.parralab.org/publications/ROAST_EMBC_forFinalSubmission.pdf), Proceedings of the 40th Annual International Conference of the IEEE Engineering in Medicine and Biology Society, Honolulu, HI, July 2018

If you use New York head to run simulation, please also cite the following:
Huang, Y., Parra, L.C., Haufe, S.,2016. [The New York Head - A precise
standardized volume conductor model for EEG source localization and tES
targeting.](https://www.sciencedirect.com/science/article/pii/S1053811915011325) NeuroImage,140, 150-162

If you also use the targeting feature (`roast_target`), please cite these:

TO BE ADDED...

ROAST was supported by NIH through grants R01MH111896, R01MH111439, R01NS095123, R44NS092144, R41NS076123, and by [Soterix Medical Inc](https://soterixmedical.com/).

## Notes

Version 3.0 is incompatible with simulation data generated by Version 2.7 and earlier versions.
If you do not have Matlab, there is [a Docker version](https://hub.docker.com/r/amiklos/roast/).
ROAST was not designed to build models for pathological heads, but there are plans to add this capability in the future versions.

## License

General Public License version 3 or later. See LICENSE.md for details.

This software uses free packages from the Internet, except Matlab, which is a proprietary software by the MathWorks. You need a valid Matlab license to run this software.

ROAST is considered as an "aggregate" rather than "derived work", based on the definitions in [GPL FAQ](http://www.gnu.org/licenses/gpl-faq.html#MereAggregation). The ROAST license only applies to the scripts, documentations, and the individual MRI data under `example/` folder in this package and excludes those programs stored in the `lib/` directory. The software under `lib/` follow their respective licenses. This software is only intended for non-commercial use.

(c) [Yu (Andy) Huang](https://www.parralab.org/people/yu-andy-huang/), Parra Lab at CCNY

yhuang16@citymail.cuny.edu

August 2019
