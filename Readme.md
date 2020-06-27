# MICROPHONE ARRAY BEAMFORMING TOOLBOX

Thank you for downloading the Microphone array beamforming toolbox! 
This is a Matlab library of microphone array designs and beamformers developed while carrying out our research. 
It is therefore a tool that we share with researchers for them to use to advance their work, not necesarily a flawless Matlab implementation.
This file explains the basic principles of the toolbox. 
The different files are organised in subfolders depending on their purpose.

## License
MIT (see LICENSE file)

## Citation 
This toolbox was implemented to generate the microphone array simulations described in  
Blanco Galindo, M., Coleman, P., Jackson P. J. B. (2020). Microphone array geometries for horizontal spatial audio object capture with beamforming. *Journal of the Audio Engineering Society*, 68(5), 324--337.  

Where appropriate, users of the toolbox are encouraged to reference the above article.

## Contact
Miguel Blanco Galindo - mblancogalindo@gmail.com  
Philip Coleman, Institue of Sound Recording, University of Surrey - p.d.coleman@surrey.ac.uk  

## TestBeamfoming script
To begin with, open the file `TestBeamforming.m` in the `Scripts` subfolder.
This file will be the starting point for the high level implementation of the toolbox.

Make sure you add the folder and all subfolders in your Matlab path. `TestBeamforming.m` does this by default.
The structure of `TestBeamforming.m` is as follows:

 1. Specify the array(s) to be used.
 2. Set up the array simulation conditions and compute the array manifold.
 3. Set up the beamforming simulation conditions and compute the beamformer weights and evaluation metrics.
 4. Plot the beampattern/evaluation metrics for analysis. 

### Array design

Array geometry:  
Choose from: `linear_bs` (broadside), `linear_ef` (endfire), `circular`, `dualcircular`, `multicircular`, `square`, 
	`rectangular`, `circular_grid`, `circular_rigid_sphere`, `circular_rigid_cylinder`, `stacked_circular` 
	(vertically stacked circular), `stacked_circular_rigid_cylinder`, `spherical`, `spherical_rigid_sphere`.

Number of microphones  
Radius  
Potential simulation errors such as positioning offset, calibration error and sensor self noise can be included.   
To override all simulatiuon errors with a single attribute, set SensorOffet = [] (as per `TestBeamforming.m`).  

Design the array based on the settings calling `configArray`. The result is a cell array named `ConfigAllArrays` that has information about the above settings. 

Alternatively you can pass the microphone positions `MicPos` of an arbitrary array by uncommeting line 52 (and commenting line 50 out). 

### Set up and compute the array manifold

Feed aditional information to `ConfigAllArrays` that we just created.   
This information is mainly concerned with ArrayManifold calculation so it is included in a field named ArrayMan in this variable: 
 - LookDistance: choose from `farfield` (for plane wave propagation) or `focalpoint` (for point source propagation).
 - SteerSpace: choose from `2Daz`, `2Del` (horizontal and vertical cylindrical sound fields respectively) or `3D` (spherical sound field)
 - SteerGrid: for `3D` `SteerSpace` you can select the `SteerGrid`. 
	- `SteerGrid.Type` sets the sampling scheme. Choose from: `Uniform` (nearly uniform), `Lebedev`, `Equiangular`, `Gaussian`, `TruncatedIcosahedron`.
	- `SteerGrid.Points` set the number of sampling points. Depending on the `SteerGrid.Type`, `SteerGrid.Points` may be adjusted to find the closest possible values, or throw an error (e.g. Lebedev if it is not within its default set). 
 - AzResolution: for `2Daz` it is easier to specify the number of sampling points in azimuth. Set `ElResolution` to 1.
 - ElResolution: for `2Del` it is easier to specify the number of sampling points in elevation. Set `AzResolution` to 1.
 - PerformerDistance: distance of the sound source (when `LookDistance` is equal to `focalpoint` only)
 - PerformerAngle: elevation or azimuth angle plane when `SteerSpace` is `2Daz` or `2Del`, respectively, and `LookDistance` is equal to `focalpoint` only.
 
Filter is another field of `Config` regarding the frequency vector (and later on the beamformer), with the following fields:
 - Fs: sampling frequency.
 - FreqMode: choose from `discrete` for bespoke frequency vector or `filter` for FFT vector from 0 to Fs/2 uniformly sampled. 
 - Nfft: number of FFT points for `filter` `FreqMode`.
 - FreqVector: frequency vector.

Calculate the array manifold based on the settings above calling `GetArrayManifold`.

### Set up and compute the beamformer weights and evaluation metrics

Additonal settings for the `Filter` field in `Config`:
 - Method: is the beamforming method(s). Use cell array for compatibility with multiple methods. Choose from: `ds` (delay and sum), `sda` (superdirective array), `mvdr` (minimum variance distortionless response), `lcmv` (linearly contrained minimum variance), `acc` (acoustic contrast control) and `ls` (least squares).
 - ls: if `Method` is `ls`, additional parameters are required:
	- ls.mode: type of target pattern. Choose from `hypercardioid`, `cardioid` or `chebyshev`.
	- ls.order: order of the target pattern as integer.
	- ls.EQ: option to equalise the response at the look direction when inversion rolls off the beampattern. 
 - RegMethod: type of regularisation approach. Choose from: `external` (explicit regularisation parameter of diagonal loading) or `wnglimit` (lower bound robustness based on white noise gain (WNG) constraint).
 - RegParam: value of RegMethod. If `RegMethod` is equal to `external` it is the regularisation value of the diagonal loading, and if equal to `wnglimit` it is the WNG constraint in dB (where positive values imply increase in robustness at the beamformer output). 
 
The field Scene specifies the following subfields:
 - Targets: values corresponding to the target sources.
	- Targets.Gain: amplitude of the equivalent target source(s) (in natural units). 
	- Targets.Angle: angles of the target source(s). `L` x 2 where `L` is the number of look directions and 2 refers to elevation and azimuth (in that order).
 - Interferers: values corresponding to the interferer sources (if required for the type of beamformer such as MVDR and LCMV). 
	- Interferers.Gain: amplitude of the equivalent interferer source(s) (in natural units). 
	- Interferers.Angle: angles of the target source(s). `L` x (2*`I`) where I is the number of interferer directions. Thus, for each look direction there will be a row of `I` elevation and azimuth interfering angles.

Compute the beamformer weights and evaluate the beampattern and associated metrics calling `calculate_beamforming_weights_and_evaluation_different_arrays`.

### Plot the beampattern/evaluation metrics for analysis

Call `rearrange_directional_response_measures` to rearrange the data for plotting purposes. 

Specify parameters in `FigConfig` to display the arrays, beamformers and metrics of interest.

Call `PlotDirectionalResponse2DAllArrays` for "waterfall" plot of the beampattern as a function of angle and frequency.

Call `PlotDirectionalResponseMeasuresAllArrays` to plot the evaluation metrics of interest as a function of frequency.

Call `PlotFrequencyRange` to plot the operating frequency range of the arrays and beamformers of interest. 

Call `PlotLSBDirRespError` to plot the error between the achieved and target beampatterns for the least-squares beamformer. 

Call `PlotBFWeights` to plot the designed weights in frequency or time domain.


`TestBeamforming3D.m` is an equivalent script to the one described above, but performing 3D design and evaluation and using the `circular_rigid_cylinder` and `spherical_rigid_cylinder`.