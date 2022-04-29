# EONS - Evaluation of Non-Sinusoidal Magnetic Fields
EONS is a MATLAB GUI for the evaluation of the neuronal response to sinusoidal and/or non-sinusoidal electromagnetic field exposure in the context of occupational safety and health. Conference papers on EONS will be presented at the BioEM 2022 Conference (https://www.bioem2022.org/) [1-2]. The neuronal response is calculated in EONS with the Spatially Extended Non-linear Node (SENN) model [3]. 

![fig_41251_1_git](https://user-images.githubusercontent.com/45427880/165744019-75f1aa61-7b41-491c-8032-819519cc0bbd.png)
**Fig. 1.** Screenshot of the EONS GUI: a magnetic flux density recording of a phase-cut sine is shown, recorded during resistance welding. 

The user can load a custom magnetic flux density waveform ('Load File') or specify a predefined electric field waveform ('Use Waveform'). Neuronal simulations are started with the 'SENN (Reilly)' or 'SENN (Ghent U.)' buttons, to run external Fortran code of Reilly & Diamant [4] (close the external window after the run is finished without saving the results) or Matlab code from Tarnaud et al. [5], respectively. The E-field is uniform and parallel to the axon, analogous to the setup used in [4] for deriving exposure guidelines.   

**1.** After a simulation is finished, exposure indices (EI) calculated with the Weighted Peak Method (WPM) and with the SENN excitation threshold are tabulated. The table also includes information on the excitation threshold level, the latency until the first action potential (AP Time) and the excitation node number (AP Node). The EI and WPM calculation is described in [6]. Note that the WPM filters change with the selected membrane model as described in [2].  \
**2.** For the SENN (Ghent U.) code, the used membrane model can be changed to Hodgkin-Huxley (HH), Frankenhaeuser-Huxley (FH), Chiu-Ritchie-Roggart-Stagg-Sweeney (CRRSS), Schwarz-Eikhof (SE) or Schwarz-Reid-Bostock (SRB) dynamics. For the SENN (Reilly) code the FH-dynamics is used. \
**3.** Solver options for the SENN (Ghent U.) code can be altered: 
- max. step and sampling/period define the temporal discretization (dt = min(max. step,(1/frequency)/(sampling/period))). 
- Order of solution (1 or 2) is the Crank-Nicolson order of accuracy. Order of solution equal to two will result in higher accuracy, but slower simulations. Note: for conditionally linear membrane dynamics (i.e., HH and CRRSS) the order of solution is always equal to two.
- Lower and upper limits can be used to define the lower and upper boundary of the excitation threshold to speed up the titration algorithm. If unknown, these limits should be set to zero. 
- Eapp: if zero, the excitation threshold is calculated. If non-zero, EONS will calculate the neuronal response for the applied electric field (Eapp [V/m]) and will display the corresponding membrane voltage plots.

**4.** A predefined waveform can be specified (sine or pulse) in the GUI.  

EONS was created using MATLAB R2021. In order to create an MLAPP file that the MATLAB AppDesigner can read, you have to put the contents of the src\ folder in a zip file and name it EONS.mlapp. After making changes to the MLAPP file you can unzip its contents back to the src\ folder and commit the changes to the repository.

The Logo.gif file and the EONS_parameters.txt must be in the MATLAB working directory when you run EONS from within the AppDesigner. 

The path for Reilly's SENN program has to be specified in the parameters file. It can be downloaded here:
https://us.artechhouse.com/Assets/reslib/reilly/reilly.html

A non-sinusoidal sample waveform is provided in the waveforms\ folder. It was recorded during welding using phase cutting.

# COPYRIGHT
Reilly's SENN nerve stimulation model was obtained from https://www.fda.gov/about-fda/cdrh-offices/senn-nerve-stimulation and is included in the repository.

# AUTHORS
Dr. Florian Soyka, Institute for Occupational Safety and Health of the German Social Accident Insurance, Florian.Soyka@dguv.de <br>
Dr. Thomas Tarnaud, INTEC-WAVES, Ghent University – IMEC, Thomas.Tarnaud@UGent.be

# REFERENCES
[1] F. Soyka, T. Tarnaud, W. Joseph, L. Martens, E. Tanghe, “The Influence of Membrane Channel Dynamics on Occupational Exposure Limit Values”, Joint Meeting of The Bioelectromagnetics Society and the European BioElectromagnetics Association BioEM2022 Nagoya, Japan, June 19-24 2022 \
[2] T. Tarnaud, F. Soyka, R. Schoeters, T. Plovie, W. Joseph, L. Martens, E. Tanghe, “EONS: Evaluation of Non-Sinusoidal Magnetic Fields for Electromagnetic Safety to Intermediate Frequencies”, Joint Meeting of The Bioelectromagnetics Society and the European BioElectromagnetics Association BioEM2022 Nagoya, Japan, June 19-24 2022 \
[3] J. P. Reilly, V. T. Freeman, W. D. Larkin (1985). Sensory effects of transient electrical stimulation-evaluation with a neuroelectric model. IEEE transactions on biomedical engineering, (12), 1001-1011. \
[4] J. P. Reilly, A. M. Diamant (2011) Electrostimulation: Theory, Applications, and Computational Models: Artech House Publishers. \
[5] T. Tarnaud, W. Joseph, L. Martens, E. Tanghe (2018). Dependence of excitability indices on membrane channel dynamics, myelin impedance, electrode location and stimulus waveforms in myelinated and unmyelinated fibre models. Medical & Biological Engineering & Computing, 56(9), 1595-1613. \
[6] F. Soyka, "Evaluation of non-sinusoidal magnetic fields: Comparing the weighted peak method with a new method using the Spatially Extended Nonlinear Node electrostimulation model", BioElectromagnetics Conference (BioEM 2021) Ghent, Belgium, September 26-30 2021 