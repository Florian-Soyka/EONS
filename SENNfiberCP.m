% Single neuron SENN solver -- Implemented SENN model: parametervalues are
% based on [McNeal1976], [Tasaki1955] and [FitzHugh1962]
function SENNfiberCP(Diamet,Ltot,x,y,z,Tp,ancat,monbi,Ibegin,Iend,SearchMode,SweepActI,ElecInj,model,conf,fileName,sinus_freq,varargin)
 LLx0s = str2double(x);LLy0s = str2double(y); LLz0s = str2double(z);
 LLTp = str2double(Tp); LLancat = str2double(ancat); LLmonbi = str2double(monbi);
 LLconf = str2double(conf); LLIbegin = str2double(Ibegin);LLIend = str2double(Iend);
 LLSearchMode=str2double(SearchMode);
 LLElecInj = str2double(ElecInj); LLSweepActI = str2double(SweepActI);
 sinus_freq = str2double(sinus_freq);
 if nargin >= 27
     error('Too many input arguments');
 elseif nargin == 18
 EONSWaveform = varargin{1};
 elseif nargin == 19
 EONSWaveform = varargin{1};
 SolverMaxStep = str2double(varargin{2});    
 elseif nargin == 20
 EONSWaveform = varargin{1};
 SolverMaxStep = str2double(varargin{2});
 OrderOfSolution = str2double(varargin{3});    
 elseif nargin == 21
 EONSWaveform = varargin{1};
 SolverMaxStep = str2double(varargin{2});
 OrderOfSolution = str2double(varargin{3});
 SamplingPeriod = str2double(varargin{4});  
 elseif nargin == 22
     error('Please also enter nThreshAPs');
 elseif nargin == 23
 EONSWaveform = varargin{1}; SolverMaxStep = str2double(varargin{2});
 OrderOfSolution = str2double(varargin{3}); SamplingPeriod = str2double(varargin{4});
 llActivationDetectionMode = str2double(varargin{5}); llThreshnAPs = str2double(varargin{6});
 elseif nargin == 24
 EONSWaveform = varargin{1}; SolverMaxStep = str2double(varargin{2});
 OrderOfSolution = str2double(varargin{3}); SamplingPeriod = str2double(varargin{4});
 llActivationDetectionMode = str2double(varargin{5}); llThreshnAPs = str2double(varargin{6});
 minAPdur = str2double(varargin{7});
 end
 if nargin == 25
     error('Please also include DynTol');
 end
 %nargin = 26, DynTol and utol defined in code.
 
tic;                            % Start stopwatch timer
Acc = 0;                        % Discretisation accuracy (0 = low, 1 = high)
reference = 1;                  % Voltage reference. 0: extracellular potential. 1: rest potential. 
SYMM = 0;                       % Symmetry for reps
DISPLAY = 0;
EONS = 1;                       % Coupling with EONS (evaluation of non-sinusoidal magnetic fields) solver 
% The program uses the extracellular potential reference to make all calculations, this option is
% only for plotting.

if ~exist('OrderOfSolution','var')
OrderOfSolution = 1; %1 or Reilly SENN = 2
end

% 1. Initialisation of parameters
if DISPLAY
disp('1. Initialising parameters and settings')
end
% All parameters are initialised as vectors: every vector-element 
% refers to a specific compartment (axon, dendrite, soma, ...)

% 1.1. Independent variables (These parameters can be varied in the
% SENN-model)
AHH = 3e-9;     % (m^2) Fitzhugh modification to HH-SENN (1962)
Douter = str2double(Diamet);                                 % Outer diameter (mm)
if str2double(Ltot) <= 0
Ltotal = 20.5;                                % Total (minimal) axon length (mm)
else
Ltotal = str2double(Ltot);                      % (mm)
end
SENNmodel = str2double(model); % Selection of neuronal model.1 = Hodgkin-Huxley,
%2 = FH-model, 3 = CRRSS model, 4 = SE model, 5 = SRB model 

% 1.2. Dependent or fixed variables in the SENN-model
% Convention: we start with the dendrites and end with the synapse
Nmyel = 75e7*Douter;                    % Number of myelin layers
%dinner = 0.64*Douter;                         % Inner diameter (mm)
dinner = 0.7*Douter;                         % Inner diameter (mm), Reilly SENN
%lnode = 0.0015;                              % length of a node (mm)
lnode = 0.0025;                              % length of a node (mm), Reilly SENN
linter = 100*Douter;                         % lenght of an internode (mm)
Re = 3;                                     % extracellular resistivity (Ohm*m)
if SYMM == 1
reps = 2*ceil((Ltotal-lnode)/(2*(lnode+linter))); % Number of necessary repetitions to obtain Ltotal
else
reps = ceil((Ltotal-lnode)/((lnode+linter))); % Number of necessary repetitions to obtain Ltotal
end
L = [repmat([lnode;linter],reps,1);lnode]; % Column vector with lengths (mm)
C = length(L);                                % Number of compartments
if Acc == 0
N = ones(C,1);
else
N = ceil((L./0.10));         % Number of segments vector -> We use 1 segments per 0.5 mm length of compartment
N(N>10) = 10;                   % Maximum 10 segments per compartment, otherwise calculating time is high
end
T = 1000*LLTp;                        % Simulation time (ms)

d = dinner*ones(C,1); % Diameter vector (mm)
Ra = 1.1*ones(C,1);                % Axial resistivity vector (Ohm*m), Reilly SENN
%Ra = 0.7*ones(C,1);                % Axial resistivity vector (Ohm*m) 
                                % --> Resistivity of 0.1 kOhm*cm is used
Cnode = 0;            % Capacity at ranvier nodes (F/m^2)
Cinter = 0.01/Nmyel;
%Cinter = (1.6*10^(-6))/(pi*dinner); % Capacity at internodes (F/m^2)
Cm = [repmat([Cnode;Cinter],reps,1);Cnode];  % Membrane capacitance per unit area vector (F/m^2); if 0 -> default value
                                % --> Universal membrane capacitance is used.
CFL = 1;                        % Courant number
c = 100;          % Travelspeed of electric impulses through the neuron (m/s)
%Temp = 37*ones(C,1);        % Temperature (°C)
Temp = 22*ones(C,1);        % Temperature (°C), Reilly SENN
tol = 1e-50;                    % Tolerance used in matlab fsolve (default is 1e-6)
DynTol = 1e-50;                  % Dynamic tolerance used in Matlab fsolve (if 0, use abs. tol)
if nargin == 25
tol = varargin{8};
DynTol = varargin{9};
end
% SI-units (All voltages are in mV)
L = L.*10^(-3);                 % Length (m)
T = T.*10^(-3);                 % Simulation time (s)    
d = d.*10^(-3);                 % Diameter (m)
% 1.1. Potentials (Depolarisation is positive, the potentialreference is the
% extracellular potential)
% If 0 -> default value 
Vr = zeros(C,1);   % Membrane resting potential (mV)
Vna = zeros(C,1);                       % Sodium Nernst potential (mV)
Vk = zeros(C,1);                        % Potassium Nernst potential (mV)
Vl = zeros(C,1);                        % Leakage Nernst potential (mV)

% 1.2. Conductance
% If 0 -> default value 
% In the unmyelinated axon and synapse (last 7 compartments) the channel density is doubled 
% by doubling all Hodgkin-Huxley leakage currents
Gna = zeros(C,1);                     % Sodium conductance per unit area (S/m^2)
Gk =  zeros(C,1);                     % Potassium conductance per unit area (S/m^2)
Gkslow = zeros(C,1);                  % Slow potassium conductance per unit area (S/m^2) 
Gl =  zeros(C,1);           % Leakage conductance per unit area (S/m^2)

% 1.2. Conductance
%Gminter = (1/290)/(pi*dinner); % Membrane conductance (S/m^2) of internodes
Gminter = 10/Nmyel;
Gm = [repmat([0;Gminter],reps,1);0];  % Membrane conductance per unit area (S/m^2)

% 1.3. Constants
F=96485;                         % Faraday constant [C/mol]
R=8.3144;                        % Gas constant [J/(mol.K)]

% 1.4. Time step and space step
dx = L./N;                       % Discr. step in space (m)
index = zeros(sum(N),1);            
index([1; cumsum(N(1:end-1))+1]) = 1;
index = cumsum(index); 
dxI = dx(index);
dxfI = (dxI+circshift(dxI,-1,1))./2; 
dt = min((dxfI(1:end-1)./c).*CFL);          % Discr. step in time  (s)
dt=25e-6;
if ~exist('SolverMaxStep','var')
    SolverMaxStep = 25;                 % (us) 
end
if ~exist('SamplingPeriod','var')
    SamplingPeriod = 40;
end
dt=min(1e-6*SolverMaxStep,SamplingPeriod/sinus_freq);
% -> scalar value: only one time step, because only one simulation time

% 2. Plots and processing
SweepActI = LLSweepActI;                  % if 1 the code will sweep for the minimal activation current -> furthermore other options (plots, ...) will be ignored
% 2.1. General plots
Plotf = 1;                      % if 1 the activating function f (Rattay) is calculated and plotted
VtotPlot = 1;                   % if 1 all membrane-voltage data is stored and plotted
mtotPlot = 1;                   % if 1 all m-gate data is stored and plotted
ntotPlot = 0;                   % if 1 all n-gate data is stored and plotted
htotPlot = 0;                   % if 1 all h-gate data is stored and plotted
ptotPlot = 0;                   % if 1 all p-gate data is stored and plotted
PhaseVPlot = 0;                 % if 1 the phase velocity of the activation pulse is plotted
GroupVPlot = 0;                 % if 1 the group velocity of the activation pulse is plotted
VexactPlot = 0;                 % if 1 a plot of group velocity and activation pulse path is shown
PlotVelocityOne = 0;            % Bool: if 1, only one velocity graph is made for group/phase velocity
fPhaseV = (0:5000:10000);       % frequencies for which the phase velocity is plotted (Hz)
fGroupV = (0:5000:10000);       % frequencies for which the group velocity is plotted (Hz)
% 2.2. Space plots (i.e. plots at a certain time t)
PlotSpaceOne = 1;              % Bool: if 1, only one graph per variable is produced
PlotSpaceV = 0;                % Bool: if 1, calculate space plots
PlotSpacem = 0;
PlotSpacen = 0;
PlotSpaceh = 0;
PlotSpacep = 0;
PlotVx = (1/7:1/7:1);   % Row matrix of relative times where space plots should be calculated
Plotmx = (1/7:1/7:1);   % values are chosen in ]0,1]
Plotnx = (1/7:1/7:1); 
Plothx = (1/7:1/7:1);
Plotpx = (1/7:1/7:1);

% 2.3. Time plots (i.e. plots at a certain location)
PlotTimeOne = 1;                % Bool: if 1, only one graph per variable is produced
PlotTimeV = 0;                  % Bool: if 1, calculate time plots
PlotTimem = 0; 
PlotTimen = 0; 
PlotTimeh = 0; 
PlotTimep = 0; 
PlotVt = [0.0023; 0.0272; 0.0522; 0.0771; 0.1020; 0.1270; 0.1519; 
0.1769;0.2018; 0.2268; 0.2517; 0.2766; 0.3016; 0.3265; 0.3515; 0.3764;
0.4014; 0.4263; 0.4512; 0.4762; 0.5011];   % Column matrix of relative positions where time plots should be calculated 
Plotmt = (1/7:1/7:1).';   % values are chosen in ]0,1]
Plotnt = (1/7:1/7:1).';
Plotht = (1/7:1/7:1).';
Plotpt = (1/7:1/7:1).';

% 2.4. Activating function
% The activating function can be defined as an electric field El or as an
% external potential distribution Ve. We differentiate between these two
% definitions, to optimalise the numerical solution.
ActF = 1;                       % The method to define the activating function: 0 = El, 1 = Ve
% 2.4.1. El
if ActF == 0
% 2.4.1.1. Test fields
A = 10^6;                     % Oscillation amplitude of the effective electric test field (mV/m)
f = 10000;                     %  Stimulation frequency of the testfield (Hz)
% 2.4.1.1.1. Sinusoidal test field (at 1 location along the axon)
SinusTest = 0;                % Bool: if 1 a sinusoidal test field is used
Pulse = 0;                    % Bool: if 1 (and SinusTest is 1) the test field is a sinusoidal pulse V(t) = 0 if t > T.
x0 = 0.5;                     % in [0,1[: location of the sinusoidal test field along the axon
% 2.4.1.1.2. S4L test field 
% A test field produced by S4L from a .mat file is used
S4Lfield = 1;                % Bool: if 1 the S4L testfield is used
filename = ['TotalField60Vx.mat';'TotalField60Vy.mat';'TotalField60Vz.mat']; % Total field projections along the 3 axis together with axis data
% is stored in this .mat file vector. 
% !! All filenames must have the same length -> pad blanks at the end if
% necessary
NGeom = [-0.001,0.0034;
         0.002, 0.002;
         0,0];               % Neuron geometry. NGeom = [r1,r2,r3,...,rN], with r1,r2,...,rN vectors
                             % that define N points in space to determine the neuron geometry. (The neuron is assumed to be 
                             % a piecewise broken-line. 
method = 'spline';           % The total field has to be interpolated on the right neuron positions
                             % determined by NGeom. Define the interpolation method: nearest, linear, cubic or spline.

else
% 2.4.2. Ve
% 2.4.2.1. Test potentials
% 2.4.2.1.1. Potential caused by point electrode 
x0 = LLx0s;                     % coordinates of the point electrode (m)
y0 = LLy0s;
z0 = LLz0s;
Ielec = (-1)^(LLancat-1)*LLElecInj;                          % injected current amplitude (mA)
Tp = LLTp;                        % Pulse time (s)
% 2.4.2.1.1.1. Constant current
% A constant current is injected by a point electrode
CCPE = 1;                       % Bool: if 1 a constant current is injected by a point electrode
Pulse = 1;                      % Bool: if 1 (and CCPE is 1) the injected current is a constant pulse
end
% 2.4.3. Sweep activation current values
    % Sweep ranges are defined as : [Begin point, Number of samples, end
    % point] 
x0s = [LLx0s, 1, LLx0s];    % Sweep range for x0 (m)
y0s = [LLy0s, 1, LLy0s];    % Sweep range for y0 (m)
z0s = [LLz0s,1,LLz0s];       			% Sweep range for z0 (m)
Tps = [LLTp, 1, LLTp];   % Sweep range for Tp (s) (only phase time for initial phase if biphasic) 
TpsLog = 0;                     % If 1, distribute Tps logarithmically
x0sLog = 0;                     % If 1, distribute x0s logarithmically
y0sLog = 0;                     % If 1, distribute y0s logarithmically
z0sLog = 0;                     % If 1, distribute z0s logarithmically
ElecsPrecision = 2;             % The precision (amount of significant numbers after the first digit)  that is calculated for Ielecs (minimal activation current (mA))
TimeScale = 3;              % If different from zero: dynamical timescale -> T = 2*ITps*timescale (timescale times biphasic pulseduration)
minTime = 0.001;         % Minimum simulation time (if TimeScale ~= 0)
TimeScale = 1/2; minTime = 0; % For SweepActI = 1: simulation duration equal to pulse duration
AnCat = LLancat;                  % 1= only anode, 2 = only cathode, 3 = both anode and cathode
MonBi = LLmonbi;                  % 1 = only monophasic, 2 = only biphasic, 3 = both monophasic and biphasic

% 2.5. Activated tissue
NeuronActivation = 1;           % Bool: if 1 check if neuron is activated.
ActivationTable = 1;            % Bool: if 1 show overview of activation locations and times
ActivationDetectionMode = llActivationDetectionMode; % If 1: activation detected according to conditions (1 AP), 
                                                     % If 2: activated if number APs at terminal exceeds ThreshnAPs
ThreshnAPs = llThreshnAPs;     % Number of APs to be reached at terminal, if ActivationDectectionMode = 2
TableName = 'CRRSS Cortex';    % The name of the neuron displayed in the activation table
Vtres = 80;                     % Treshold voltage (mV) for activation. % Soyka: from 50 to 80
maxc = 10000000;                     % maximum propagation speed (m/s)
Conditions = [0, 4;                  % Conditions for neuron propagation: 
              2, 2;             % -> 1 condition/row: every number refers to a number of active compartments
              4, 0;             % -> 1st column: forward propagation; 2nd column; backward propagation.
                  ];            % F.i. Conditions = [1, 5]: the activation signal is propagating, if it
                               % propagates over 1 active compartment in the forward direction, and 5 compartments in the
                                % backward direction.
Conditions = [0,2; 2,0; 1,1];  % Soyka: that should correspond to Reilly's criterium for 3 consecutive nodes
fExact = [];                   % Frequencies for which an exact estimation is made of the MSOAf and MSOAb (Hz).

% 3. Neuronal model 
Model = [repmat([SENNmodel;0],reps,1);SENNmodel];
% Selection of neuronal model. 0 = Passive TL, 1 = Hodgkin-Huxley,
%2 = FH-model, 3 = CRRSS model, 4 = SE model, 5 = SRB model -> This is also
%a vector: different neuronal models on the same neuron are possible
BoundaryConditions = 0;   % Specify the boundary conditions to finish the neuronal model
% 0 = Sealed-end boundary conditions, 1 = Voltage-clamped on rest voltage,
% 2 = voltage-clamped on specified voltage
VClamp = 0; % If BoundaryConditions = 2, specify the voltage to which the neuron is clamped

% 3.1. Neuronal model general parameters
% Here parameters for voltages, conductance and capacitance from table 1.1
% in "Functional Electrical Stimulation of the Central Nervous System (Rattay)" are used.
% All voltages are defined with respect to the extracellular potential

% 3.1.1. General 
Model0l = (Model==0);Model1l = (Model==1);Model2l = (Model==2);Model3l = (Model==3);Model4l = (Model==4);
Model5l = (Model==5);           %Logical model values (last l = logical)
Model0 = double(Model0l); Model1 = double(Model1l); Model2 = double(Model2l); Model3 = double(Model3l); 
Model4 = double(Model4l); Model5 = double(Model5l);       % Double model values
S = sum(N);            % Number of segments

Vr(Vr==0&(Model0l|Model1l|Model2l)) = -70;
Vr(Vr==0&Model3l) = -80; 
Vr(Vr==0&Model4l) = -78;
Vr(Vr==0&Model5l) = -84;

Vna(Vna==0&Model1l) = 45;
Vna(Vna==0&Model3l) = 35;

Vk(Vk==0&Model1l) = -82;
Vk(Vk==0&Model5l) = -84;

Vl(Vl==0&Model1l)= -59.4;
Vl(Vl==0&Model2l) = -69.974;
Vl(Vl==0&Model3l) = -80.01;
Vl(Vl==0&Model4l) = -78;
Vl(Vl==0&Model5l) = -84;

Gna(Gna==0&Model1l) = 1200;
Gna(Gna==0&Model3l) = 14450;

Gk(Gk==0&Model1l) = 360;
Gk(Gk==0&Model5l) = 300;

Gkslow(Gkslow==0&Model5l) = 600;

Gl(Gl==0&Model1l) = 3;
Gl(Gl==0&Model2l) = 303;
Gl(Gl==0&Model3l) = 1280;
Gl(Gl==0&Model4l) = 860;
Gl(Gl==0&Model5l) = 600;

Cm(Cm==0&Model0l) = 0.01;          % This value is not based on Rattay.
Cm(Cm==0&Model1l) = 0.01;
Cm(Cm==0&Model2l) = 0.02;
%Cm(Cm==0&Model3l) = 0.025;
Cm(Cm==0&Model3l) = 0.01;
Cm(Cm==0&Model4l) = 0.028;
Cm(Cm==0&Model5l) = 0.028;
            
%Additional parameters in FH-model, SE-model and SRB-model:
Pna = zeros(C,1); Pk=zeros(C,1); Pp = zeros(C,1);
cNao = zeros(C,1); cNai = zeros(C,1); cKo = zeros(C,1);
cKi = zeros(C,1);

Pna(Model2l) = 0.00008;        % Pna (m/s)
Pk(Model2l) = 0.000012;        % Pk (m/s)
Pp(Model2l) = 0.0000054;       % Pp (m/s)
cNao(Model2l) = 114.5;        % [Na]o [mol/m^3]
cNai(Model2l) = 13.7;        % [Na]i [mol/m^3]
cKo(Model2l) = 2.5;         % [K]o [mol/m^3]
cKi(Model2l) = 120;          % [K]i [mol/m^3]

Pna(Model4l) = 0.0000328;
Pk(Model4l) = 0.00000134;
cNao(Model4l) = 154;
cNai(Model4l) = 8.71;
cKo(Model4l) = 5.9;
cKi(Model4l) = 155;

Pna(Model5l) = 0.0000704;
cNao(Model5l) = 154;
cNai(Model5l) = 30;
 
% Temperature coefficients
Qam = 3.*ones(C,1); Qbm = 3.*ones(C,1); Qan = 3.*ones(C,1);  
Qbn = 3.*ones(C,1); Qah = 3.*ones(C,1); Qbh = 3.*ones(C,1);
Qap = 3.*ones(C,1); Qbp = 3.*ones(C,1);

Qam(Model2l) = 1.8;
Qbm(Model2l) = 1.7;
Qan(Model2l) = 3.2;
Qbn(Model2l) = 2.8;
Qah(Model2l) = 2.8;
Qbh(Model2l) = 2.9;

Qam(Model4l) = 2.2;
Qbm(Model4l) = 2.2;
Qah(Model4l) = 2.9;
Qbh(Model4l) = 2.9;

Qam(Model5l) = 2.2;
Qbm(Model5l) = 2.2;
Qah(Model5l) = 2.9;
Qbh(Model5l) = 2.9;
 
% 4. Associated parameters
% 4.1. Compartmental parameters (c = compartmental)
Cmc = (pi.*d.*dx.*Cm).*double(Model~=1)+((AHH/20).*Cm).*double(Model==1);              % Compartmental capacitance vector (F)
Gac = (pi.*d.^2)./(4.*dx.*Ra);    % Compartmental conductance vector (S)

NeuronActivated = 0; TerminalActivated = 0;
NumberSimulations = 0; teller = 0; reverseStr = '';

if x0s(1) == x0s(3)
x0sValues = (x0s(1));  
else
if x0sLog 
x0sValues =(x0s(3)/x0s(1)).^((0:1:(x0s(2)-1))/(x0s(2)-1))*x0s(1);
else
xStep = (x0s(3)-x0s(1))/(x0s(2)-1);
x0sValues = (x0s(1):xStep:x0s(3));
end
end
NxStep = length(x0sValues);

if y0s(1) == y0s(3)
y0sValues = (y0s(1));
else
if y0sLog
y0sValues = (y0s(3)/y0s(1)).^((0:1:(y0s(2)-1))/(y0s(2)-1))*y0s(1);
else
yStep = (y0s(3)-y0s(1))/(y0s(2)-1); 
y0sValues = (y0s(1):yStep:y0s(3));
end
end
NyStep = length(y0sValues);

if z0s(1) == z0s(3)
z0sValues = (z0s(1));
else
if z0sLog
z0sValues = (z0s(3)/z0s(1)).^((0:1:(z0s(2)-1))/(z0s(2)-1))*z0s(1);
else
zStep = (z0s(3)-z0s(1))/(z0s(2)-1); 
z0sValues = (z0s(1):zStep:z0s(3));
end
end
NzStep = length(z0sValues);

if Tps(1) == Tps(3)
TpsValues = (Tps(1));
else
if TpsLog
TpsValues = (Tps(3)/Tps(1)).^((0:1:(Tps(2)-1))/(Tps(2)-1))*Tps(1);
else
TStep = (Tps(3)-Tps(1))/(Tps(2)-1);
TpsValues = (Tps(1):TStep:Tps(3));
end
end
NTStep = length(TpsValues);

if SweepActI
dts = min(dt,TpsValues./75);        % dt must be small enough w.r.t. ITp
else if ActF == 1 && Pulse == 1 && CCPE == 1
    dts = min(dt,Tp/75);
    else
    dts = dt;
    end
end
if SweepActI && TimeScale ~= 0
Nts = max(floor(2*TimeScale.*TpsValues./dts),floor(minTime./dts));
else
Nts = floor(T./dts);
end
T = Nts.*dts;

SweepActImat = zeros(NxStep,NyStep,NzStep,NTStep,2,2);    
if SweepActI
    if TimeScale ~= 0
TotalIt = 4*NxStep*NyStep*NzStep*sum(Nts);
    else
TotalIt = 4*NxStep*NyStep*NzStep*NTStep*sum(Nts);
    end
else
TotalIt = sum(Nts);
end

if AnCat == 1
AnCatLims = [1, 1];
else if AnCat == 2
     AnCatLims = [2, 2];
    else if AnCat == 3
            AnCatLims = [1, 2];
        end
    end
end

if MonBi == 1
MonBiLims = [1, 1];
else if MonBi == 2
        MonBiLims = [2, 2];
    else if MonBi == 3
            MonBiLims = [1, 2];
        end
    end
end

for Inx0s = 1:NxStep
    Ix0s = x0sValues(Inx0s);    
for Iny0s = 1:NyStep
    Iy0s = y0sValues(Iny0s);  
for Inz0s = 1:NzStep
    Iz0s = z0sValues(Inz0s);
for InTps = 1:NTStep
    ITps = TpsValues(InTps);
    Nt = Nts(InTps);
    dt = dts(InTps);
for Iancat = AnCatLims(1):AnCatLims(2)
for Imb = MonBiLims(1):MonBiLims(2)  
flg2 = 0;               % set flag 2
Ielecs = [LLIbegin LLIend]; PrecisionCheck = 0; SearchMode = LLSearchMode; % Ielecs in mA
if SearchMode % SearchMode = 1
IIelecs = (Ielecs(1)+Ielecs(2))/2;
else % SearchMode = 0
if ~any(Ielecs)
IIelecs = 1;
else
IIelecs = Ielecs(~~Ielecs)*10^(3-2*find(Ielecs));
end
end
while ~SearchMode || ~PrecisionCheck 
NumberSimulations = NumberSimulations+1;

% 4.2. Recorders
% For vectors and matrices the convention is used that rows refer to space
% and columns to time

%4.2.1. Time recorders
timeflow = zeros(1,Nt+1);
if ~SweepActI
if PlotTimeV
recorderV = zeros(length(PlotVt),Nt+1);      % Recorders in time for V,h,p,n,m
end
if PlotTimen
recordern = zeros(length(Plotnt),Nt+1);
end
if PlotTimeh
recorderh = zeros(length(Plotht),Nt+1);
end
if PlotTimem
recorderm = zeros(length(Plotmt),Nt+1); %#ok<*UNRCH>
end
if PlotTimep
recorderp = zeros(length(Plotpt),Nt+1);
end
end
% 4.2.2. Space recorders
if ~SweepActI
if PlotSpaceV
recorderVx = zeros(S,length(unique(ceil((Nt+1).*PlotVx))));      % Recorders in space for V,h,p,n,m
end                                                              % Here it is necessary that we don't use 2 times the same time (-> unique)
if PlotSpacen
recordernx = zeros(S,length(unique(ceil((Nt+1).*Plotnx))));
end
if PlotSpaceh
recorderhx = zeros(S,length(unique(ceil((Nt+1).*Plothx))));
end
if PlotSpacem
recordermx = zeros(S,length(unique(ceil((Nt+1).*Plotmx))));
end
if PlotSpacep
recorderpx = zeros(S,length(unique(ceil((Nt+1).*Plotpx))));
end
end

% 4.2.3. Total field recorders
if VtotPlot || NeuronActivation || SweepActI
recorderVtot = zeros(S,Nt+1);     % Stores all data of the membrane-voltage
end
if ~SweepActI
if mtotPlot 
    recordermtot = zeros(S,Nt+1);
end
if htotPlot 
    recorderhtot = zeros(S,Nt+1);
end
if ntotPlot
    recorderntot = zeros(S,Nt+1);
end
if ptotPlot 
    recorderptot = zeros(S,Nt+1);
end
end
% 5. Initialisation of variables
% To initialise a variable as V, we need a Sx1 vector, with Vr(1) at
% the first N(1) entries, Vr(2) at the following N(2) entries of the
% vector, ...
index = zeros(S,1);            % 1) An index vector with the same dimension as V is initialised
index([1; cumsum(N(1:end-1))+1]) = 1; % 2) 1 is declared at each index where a new sequence of the same V value will start
index = cumsum(index);
V = Vr(index);              % 3) Membrane potential is defined (mV)             

% 5.1. Definition 'indexed' variables and parameters: 
Model0I = Model0(index); Model1I = Model1(index); Model2I = Model2(index); Model3I = Model3(index); 
Model4I = Model4(index); Model5I = Model5(index);

GacI = Gac(index); CmcI = Cmc(index); dxI = dx(index); GmI = Gm(index); dI = d(index);

VrI = Vr(index); VnaI =  Vna(index); VkI = Vk(index); VlI = Vl(index);                                 
GnaI = Gna(index); GkI = Gk(index); GkslowI = Gkslow(index); GlI = Gl(index);  
PnaI = Pna(index); PkI=Pk(index); PpI = Pp(index);
cNaoI = cNao(index); cNaiI = cNai(index); cKoI = cKo(index);
cKiI = cKi(index);

QamI = Qam(index); QbmI = Qbm(index); QanI = Qan(index); QbnI = Qbn(index);QahI = Qah(index); 
QbhI = Qbh(index); QapI = Qap(index); QbpI = Qbp(index);
TempI = Temp(index);

% Theoretical phase and group velocities
Vphase = zeros(S,Nt+1,length(fPhaseV)); % initialising Vphase (m/s) with dimensions space, time and frequency
Vgroup = zeros(S,Nt+1,length(fGroupV)); % initialising Vgroup (m/s)
Vexact = zeros(S,Nt+1,length(fExact)); % Initialising Vexact (m/s) (= Vgroup at the exact frequencies of activation table)

% The compartmental conductance vector has to be calculated at the
% staggered grid locations, we define:
dxfI = (dxI+circshift(dxI,-1,1))./2;               % Forward spatial discr vector (m)
dxbI = (dxI+circshift(dxI,1,1))./2;                % Backward spatial discr vector (m)

% Staggered grid locations:
GacfI = 2./(1./GacI+1./circshift(GacI,-1,1));  % Forward compartmental conductance vector (S)
GacbI = 2./(1./GacI+1./circshift(GacI,1,1));   % Backward compartmental conductance vector (S)
GacmI = circshift(GacfI,1,1)+circshift(GacbI,-1,1); % Middle compartmental conductance vector (S)

% x-axis
x = dxI./2+[0; cumsum(dxI(1:end-1))];

% 5.2 Neuron model parameters
% 5.2.1 Rate parameters
% Here the rate parameters of the gate-parameters m,n,h,p are defined, for each neuron model: am,
% bm, an, bn, ah, bh, ap, bp.
% We omit a factor 10^3 at each rate constant, which we add again in the end. (this is the conversion factor ms -> s)            
am1I = @(V) QamI.^(0.1.*TempI-0.63).*(2.5-0.1.*(V-VrI))./(exp(2.5-0.1.*(V-VrI))-1);
bm1I = @(V) QbmI.^(0.1.*TempI-0.63).*4.*exp(-(V-VrI)./18);
an1I = @(V) QanI.^(0.1.*TempI-0.63).*(0.1-0.01.*(V-VrI))./(exp(1-0.1.*(V-VrI))-1);
bn1I = @(V) QbnI.^(0.1.*TempI-0.63).*0.125.*exp(-(V-VrI)./80);
ah1I = @(V) QahI.^(0.1.*TempI-0.63).*0.07.*exp(-(V-VrI)./20);
bh1I = @(V) QbhI.^(0.1.*TempI-0.63).*(1./(exp(3-0.1.*(V-VrI))+1));

am2I = @(V) QamI.^(0.1.*TempI-2).*(0.36.*((V-VrI)-22))./(1-exp((22-(V-VrI))./3));
bm2I = @(V) QbmI.^(0.1.*TempI-2).*(0.4.*(13-(V-VrI)))./(1-exp(((V-VrI)-13)./20));
an2I = @(V) QanI.^(0.1.*TempI-2).*(0.02.*((V-VrI)-35))./(1-exp((35-(V-VrI))./10));
bn2I = @(V) QbnI.^(0.1.*TempI-2).*(0.05.*(10-(V-VrI)))./(1-exp(((V-VrI)-10)./10));
ah2I = @(V) -QahI.^(0.1.*TempI-2).*(0.1.*((V-VrI)+10))./(1-exp(((V-VrI)+10)./6));
bh2I = @(V) QbhI.^(0.1.*TempI-2).*(4.5)./(1+exp((45-(V-VrI))./10));
ap2I = @(V) QapI.^(0.1.*TempI-2).*(0.006.*((V-VrI)-40))./(1-exp((40-(V-VrI))./10));
bp2I = @(V) -QbpI.^(0.1.*TempI-2).*0.09.*((V-VrI)+25)./(1-exp(((V-VrI)+25)./20));
    
am3I = @(V) QamI.^(0.1.*TempI-3.7).*(97+0.363.*(V-VrI))./(1+exp((31-(V-VrI))./5.3));
bm3I = @(V) QbmI.^(0.1.*TempI-3.7).*(97+0.363.*(V-VrI))./((1+exp((31-(V-VrI))./5.3)).*(exp(((V-VrI)-23.8)./4.17))); 
bh3I = @(V) QbhI.^(0.1.*TempI-3.7).*15.6./(1+exp((24-(V-VrI))./10));  
ah3I = @(V) QahI.^(0.1.*TempI-3.7).*15.6./((1+exp((24-(V-VrI))./10)).*(exp(((V-VrI)-5.5)./5)));

am4I = @(V) QamI.^(0.1.*TempI-3.7).*1.87.*((V-VrI)-25.41)./(1-exp((25.41-(V-VrI))./6.06));
bm4I = @(V) QbmI.^(0.1.*TempI-3.7).*3.97.*(21-(V-VrI))./(1-exp(((V-VrI)-21)./9.41));
an4I = @(V) QanI.^(0.1.*TempI-3.7).*0.13.*((V-VrI)-35)./(1-exp((35-(V-VrI))./10));
bn4I = @(V) QbnI.^(0.1.*TempI-3.7).*0.32.*(10-(V-VrI))./(1-exp(((V-VrI)-10)./10));
ah4I = @(V) -QahI.^(0.1.*TempI-3.7).*0.55.*((V-VrI)+27.74)./(1-exp(((V-VrI)+27.74)./9.06));
bh4I = @(V) QbhI.^(0.1.*TempI-3.7).*22.6./(1+exp((56-(V-VrI))./12.5));
    
am5I = @(V) QamI.^(0.1.*TempI-3.7).*4.6.*((V-VrI)-65.6)./(1-exp((-(V-VrI)+65.6)./10.3));
bm5I = @(V) QbmI.^(0.1.*TempI-3.7).*0.33.*(61.3-(V-VrI))./(1-exp(((V-VrI)-61.3)./9.16));
an5I = @(V) QanI.^(0.1.*TempI-3.7).*0.0517.*((V-VrI)+9.2)./(1-exp((-(V-VrI)-9.2)./1.1));
bn5I = @(V) QbnI.^(0.1.*TempI-3.7).*0.092.*(8-(V-VrI))./(1-exp(((V-VrI)-8)./10.5));
ah5I = @(V) -QahI.^(0.1.*TempI-3.7).*0.21.*((V-VrI)+27)./(1-exp(((V-VrI)+27)./11));
bh5I = @(V) QbhI.^(0.1.*TempI-3.7).*14.1./(1+exp((55.2-(V-VrI))./13.4));
ap5I = @(V) QapI.^(0.1.*TempI-3.7).*0.0079.*((V-VrI)-71.5)./(1-exp((71.5-(V-VrI))./23.6));
bp5I = @(V) -QbpI.^(0.1.*TempI-3.7).*0.00478.*((V-VrI)-3.9)./(1-exp(((V-VrI)-3.9)./21.8));

amI = @(V) 10^3.*(am1I(V).*Model1I+am2I(V).*Model2I+am3I(V).*Model3I+am4I(V).*Model4I+am5I(V).*Model5I);
bmI = @(V) 10^3.*(bm1I(V).*Model1I+bm2I(V).*Model2I+bm3I(V).*Model3I+bm4I(V).*Model4I+bm5I(V).*Model5I);
anI = @(V) 10^3.*(an1I(V).*Model1I+an2I(V).*Model2I+an4I(V).*Model4I+an5I(V).*Model5I);
bnI = @(V) 10^3.*(bn1I(V).*Model1I+bn2I(V).*Model2I+bn4I(V).*Model4I+bn5I(V).*Model5I);
ahI = @(V) 10^3.*(ah1I(V).*Model1I+ah2I(V).*Model2I+ah3I(V).*Model3I+ah4I(V).*Model4I+ah5I(V).*Model5I);
bhI = @(V) 10^3.*(bh1I(V).*Model1I+bh2I(V).*Model2I+bh3I(V).*Model3I+bh4I(V).*Model4I+bh5I(V).*Model5I);
apI = @(V) 10^3.*(ap2I(V).*Model2I+ap5I(V).*Model5I);
bpI = @(V) 10^3.*(bp2I(V).*Model2I+bp5I(V).*Model5I);

% 5.2.2. Initial conditions
m = amI(VrI)./(amI(VrI)+bmI(VrI));       % Remark: everywhere where m0 isn't defined, the division will result in NaN
h = ahI(VrI)./(ahI(VrI)+bhI(VrI));
n = anI(VrI)./(anI(VrI)+bnI(VrI));       % n0 isn't defined on model 3 -> on the compartments with model 3 the division will result in NaN
p = apI(VrI)./(apI(VrI)+bpI(VrI));

m(isnan(m)|isinf(m))=0;           % Replace NaN -> 0
h(isnan(h)|isinf(h))=0;
n(isnan(n)|isinf(n))=0;
p(isnan(p)|isinf(p))=0;
m0 = m; p0 = p; n0 = n; h0 = h;         % The variables m0, p0, n0,h0 are defined for debugging only 
% 6.1. The first values of V,m,n,h,p are recorded (-> the plots will show V,m,h,n,p
% starting from their initial values)

timeflow(1,1) = 0;
if ~SweepActI
if PlotTimeV 
recorderV(:,1) = V(ceil(S.*PlotVt(:,1)));
end

if PlotSpaceV 
    sortPlotVx = sort(PlotVx);
     if ceil(sortPlotVx(1)*(Nt+1))==1
recorderVx(:,1) = V;         
    end
end
end
if VtotPlot || NeuronActivation || SweepActI
    recorderVtot(:,1) = V;
end

if ~SweepActI
if sum(Model) ~= 0              % If Model is always 0, plotting gate-parameters is meaningless
if PlotTimem 
recorderm(:,1) = m(ceil(S.*Plotmt(:,1)));
end

if PlotSpacem 
    sortPlotmx = sort(Plotmx);
    if ceil(sortPlotmx(1)*(Nt+1)) == 1
        recordermx(:,1) = m;
    end
end

if mtotPlot 
    recordermtot(:,1) = m;
end

if sum(~(Model == 3)) ~= 0          % Only plot n, if there is at least one different from Model 3
    if PlotTimen 
recordern(:,1) = n(ceil(S.*Plotnt(:,1)));
    end

if PlotSpacen
    sortPlotnx = sort(Plotnx);
    if ceil(sortPlotnx(1)*(Nt+1))==1
        recordernx(:,1) = n;
    end
end

    
if ntotPlot
    recorderntot(:,1) = n;
end
end

if PlotTimeh 
recorderh(:,1) = h(ceil(S.*Plotht(:,1)));
end

if PlotSpaceh 
    sortPlothx = sort(Plothx);
    if ceil(sortPlothx(1)*(Nt+1))==1
        recorderhx(:,1) = h;
    end
end

if htotPlot
    recorderhtot(:,1) = h;
end

if sum((Model == 2)|(Model == 5)) ~= 0                % Only plot, if there is at least one model 2 or 5
    if PlotTimep 
recorderp(:,1) = p(ceil(S.*Plotpt(:,1)));
    end
    
if PlotSpacep 
    sortPlotpx = sort(Plotpx);
    if ceil(sortPlotpx(1)*(Nt+1))==1
        recorderpx(:,1) = p;
    end
end    
    
if ptotPlot 
    recorderptot(:,1) = p;
end
end
end
end

if SweepActI 
Ve = zeros(S,Nt);
else
if ActF == 0
E = zeros(S-1,Nt);                 % Effective electric field (mV/m) (from t = dt/2 -> t = Nt*dt-dt/2 = T-dt/2)
% E,m,n,p,h ->  shifted time steps  
% V -> whole time steps
else
Ve = zeros(S,Nt);               % External potential (mV) (from t=dt/2 -> t=Nt*dt-dt/2 = T-dt/2 -> Ve and V are not allocated in time!)
end
end

% 6.2.Test fields
% 6.2.1. El
if ~SweepActI
if ActF == 0
% Here test fields E are defined to test the code. 
% 6.2.1.1. Sinusoidal test field
if SinusTest 
% A sinusoidal activation function is used
% while El = 0 at other locations
E(floor(x0*S),:) = A.*sin(2.*pi.*f.*dt.*((0:1:Nt-1)+1/2));    % Sinusoidal excitation 
if Pulse
    Tp = 1/f;
    Np = ceil(Tp/dt-1/2);       % Np*dt+dt/2 = Tp
    E(floor(x0*S),Np:end) = 0;
end
end 

TIME1 = toc;
% 6.2.1.2. S4L field
if S4Lfield 
if DISPLAY    
disp('1.1. Loading S4L electric field...');
end
filenameCell = cellstr(filename);      % Convert filename to cell of strings (removes blanks at the end)
dots = size(NGeom,2); edges = dots-1;  % number of edges, number of endpoint dots

for i = 1:3
TotalFieldstruct = load(filenameCell{i});         % Load the total field from the S4L .mat file. 
if i == 1  
    Snapshot=zeros(length(real(TotalFieldstruct.Snapshot0)),3); % preallocate in first iteration step 
end
Snapshot(:,i) = real(TotalFieldstruct.Snapshot0);  % Snapshot holds the x-, y- and z-projection of the S4L field
end
Axis0 = TotalFieldstruct.Axis0; Axis1 = TotalFieldstruct.Axis1; Axis2 = TotalFieldstruct.Axis2;
% Axisi (i=0:2) determines the x,y and z-axis. Snapshot is the 
% total field discretised according to Axisi. The same S4L
% discretisation for all field projections is assumed.
% Snapshot is expressed in  V/m instead of mV/m. This is corrected at the end, 
% instead of right away, to increase the stability of the interp3 matlab function.
if DISPLAY
disp('-> S4L file has loaded.'); 
end
S4LTIME = toc-TIME1;
tic;

% 6.2.1.2.1. Determination of the query points Xq, Yq, Zq to which the effective electric S4L field has to be interpolated                            
Dir = circshift(NGeom,-1,2)-NGeom;  % Definition of the matrix with direction vectors along the neuron
Dir = Dir(:,1:end-1);
Neuronl = sqrt(sum(Dir.^2,1));       % lengths between the endpoints of a segment (l(i) = d(x(i),x(i+1)))
Dir = Dir./repmat(Neuronl,3,1);                  % Normalisation of the directionvector.
CNeuronl = cumsum(Neuronl); CNeuronl = [0, CNeuronl]; %#ok<AGROW> % total length to begin point
Edge = (repmat(CNeuronl(1:end-1),S,1)<repmat(x,1,edges))&(repmat(x,1,edges)<=repmat(CNeuronl(2:end),S,1));             
% Edge is a S x edges logical matrix: Edge(i,j) = 1 iff segment i is on
% edge j

NEURONcoord = zeros(3,S);
for i=1:dots-1
NEURONcoord = NEURONcoord+double(repmat(Edge(:,i)',3,1)).*...
    (repmat(NGeom(:,i),1,S)+Dir(:,i)*(x'-CNeuronl(i)));
end
Xq = NEURONcoord(1,1:end-1)'; Yq = NEURONcoord(2,end-1)'; Zq = NEURONcoord(3,end-1)';

% 6.2.1.2.2. interpolate Axis0, Axis1, Axis2 to gridpositions of Snapshot0 -> AxisX,
% AxisY, AxisZ
AxisX = (Axis0+circshift(Axis0,-1,2))./2; AxisX = AxisX(1:length(Axis0)-1);
AxisY = (Axis1+circshift(Axis1,-1,2))./2; AxisY = AxisY(1:length(Axis1)-1);
AxisZ = (Axis2+circshift(Axis2,-1,2))./2; AxisZ = AxisZ(1:length(Axis2)-1);

[AXISX,AXISY,AXISZ] = meshgrid(AxisX,AxisY,AxisZ);          % Define gridmatrices (meshgrid is used (<-> ndgrid) because interp3 is used later on)
nx = length(AxisX); ny = length(AxisY); nz = length(AxisZ);
% 6.2.1.2.3. Interpolate the effective electric field from the 3D S4L Snapshot field
% -> 2D neuron gridlocations determined by Xq, Yq, Zq.

% The S4L electric field stored in Snapshot is a linear indexed column
% matrix of the 3D field -> transform to matrix.
SnapshotM = zeros(nx,ny,nz,3);
EP = zeros(S-1,3); corr = zeros(3,1);
if DISPLAY
disp('1.2. Interpolating S4L electric field');
end
for i = 1:3
SnapshotM(:,:,:,i) = permute(reshape(Snapshot(:,i),nx,ny,nz),[2 1 3]);    % We need to permute because meshgrid defines x-axis along increasing columns      
if i == 1 
proj = 'x';
else if i == 2
        proj = 'y';
    else proj = 'z';
    end
end
if DISPLAY
disp(['1.2.' num2str(i) ' ' proj '-projection']);
end
Log = 1; k = 0; MaxIt = 10; 
while Log && k<MaxIt  % Keep iterating as interpolated E-field is singular
SnapshotMt=SnapshotM(:,:,:,i)./(10^k);
EP(:,i) = interp3(AXISX,AXISY,AXISZ,SnapshotMt,Xq,Yq,Zq,method);
Log = logical(sum(isnan(EP(:,i))|isinf(EP(:,i))));
k=k+1;
end
corr(i) = k;
if Log                                % Stop iteration if a singularity occurs (to simplify debugging)
error(['Interpolated S4L field is singular after ' num2str(k) ' iteration(s)']);
else
if DISPLAY
disp(['-> Interpolation succesful after ' num2str(k) ' iteration(s)']);
end
end
end
% 1) EPc = gmultiply(1000.*10.^(corr-1)',EP); % Correct the electric field (EPc (mV/m))
% 2) Projection = (EPc*Dir); % S-1x edges matrix : Projection_ij -> segment i projected according to edge j
% 3) filter = double(Edge).*(Projection); % filter the right edge (Sxedges matrix)
% 4) E1time = sum(filter,2); % sum over columns to create the right projected electric field for 1 time
% 5) E = repmat(E1time,1,Nt); % calculate E-field. 
% % We do previous steps in one line, to reduce the amount of stored
% % variables:
E = repmat(sum(double(Edge(1:end-1,:)).*((gmultiply(1000.*10.^(corr-1)',EP))*Dir),2),1,Nt);
end
else
% 6.2.2. Ve
if CCPE
    if exist(fileName) %#ok<EXIST>
fnameStruct = load(fileName);
Waveform = fnameStruct.s.Waveform; 
Waveform = interp1(Waveform(:,1),Waveform(:,2),(0:dt:Nt*dt)','linear','extrap');
Ve = Ielec.*gmultiply(fnameStruct.s.V{LLconf},Waveform');
if Pulse
Np = ceil(Tp/dt-1/2);              % Np*dt+dt/2 = Tp
Ve(:,Np:end) = 0;
end
    else
% A constant current is applied by a point electrode at (x0,y0,z0)
% Ve = Re*Ielec/(4*pi*r), we assume the neuron to be aligned with the
% x-axis
r = @(x) sqrt((x-x0).^2+y0^2+z0^2);
VeP = Re.*Ielec./(4.*pi.*r(dxI./2+[0; cumsum(dxI(1:end-1))]));
Ve = repmat(VeP,1,Nt);
if Pulse
Np = ceil(Tp/dt-1/2);              % Np*dt+dt/2 = Tp
Ve(:,Np:end) = 0;
end
    end
end
if EONS
Waveform = interp1(EONSWaveform(:,1),EONSWaveform(:,2)./max(abs(EONSWaveform(:,2))),(0:dt:Nt*dt),'linear',0);
Ve = -(-1)^(Iancat-1).*Ielec.*gmultiply((dxI./2+[0; cumsum(dxI(1:end-1))]),Waveform);
end
end
else 			% SweepActI = 1
    if ~EONS  % no EONS
    if exist(fileName) %#ok<EXIST>    
fnameStruct = load(fileName);
Waveform = fnameStruct.s.Waveform; 
Waveform = interp1(Waveform(:,1),Waveform(:,2),(0:dt:Nt*dt)','linear','extrap');
Ve = (-1)^(Iancat-1).*IIelecs.*gmultiply(fnameStruct.s.V{LLconf},Waveform');
if Pulse
Nps = ceil(ITps/dt-1/2);              % Np*dt+dt/2 = Tp
if Nps > Nt
    error('Insufficient simulation time!');
end
if (Imb-1)==0
Ve(:,Nps:end) = 0; 
else
Ve(:,Nps:2*Nps-2) = -Ve(:,Nps:2*Nps-2);          % Biphasic pulse: initial phasetime matches with monophasic stimulation -> double total stimulationtime
Ve(:,2*Nps-1:end) = 0;
end
end
else
r = @(x) sqrt((x-Ix0s).^2+Iy0s^2+Iz0s^2);
VeP = (-1)^(Iancat-1).*Re.*IIelecs./(4.*pi.*r(dxI./2+[0; cumsum(dxI(1:end-1))]));
Nps = ceil(ITps/dt-1/2);
Ve = repmat(VeP,1,Nt);
if Nps > Nt
    error('Insufficient simulation time!');
end
if (Imb-1)==0
Ve(:,Nps:end) = 0; 
else
Ve(:,Nps:2*Nps-2) = -Ve(:,Nps:2*Nps-2);          % Biphasic pulse: initial phasetime matches with monophasic stimulation -> double total stimulationtime
Ve(:,2*Nps-1:end) = 0;
end
    end
    else % EONS 
Waveform = interp1(EONSWaveform(:,1),EONSWaveform(:,2)./max(abs(EONSWaveform(:,2))),(0:dt:Nt*dt),'linear',0);
Ve = -(-1)^(Iancat-1).*IIelecs.*gmultiply((dxI./2+[0; cumsum(dxI(1:end-1))]),Waveform);
    end
end

% 6.3. Crank-Nicholson matrices
DiagA0 = [-GacfI GacmI -GacbI];   % Interpolated axial conductance values are used
A0 = spdiags(DiagA0,[-1 0 1],S,S); 

if BoundaryConditions == 0
% Adding Neumann boundary conditions:
% dV/dx(0) = dV/dx(L) = 0 -> V(-1) = V(1) and V(L+1)=V(L-1) 

%Reilly SENN: remove factor 2: factor 2 corresponds with halving of the capacitance if the sealed end is at the center of the Ranviernode
% begin point:
%A0(1,1) = 2*GacfI(1); %#ok<*SPRIX>
%A0(1,2) = -2*GacbI(2);
A0(1,1) = GacfI(1); %#ok<*SPRIX>
A0(1,2) = -1*GacbI(2);

% end point:
%A0(S,S-1) = -2*GacfI(S-1);
%A0(S,S) = 2*GacbI(S);
A0(S,S-1) = -1*GacfI(S-1);
A0(S,S) = GacbI(S);

else
if BoundaryConditions == 1 || BoundaryConditions == 2
% Adding Dirichlet boundary conditions:
% begin point:
A0(1,1) = 0;
A0(1,2) = 0;
% end point:
A0(S,S-1) = 0;
A0(S,S) = 0;    
end
end

B = diag(2*CmcI./dt); 
if BoundaryConditions == 1 || BoundaryConditions == 2
    B(1,1) = 1; B(S,S) = 1;
end
A1 = A0+B;

if ~SweepActI
if ActF == 0                % Activating function = El
diagC = [GacfI.*dxfI -GacfI.*dxfI zeros(S,1)];
C = spdiags(diagC,[-1 0 1],S,S-1);
if BoundaryConditions == 0
% Adding Neumann boundary conditions
C(1,1) = -2*GacfI(1)*dxfI(1);                       
C(S,S-1) = 2*GacfI(S-1)*dxfI(S-1);   
else if BoundaryConditions == 1 || BoundaryConditions == 2
C(1,1) = 0;
C(S,S-1) = 0;
    end
end
else                        % Activating function = Ve
C = -A0;
end
else
C = -A0;
end

% 6.4. Calculate activating function f (Rattay)
if ~SweepActI
if Plotf&&ActF
Ve1 = Ve(:,1);          % External potential at t = dt/2 (the 1 subscript refers to the time, begin/end refer to space)
f = ((-GacfI-GacbI).*Ve1+GacbI.*circshift(Ve1,1,1)+GacfI.*circshift(Ve1,-1,1))./CmcI;  %f (mV/s)
% Begin and endpoint
f(1) = 2*GacfI(1)*(Ve1(2)-Ve1(1))/CmcI(1);  %f(1) (mV/s)
f(S) = 2*GacbI(S)*(Ve1(S-1)-Ve1(S))/CmcI(S);  %f(S) (mV/s)
end
if Plotf&&~ActF
E1 = E(:,1);
f = [(-GacfI(1:end-1).*dxfI(1:end-1).*E1+GacbI(1:end-1).*dxbI(1:end-1).*circshift(E1,1,1))./CmcI(1:end-1);0];  % f (mV/s)
% Begin and endpoint
f(1) = -GacfI(1)*dxfI(1)*E1(1)/CmcI(1);
f(S) = GacfI(S-1)*dxfI(S-1)*E1(S-1)/CmcI(S);
f = (10^(-3)).*f; % f (mV/ms) or f (V/s)
end
end 

% 7. Calculate potential and gate parameters
if ~SweepActI
if DISPLAY
disp('2. Calculating membrane potential and gate parameters');
end
end
% A) V is recorded from t = 0*dt -> Nt*dt = T, i.e. Nt+1 time steps, with V(0)
% = Vr.
% B) n,m,h,p are recorded from t = dt/2 -> T+dt/2 (Nt+1 time steps), 
% with n(dt/2) = n0 (the initial value), m(dt/2) = m0,...
% C) El is defined on Nt time steps: from t = dt/2 -> T-dt/2
% D) All velocities (Vphase,Vgroup,Vexact) are calculated on the same times
% as the gate parameters.

SingularityPass = 0;        % If signals become singular, this Pass will become 1, to skip processing code
for it = 1:Nt
    if ~flg2 
   teller = teller+1;
    end
   if DISPLAY
   progress = (teller/TotalIt)*100;      % Shows progress (%)
   msg = sprintf('Progress: %3.1f', progress); 
   fprintf([reverseStr, msg]);
   reverseStr = repmat(sprintf('\b'), 1, length(msg)); 
   end
   t = dt*it;
timeflow(1,it+1) = t;             % timeflow follows the timeflow of V, the timeflow of n,m,... is shifted with dt/2
if ~SweepActI
if ActF == 0
AF = E(:,it);                   % General activating function (AF) 
else
AF = Ve(:,it);
end
else
AF = Ve(:,it);
end

% 7.1. Update membrane voltage
% Membrane voltage is updated in 2 steps
% 7.1.1. Update voltage t -> t+dt/2
% 7.1.1.1. Compartment membrane-conductance and compartment voltage matrix
                               
Gmc0 = pi.*dI.*dxI.*GmI;             % Passive compartment membrane conductance 
Vc0 = pi.*dI.*dxI.*GmI.*VrI;       % Passive compartment voltage matrix
Gmc1 = AHH.*(GnaI.*m.^3.*h+GkI.*n.^4+GlI);  %HH compartment membrane conductance
Vc1 = AHH.*((GnaI.*VnaI).*m.^3.*h+(GkI.*VkI).*n.^4+(GlI.*VlI));  % HH compartment voltage matrix
Gmc3 = pi.*dI.*dxI.*(GnaI.*m.^2.*h+GlI);  %CRRSS compartment membrane conductance
Vc3 = pi.*dI.*dxI.*((GnaI.*VnaI).*m.^2.*h+(GlI.*VlI));  % CRRSS compartment voltage matrix

Gmc = Gmc0.*Model0I+Gmc1.*Model1I+Gmc3.*Model3I; 
Vc = Vc0.*Model0I+Vc1.*Model1I+Vc3.*Model3I;

% Calculate theoretical propagation velocities of the activation pulse
% (necessarily same timesteps as gate parameters)
% column dimension = space, row dimension = time, page dimension =
% frequency
if ~SweepActI
if PhaseVPlot % Phase velocity of activation pulse
omegaPhaseV = 2.*pi.*reshape(fPhaseV,1,1,[]);  % Redefine the frequency vector as page vector
if logical(sum(omegaPhaseV==0))     % if there is a zero frequency
Vphase(:,it,omegaPhaseV==0)=((2.*dxI)./CmcI).*sqrt(Gmc.*GacI);
end
omegaPhaseV0 = omegaPhaseV(omegaPhaseV~=0);     % Delete 0 from the page vector omega
if ~isempty(omegaPhaseV0)
Vphase(:,it,omegaPhaseV~=0) = gmultiply(dxI,gdivide(sqrt(2).*omegaPhaseV0,...
    sqrt(gadd(-Gmc./GacI,gmultiply(1./GacI,sqrt(gadd(Gmc.^2,...
    gmultiply(omegaPhaseV0.^2,CmcI.^2)))))))); % Vphase = omega/k
end
end
if GroupVPlot % Group velocity of activation pulse
omegaGroupV = 2.*pi.*reshape(fGroupV,1,1,[]);   % Redefine the frequency vector as a page vector
if logical(sum(omegaGroupV==0))
Vgroup(:,it,omegaGroupV==0)=((2.*dxI)./CmcI).*sqrt(Gmc.*GacI);
end
omegaGroupV0 = omegaGroupV(omegaGroupV~=0);   % Delete 0 from the page vector omega
if ~isempty(omegaGroupV0)
Vgroup(:,it,omegaGroupV~=0) = gmultiply(dxI,(2.*sqrt(2)).*sqrt(gadd(-Gmc./GacI,gmultiply(1./GacI,...
    sqrt(gadd(gmultiply(CmcI.^2,omegaGroupV0.^2),Gmc.^2))))).*sqrt(gadd(gmultiply(CmcI.^2,...
    omegaGroupV0.^2),Gmc.^2))./gmultiply(CmcI.^2./GacI,omegaGroupV0)); % Vgroup = domega/dk
end
end
if ActivationTable % Exact speed in activation table is based on group velocity on fExact
omegaExact = 2.*pi.*reshape(fExact,1,1,[]); % Redefine the frequency vector as a page vector  
if logical(sum(omegaExact==0))
Vexact(:,it,omegaExact==0)=((2.*dxI)./CmcI).*sqrt(Gmc.*GacI);
end
omegaExact0 = omegaExact(omegaExact~=0);    % Delete 0 from the page vector omega
if ~isempty(omegaExact0)
Vexact(:,it,omegaExact~=0) = gmultiply(dxI,(2.*sqrt(2)).*sqrt(gadd(-Gmc./GacI,gmultiply(1./GacI,...
    sqrt(gadd(gmultiply(CmcI.^2,omegaExact0.^2),Gmc.^2))))).*sqrt(gadd(gmultiply(CmcI.^2,...
    omegaExact0.^2),Gmc.^2))./gmultiply(CmcI.^2./GacI,omegaExact0)); % Vexact = domega/dk
end
end
end

if logical(sum(isnan(Gmc)|isinf(Gmc))) % Stop iteration if a singularity occurs (to simplify debugging)
if ~SweepActI
if DISPLAY
disp(' ');
disp(['Membrane conductance vector became singular at ' num2str(progress) '%']);
disp(['Iteration is stopped after ' num2str(it) ' iteration(s)']);
end
end
NeuronActivated = 1;
SingularityPass = 1;
break
end

GmcMAT = diag(Gmc);
if BoundaryConditions == 1 || BoundaryConditions == 2 % Keep voltage clamped
GmcMAT(1,1) = 0; GmcMAT(S,S) = 0; Vc(1) = 0; Vc(S) = 0;
end

A = A1+GmcMAT;
Vold = V;                           %Store Vold = V(t)
b = B*V+C*AF+Vc;

% To update the voltage, the Crank-Nicholson equation is solved:
% 7.1.1.2. if there is no model 2 (FH-model), model 4 (SE-model) or model 5 (SRB-model) 
% used, the Crank-Nicholson voltage step equation is linear and is solved
% by direct inverse:
% AV(t+dt/2) = BV(t)+CE(t+dt/2)+Vc

if sum(Model2l|Model4l|Model5l) == 0 % i.e. if no model 2,model4 or model 5 is used. 
V = A\b;                %V(t) -> V(t+dt/2)
else 
    
% 7.1.1.3. for FH, SE and SRB-model, the Crank-Nicholson equation becomes:
% A1*V(t+dt/2)+I(t+dt/2)= B*V(t)+C*E(t+dt/2) := b with I(t+dt/2) non-linear in
% V(t+dt/2) -> non-linear equation in V -> needs to be solved iteratively

%Model 2: FH-model
Ina2 = @(V) ((PnaI.*F^2)./(R.*(TempI+273))).*m.^2.*h.*V.*(cNaoI-cNaiI.*exp(10^(-3).*V.*(F./(R.*(TempI+273)))))./(1-exp(10^(-3).*V.*(F./(R.*(TempI+273)))));
Ik2 = @(V) ((PkI.*F^2)./(R.*(TempI+273))).*n.^2.*V.*(cKoI-cKiI.*exp(10^(-3).*V.*(F./(R.*(TempI+273)))))./(1-exp(10^(-3).*V.*(F./(R.*(TempI+273)))));
Ip2 = @(V) ((PpI.*F^2)./(R.*(TempI+273))).*p.^2.*V.*(cNaoI-cNaiI.*exp(10^(-3).*V.*(F./(R.*(TempI+273)))))./(1-exp(10^(-3).*V.*(F./(R.*(TempI+273)))));
Il2 = @(V) GlI.*(V-VlI);
I2 = @(V) pi.*dI.*dxI.*(Ina2(V)+Ik2(V)+Ip2(V)+Il2(V));

%Model 4: SE-model
Ina4 = @(V) ((PnaI.*F^2)./(R.*(TempI+273))).*m.^3.*h.*V.*(cNaoI-cNaiI.*exp(10^(-3).*V.*(F./(R.*(TempI+273)))))./(1-exp(10^(-3).*V.*(F./(R.*(TempI+273)))));
Ik4 =   @(V) ((PkI.*F^2)./(R.*(TempI+273))).*n.^2.*V.*(cKoI-cKiI.*exp(10^(-3).*V.*(F./(R.*(TempI+273)))))./(1-exp(10^(-3).*V.*(F./(R.*(TempI+273)))));
Il4 = @(V) GlI.*(V-VlI);
I4 = @(V) pi.*dI.*dxI.*(Ina4(V)+Ik4(V)+Il4(V));

%Model 5: SRB-model
Ina5 = @(V) ((PnaI.*F^2)./(R.*(TempI+273))).*m.^3.*h.*V.*(cNaoI-cNaiI.*exp(10^(-3).*V.*(F./(R.*(TempI+273)))))./(1-exp(10^(-3).*V.*(F./(R.*(TempI+273)))));
Ikfast5 = @(V) GkI.*n.^4.*(V-VkI);
Ikslow5 = @(V) GkslowI.*p.*(V-VkI);
Il5 = @(V) GlI.*(V-VlI);
I5 = @(V) pi.*dI.*dxI.*(Ina5(V)+Ikslow5(V)+Ikfast5(V)+Il5(V));    


I = @(V) I2(V).*Model2I+I4(V).*Model4I+I5(V).*Model5I;    % general I vector

% A1*V+I(V) = b has to be solved (with V = V(t+dt/2)) where model 2,4 or 5
% is applied. The equation is solved in general, together with the
% compartments where model 0,1 or 3 are used.
% Because this equation is now non-linear due to I(V), in the compartments of
% model 2,4 or 5, we need to solve it iteratively

% The non-linear cranck-nicholson (NLCN) equation becomes:
% (A1*V+diag(Gmc)*V+I(V)=) A*V+I(V) = b (=B*V+C*El+Vc) 
% -> A*V+I(V)-b = 0 -> F(V) = 0
if OrderOfSolution == 2
FUN = @(V) A*V+I(V)-b;
if DynTol == 0
   utol = tol;
else
   utol = DynTol*min(abs(b));
end
options = optimoptions('fsolve','TolFun',utol,'Display','none');
V = fsolve(FUN,Vold,options);
else if OrderOfSolution == 1
V = A\(b-I(Vold));
    end
end
end


% 7.1.2. Update voltage t -> t+dt
% V(t+dt) = 2*V(t+dt/2)-V(t)
V = 2.*V-Vold;

if logical(sum(isnan(V)|isinf(V)))      %Stop iterations at singularity
if ~SweepActI
if DISPLAY
disp(' ');
disp(['Membrane potential became singular at ' num2str(progress) '%']);
disp(['Iteration is stopped after ' num2str(it) ' iteration(s)']);
end
end
NeuronActivated = 1;
SingularityPass = 1;
break
end
if ~SweepActI
if PlotTimeV
recorderV(:,it+1) = V(ceil(S.*PlotVt));
end
if PlotSpaceV 
if logical(sum(unique(ceil((Nt+1).*PlotVx)) == it+1))
recorderVx(:,(unique(ceil((Nt+1).*PlotVx)) == it+1)) = V;
end
end
end

if VtotPlot || NeuronActivation || SweepActI
recorderVtot(:,it+1) = V;
end


% 7.2 Update gate-parameters
% Gate parameters are updated t+dt/2 -> t+3dt/2 
if sum(Model~=0) ~= 0            % No gate parameters for passive model
m = (amI(V)+(1/dt).*m-((amI(V)+bmI(V))./2).*m)./((1/dt).*ones(S,1)+(amI(V)+bmI(V))./2);
h = (ahI(V)+(1/dt).*h-((ahI(V)+bhI(V))./2).*h)./((1/dt).*ones(S,1)+(ahI(V)+bhI(V))./2);
m(m>1) = 1; m(m<0)=0; h(h>1) = 1; h(h<0)=0;
if sum(Model~=3) ~= 0
n = (anI(V)+(1/dt).*n-((anI(V)+bnI(V))./2).*n)./((1/dt).*ones(S,1)+(anI(V)+bnI(V))./2);
n(n>1) = 1; n(n<0) = 0;
end 
end
if sum((Model==2)|(Model==5))~=0
p = (apI(V)+(1/dt).*p-((apI(V)+bpI(V))./2).*p)./((1/dt).*ones(S,1)+(apI(V)+bpI(V))./2);
p(p>1)=1;p(p<0)=0;
end

%Stop iterations at singularity:
if logical(sum(isnan(n)|isinf(n)|isnan(h)|isinf(h)|isnan(p)|isinf(p)|isnan(m)|isinf(m))) 
if ~SweepActI
if DISPLAY
disp(' ');
disp(['Gate parameter became singular at ' num2str(progress) '%']);
disp(['Iteration is stopped after ' num2str(it) ' iteration(s)']);
end
end
NeuronActivated = 1;
SingularityPass = 1;
break
end

if ~SweepActI
if sum(Model~=0) ~= 0
if PlotTimem
recorderm(:,it+1) = m(ceil(S.*Plotmt));
end

if PlotSpacem
if logical(sum(unique(ceil((Nt+1).*Plotmx)) == it+1))
recordermx(:,(unique(ceil((Nt+1).*Plotmx)) == it+1)) = m;
end
end

if mtotPlot 
    recordermtot(:,it+1) = m;
end

if sum(Model~=3)~=0
if PlotTimen 
recordern(:,it+1) = n(ceil(S.*Plotnt));    
end

if PlotSpacen
if logical(sum(unique(ceil((Nt+1).*Plotnx)) == it+1))
recordernx(:,(unique(ceil((Nt+1).*Plotnx)) == it+1)) = n;
end
end

if ntotPlot 
    recorderntot(:,it+1) = n;
end
end
if PlotTimeh 
recorderh(:,it+1) = h(ceil(S.*Plotht));
end

if PlotSpaceh
if logical(sum(unique(ceil((Nt+1).*Plothx)) == it+1))
recorderhx(:,(unique(ceil((Nt+1).*Plothx)) == it+1)) = h;
end
end

if htotPlot
    recorderhtot(:,it+1) = h;
end

if sum((Model==2)|(Model==5)) ~= 0
    
if PlotTimep
recorderp(:,it+1) = p(ceil(S.*Plotpt));
end

if PlotSpacep
if logical(sum(unique(ceil((Nt+1).*Plotpx)) == it+1))
recorderpx(:,(unique(ceil((Nt+1).*Plotpx)) == it+1)) = p;
end
end
  
if ptotPlot 
    recorderptot(:,it+1) = p;
end
end

end
end
end

if ~SingularityPass
TIME2 = toc;
% Velocities at t = dt*Nt+dt/2 have yet to be defined after the loop
% -> first: define the conductivities (depending on gate parameters at t =
% Nt*dt+dt/2)
Gmc1 = pi.*dI.*dxI.*(GnaI.*m.^3.*h+GkI.*n.^4+GlI);  %HH compartment membrane conductance at t = T+dt/2
Gmc3 = pi.*dI.*dxI.*(GnaI.*m.^2.*h+GlI);  %CRRSS compartment membrane conductance at t = T+dt/2
Gmc = Gmc0.*Model0I+Gmc1.*Model1I+Gmc3.*Model3I; 

% -> Secondly: Calculate theoretical propagation velocities of the activation pulse
% at T=Nt*dt+dt/2
if ~SweepActI
if PhaseVPlot % Phase velocity of activation pulse at t = T+dt/2
if logical(sum(omegaPhaseV==0))    
Vphase(:,end,omegaPhaseV==0)=((2.*dxI)./CmcI).*sqrt(Gmc.*GacI);
end
if ~isempty(omegaPhaseV0)
Vphase(:,end,omegaPhaseV~=0) = gmultiply(dxI,gdivide(sqrt(2).*omegaPhaseV0,...
    sqrt(gadd(-Gmc./GacI,gmultiply(1./GacI,sqrt(gadd(Gmc.^2,...
    gmultiply(omegaPhaseV0.^2,CmcI.^2)))))))); % Vphase = omega/k
end
end
if GroupVPlot % Group velocity of activation pulse at t = T+dt/2
if logical(sum(omegaGroupV==0))
Vgroup(:,end,omegaGroupV==0)=((2.*dxI)./CmcI).*sqrt(Gmc.*GacI);
end
if ~isempty(omegaGroupV0)
Vgroup(:,end,omegaGroupV~=0) = gmultiply(dxI,(2.*sqrt(2)).*sqrt(gadd(-Gmc./GacI,gmultiply(1./GacI,...
    sqrt(gadd(gmultiply(CmcI.^2,omegaGroupV0.^2),Gmc.^2))))).*sqrt(gadd(gmultiply(CmcI.^2,...
    omegaGroupV0.^2),Gmc.^2))./gmultiply(CmcI.^2./GacI,omegaGroupV0)); % Vgroup = domega/dk
end
end
if ActivationTable % Exact speed in activation table is based on group velocity on fExact (defined here at t=T+dt/2)
if logical(sum(omegaExact==0))
Vexact(:,end,omegaExact==0)=((2.*dxI)./CmcI).*sqrt(Gmc.*GacI);
end
if ~isempty(omegaExact0)
Vexact(:,end,omegaExact~=0) = gmultiply(dxI,(2.*sqrt(2)).*sqrt(gadd(-Gmc./GacI,gmultiply(1./GacI,...
    sqrt(gadd(gmultiply(CmcI.^2,omegaExact0.^2),Gmc.^2))))).*sqrt(gadd(gmultiply(CmcI.^2,...
    omegaExact0.^2),Gmc.^2))./gmultiply(CmcI.^2./GacI,omegaExact0)); % Vexact = domega/dk
end
end
end

% 8. Determination of neuronactivation
if ~SweepActI
if DISPLAY
disp(' ');
disp('3. Determining neuronactivation');
end
end
if NeuronActivation || SweepActI                   % We determine from recorderVtot if the neuron is activated 
recorderVtotRef = recorderVtot-repmat(VrI,1,Nt+1);   % We refer all voltages with respect to the rest potential 
Vactive = recorderVtotRef(~Model0I,:);  % To determine activation, we only use the active regions (mV)
AS = sum(~Model0I);                      % Number of active segments (A.S.)
xAct = x(~Model0I);                    % Associated active positions (m)
distb = xAct-circshift(xAct,1,1);           % backward distance between active segments
distf = circshift(xAct,-1,1)-xAct;       % forward distance between active segments
% -> Vactive is referred to the rest potential, to simplify comparison with Vtres.

% Convention: a neuron is activated at a location xAct and time tAct if 2 conditions hold:
% A) The neuron is depolarised by Vtres at xAct at time tAct. (Local
% activation)
% B) The activation pulse propagates to other segments (i.e. x ~= xAct at
% times t > tAct, such that abs(x-xAct)/(t-tAct) < c). (Propagation)

% 8.1. Determine local activation peaks
% We define the local maxima in time at each segment as potential
% activation pulses. -> for local activation, these maxima have to be
% higher than Vtres.
pkscell = cell(AS,1); timescell = cell(AS,1); % pkscell = all peaks (local maxima) at each active segment. timescell = corresponding times
propf = cell(AS,1); propb = cell(AS,1);     % propf = propagating forward, propb = propagating backward (for later use, see 8.2)
pksAPdurcell = cell(AS,1); timesAPdurcell = cell(AS,1);  % These cells take into account a minimal AP duration
for i = 1:AS
[pks,locs] = findpeaks(Vactive(i,:));   % find peaks and corresponding location numbers of Vactive.
times = timeflow(locs(pks>Vtres)); pks = pks(pks>Vtres);  % Only peaks > Vtres are local activation peaks. We also determine the associated times.
pkscell{i} = pks; timescell{i} = times; % All values are stored in cells: f.i. pkscell{j} = rowmatrix with all local activation peaks of active segment j.
end
SynapseActivated = 1;              % Check for activation of the synapse
if isempty(pkscell{AS})
    SynapseActivated = 0;
end

for i = 1:AS
    timesAPduri = timescell{i};   % APdur: take into account minimal AP duration (minAPdur)
    pksAPduri = pkscell{i}; 
    temptimesAPduri = timesAPduri;
    if ~isempty(temptimesAPduri)
    indAPduri = [1]; %#ok<NBRAK>
    else 
    indAPduri = [];
    end
    while ~isempty(temptimesAPduri)
    indNextAP = find(diff(temptimesAPduri)>=minAPdur,1)+indAPduri(end);
    % If indNextAP is empty, it will have no effect on indAPduri
    % andtimesAPduri will be empty. 
    indAPduri = horzcat(indAPduri,indNextAP); %#ok<AGROW> 
    temptimesAPduri = temptimesAPduri(indNextAP:end);
    end
    pksAPdurcell{i} = pksAPduri(indAPduri);
    timesAPdurcell{i} = timesAPduri(indAPduri);    
end

nTermAPs = max([length(pksAPdurcell{1}) length(pksAPdurcell{AS})]);     % Number of APs at the most active terminal
TerminalActivated = double(nTermAPs >= ThreshnAPs);              % Terminal activated w.r.t. ThreshnAPs     
% 8.2. Determine propagation numbers

% We have determined the local activation peaks and times. Now we
% determine if these peaks propagate (defined by the Conditions matrix): if
% they do, the neuron is activated. To do this two cells are defined:
% propf (propagate forward) and propb (propagate backward). They have the
% same structure as pkscell/timescell and each entry is an integer that
% indicates the amount of active segments that the peak has propagated
% forward/backward. 
if isempty(timescell{1})
propb{1} = [];  % This if-loop is introduced, because Matlab distinguishes between 1x0 matrices and empty [] (0x0) matrices
else
propb{1} = zeros(1,length(timescell{1}));  % At the edge, no backward propagation is possible
end
for i = 2:AS
% Propating condition: 0< distb(i)/(Backt-Currentt) < c --> if this
% condition is satisfied, the peak at the previous active segment can be
% reached by the peak of the current active segment.
Currentpropb = (propb{i-1}+1)';         % Increment previous propagating numbers (the peak propagated over +1 active segment)
Backt = timescell{i-1}; Currentt = timescell{i}; % The peak times at the previous (backt) and current (currentt) active segment
if isempty(Backt)||isempty(Currentt)
vb = [];
else
vb = distb(i)./(gsubtract(Backt',Currentt));
end
% Propagating speed: v_ij is the speed with which the current peak j would
% propagate to backward (at previous active segment) peak i.
Reachableb = double(0<vb&vb<maxc);   %(empty if any of the matrices Backt/Curentt is empty)
% if Reachable_ij = 1, peak i at the previous active segment, is caused by the propagation 
% of peak j at the current active segment
if isempty(Currentt)
propb{i} = [];
else if isempty(Backt)
propb{i} = zeros(1,length(Currentt));
    else
propb{i} = max(gmultiply(Reachableb,Currentpropb),[],1); 
    end
end
end

if isempty(timescell{AS})
propf{AS} = [];  % This if-loop is introduced, because Matlab distinguishes between 1x0 matrices and empty [] (0x0) matrices
else
propf{AS} = zeros(1,length(timescell{AS}));  % At the edge, no forward propagation is possible
end
for i = AS-1:-1:1
% Propating condition: 0< distf(i)/(Forwardt-Currentt) < c --> if this
% condition is satisfied, the peak at the next active segment can be
% reached by the peak the the current active segment.
Currentpropf = (propf{i+1}+1)';         % Increment previous propagating numbers (the peak propagated over +1 active segment)
Forwardt = timescell{i+1}; Currentt = timescell{i}; % The peak times at the previous (backt) and current (currentt) active segment
if isempty(Forwardt)||isempty(Currentt)
vf = [];
else
vf = distf(i)./(gsubtract(Forwardt',Currentt));
end
% Propagating speed: v_ij is the speed with which the current peak j would
% propagate to forward (at next active segment) peak i.
Reachablef = double(0<vf&vf<maxc);   %(empty if any of the matrices Nextt/Curentt is empty)
% if Reachable_ij = 1, peak i at the next active segment, is caused by the propagation 
% of peak j at the current active segment
if isempty(Currentt)
propf{i} = [];
else if isempty(Forwardt)
propf{i} = zeros(1,length(Currentt));
    else
propf{i} = max(gmultiply(Reachablef,Currentpropf),[],1); 
    end
end
end

% 8.3. Activated positions and times
% Positions and times of activation are determined by matching propf and
% propb with Conditions.
NeuronActivated = 0;                % Bool: 1 if the neuron is activated, according to one of the conditions
TableTime = []; TablePos = []; TablePropf = []; TablePropb = [];
for i = 1:length(Conditions(:,1))     % Iterate over all possible conditions
geFC = @ (x) x >= Conditions(i,1); geBC = @(x) x >= Conditions(i,2); % Function handles for >= forward condition (geFC) and >= backward conditions (geBC)
ForwardActivated = cellfun(geFC,propf,'UniformOutput',false);   % Logical cell: 1 where the neuron is activated forward according to condition i
BackwardActivated = cellfun(geBC, propb,'UniformOutput',false); % Logical cell: 1 where the neuron is actived backward according to condition i
Activated = cellfun(@and,ForwardActivated,BackwardActivated,'UniformOutput',false); % Logical cell: the neuron is activated (1) if it is activated forward and backward according to condition i
NeuronActivated = NeuronActivated||logical(sum(cellfun(@sum,Activated))); % The neuron is activated, if NeuronaActivated == 1 for one condition.
if ActivationTable % For the activation table, we need all activation times (TableTime) and the corresponding conditions (TablePos)
% 1) Tabulate activation times
TableTimeCell =  cellfun(@(x,y) y(x),Activated,timescell,'UniformOutput',false); % time (s)
TableTimeCell(cellfun(@isempty,TableTimeCell)) = []; % Omit empty cells -> because empty cells may have different column size -> can't be concatenated in future matlab version
TableTimeCell = cellfun(@transpose,TableTimeCell,'UniformOutput',false); % Tranpose to apply vertcat later
TableTime = vertcat(TableTime,vertcat(TableTimeCell{:})); %#ok<AGROW> All activated times are tabulated (s)
% 2) Tabulate activation positions
N2 = cellfun(@sum,Activated);   % Holds the number of time the active position has to be replicated           
N2o = N2(logical(N2));          % -> the same indexing method as applied previously is used (-> subscript 2)
if ~isempty(N2o)
xActo = xAct(logical(N2));      % -> to do this, we omit the 0 values of N2 in N2 and in Act -> N20 and xActo
index2 = zeros(sum(N2o),1);       % 1) An index vector with the same dimension as TablePos is initialised
index2([1; cumsum(N2o(1:end-1))+1]) = 1; % 2) 1 is declared at each index where a new sequence of the same TablePos value will start
index2 = cumsum(index2);
TablePos = vertcat(TablePos,xActo(index2));              %#ok<AGROW> % 3) TablePos positions are defined (m)            
end
% 3) Tabulate forward propagation  (we will use this to calculate MSOA)
TablePropfCell= cellfun(@(x,y) y(x),Activated,propf,'UniformOutput',false); % forward propagations
TablePropfCell(cellfun(@isempty,TablePropfCell)) = []; % Omit empty cells
TablePropfCell = cellfun(@transpose,TablePropfCell,'UniformOutput',false); % Transpose to apply vertcat later 
TablePropf = vertcat(TablePropf,vertcat(TablePropfCell{:})); %#ok<AGROW> All forward propagations are tabulated.
% 4) Tabulate backward propagation
TablePropbCell= cellfun(@(x,y) y(x),Activated,propb,'UniformOutput',false); % backward propagations
TablePropbCell(cellfun(@isempty,TablePropbCell)) = []; % Omit empty cells
TablePropbCell = cellfun(@transpose,TablePropbCell,'UniformOutput',false); % Transpose to apply vertcat later 
TablePropb = vertcat(TablePropb,vertcat(TablePropbCell{:})); %#ok<AGROW> All backward propagations are tabulated.
end
end

% 8.4. Calculate overview table of activation times and locations
% example of table: 
%            Activated   Synaps activated    nr.      pos.(mm)         time (ms)   Mean velocity of activationpulse (m/s) 
% Neuron 1    yes                 yes         1.       5                200                      100   
%                                             2.       8                150                      50
% Neuron 2    no                  no          -        -                -                         -
  
if ActivationTable && ~SweepActI
% 8.4.1. Check for causality
% The activation times and positions in TablePos and TableTime may be
% causally related or even equal, so we need to check for causality first:
% position i and j are causally related, if
% abs((TablePos(i)-TablePos(j))/(TableTime(i)-TableTime(j)) < maxc ->
% earlier time has to be deleted

% 1   1   We circshift the second tablepos/tabletime to check for causality.   
% 2   2   To determine the minimal amount of iterations it, we use:
% 3   3   even: 2it=Ntable -> it = Ntable/2 = floor(Ntable/2)
% 4   4   odd: largest it: it+1<Ntable+1-it -> 2it < Ntable -> it < Ntable/2 -> it = floor(Ntable/2)
% 5   5
Ntable = length(TablePos); debugTableTime = TableTime; debugTablePos = TablePos; % The debug values are defined for debugging only
delete = []; debugTablePropf = TablePropf; debugTablePropb = TablePropb;
for i = 1:floor(Ntable/2)
   Causal = abs((TablePos-circshift(TablePos,-i,1))./(TableTime-circshift(TableTime,-i,1)));
   % Remark: equal positions/times will result in 0/0 = NaN and will pass
   % through the check -> we have to filter for them later
   Check = (Causal<maxc);       % Check(k) = 1 iff position k and k+i (circularily; i.e. :k+i-1 mod Ntable +1)  are causally related 
   CausalLoc1 = find(Check); % CausalLoc1 = k
   CausalLoc2 = mod(CausalLoc1+i-1,Ntable)+1; % CausalLoc2 = (k+i-1) mod Ntable +1
   causal12 = double(TableTime(CausalLoc1) < TableTime(CausalLoc2)); % time 1 caused time 2 (locations)
   causal21 = double(TableTime(CausalLoc1) > TableTime(CausalLoc2)); % time 2 caused time 1 (locations)
   delete = vertcat(delete,causal12.*CausalLoc2+causal21.*CausalLoc1);  %#ok<AGROW> % indexes to be deleted from TableTime and TablePos
end
  TableTime(delete) = []; TablePos(delete) = []; TablePropf(delete) = []; TablePropb(delete) = []; %#ok<SAGROW>
  [TablePosandTime,ia,~] = unique([TablePos TableTime],'rows'); % Filter the equal rows that passed through the "Check" due to NaN
  TablePropf = TablePropf(ia); TablePropb = TablePropb(ia);
  if size(TablePosandTime,1) ~= 0
  TablePos = TablePosandTime(:,1); TableTime = TablePosandTime(:,2);
  else 
  TablePos = []; TableTime = [];
  end

  % 8.4.2. Calculate MSOA (= mean spead of activation pulse)
  % 8.4.2.1. Simulation of the MSOA
  % The mean speed is calculated by calculating the mean speed of the
  % forward propagating activation pulse (MSOAf) and the mean speed of the
  % backward propagating pulse (MSOAb) -> then MSOA is the weighted average
  % of MSOAf and MSOAb, with the propagating lengths as weights.
  if ~isempty(TablePos)
  [~,IndexTablePos] = ismembertol(TablePos,xAct,1e-10);          % Vector with the indices of actived segments
  DistPropf = (xAct(IndexTablePos+TablePropf)-xAct(IndexTablePos));  % Vector with propagated distances (forward)
  DistPropb = -(xAct(IndexTablePos-TablePropb)-xAct(IndexTablePos)); % Vector with propagated distances (backward)
  TimePropf = cellfun(@minus,timescell(IndexTablePos+TablePropf),num2cell(TableTime),'UniformOutput',false); % Cell with elapsed times for propagating forward
  TimePropb = cellfun(@minus,timescell(IndexTablePos-TablePropb),num2cell(TableTime),'UniformOutput',false); % Cell with elapsed times for propagating backward
  Vf = cellfun(@rdivide,num2cell(DistPropf),TimePropf,'UniformOutput',false); % Cell with all forward velocities
  InfDeleteVf = cellfun(@not,cellfun(@isinf,Vf,'UniformOutput',false),'UniformOutput',false); % Mask to delete all inf values from the Vf cell
  Vf = cellfun(@(x,y) y(x),InfDeleteVf,Vf,'UniformOutput',false); % Delete all Inf values from Vf
  Vb = cellfun(@rdivide,num2cell(DistPropb),TimePropb,'UniformOutput',false); % Cell with all backward velocities
  InfDeleteVb = cellfun(@not,cellfun(@isinf,Vb,'UniformOutput',false),'UniformOutput',false); % Mask to delete all inf values from the Vb cell
  Vb = cellfun(@(x,y) y(x),InfDeleteVb,Vb,'UniformOutput',false); % Delete all Inf values from Vb
  MaskVf = cellfun(@(x) x > 0 & x < maxc,Vf,'UniformOutput',false); % Maskcell for forward velocity 
  MaskVb = cellfun(@(x) x > 0 & x < maxc,Vb,'UniformOutput',false); % Maskcell for backward velocity
  mult = @(x,y) x.*y; % multiplication handle (same as @times, but this doesn't work anymore because 'times' is already preallocated for the activation times)
  MSOAsf = cellfun(mult,MaskVf,Vf,'UniformOutput',false); % forward MSOAs cell -> 's' because there may be more than one 
  MSOAsb = cellfun(mult,MaskVb,Vb,'UniformOutput',false); % backward MSOAs cell
  MSOAf = cellfun(@max,MSOAsf); % forward MSOA vector -> here we only consider the activation pulse with the highest (physical, i.e. < maxc) speed, because this is the most important pulse 
  MSOAb = cellfun(@max,MSOAsb); % backward MSOA vector
%   MSOA = (MSOAf.*DistPropf+MSOAb.*DistPropb)./(DistPropf+DistPropb); % Calculate MSOA
%   MSOA(isnan(MSOAf)) = MSOAb; MSOA(isnan(MSOAb)) = MSOAf; % No propagation in one direction results in NaN
  MSOAfcell = num2cell(MSOAf); MSOAfcell(isnan(MSOAf)) = cellstr('No forward prop.');
  MSOAbcell = num2cell(MSOAb); MSOAbcell(isnan(MSOAb)) = cellstr('No backward prop.');
  end

% 8.4.2.2. Theoretical assesment of the MSOA
NPTs = length(TablePos);        % number of paths of activation pulses
NFs = length(fExact);           % number of frequencies
if NPTs ~= 0 && NFs ~=0
[~,Indexf] = ismembertol(xAct(IndexTablePos+TablePropf),x,1e-10); % NPTs x 1 vectors of the indices of the end of the path in forward
[~,Index0] = ismembertol(xAct(IndexTablePos),x,1e-10);            % direction, initial location, backward direction
[~,Indexb] = ismembertol(xAct(IndexTablePos-TablePropb),x,1e-10);
PTLf = Indexf-Index0;            % NPTs x 1 vector: the forward pathlength 
PTLb = Index0-Indexb;            % NPTs x 1 vector: the backward pathlength
xPathf = Czeros(PTLf+1);              % We store the positions along the path of the activationpulse
xPathb = Czeros(PTLb+1);          % in a cell: rows -> different activationpulse, columns -> different frequencies
tPathf = repmat(Czeros(PTLf+1),1,NFs);  % Also the times are stored.
tPathb = repmat(Czeros(PTLb+1),1,NFs);   % -> NPTs x NFs cells and NPTs x 1 cells

% Initialisation of locations and times as NPTs x NFs matrices (1 matrix
% per iteration step)
Locf = TablePos;  % Initial times and positions correspond to TablePos and TableTime 
Locb = TablePos;
Timf = repmat(TableTime,1,NFs);
Timb = repmat(TableTime,1,NFs);
xPathf = DeclareCell(xPathf,Locf,1);    % Declare positions and times to the cells
xPathb = DeclareCell(xPathb,Locb,1);
tPathf = DeclareCell(tPathf,Timf,1);
tPathb = DeclareCell(tPathb,Timb,1);

for iINTf = 2:max(PTLf)+1     % Integration of Vexact in forward direction
ExtinctMaskf = (iINTf <= PTLf+1);  % mask is 0 when path is extinct
MaskedLocf = Locf(ExtinctMaskf);
MaskedTimf = Timf(repmat(ExtinctMaskf,1,NFs));
[~,PositionIndex] = ismembertol(MaskedLocf,x,1e-10);    % Determine positions indices
VInterp = permute(Vexact(PositionIndex,:,:),[2 1 3]); % Setting velocity and times up for interpolation with interp1
TimfInterp = permute(MaskedTimf,[3 1 2]);
tInterp = timeflow'+dt/2;
for i=1:NPTs
    for j=1:NFs
Vprop(i,j) = interp1(tInterp,VInterp(:,i,j),TimfInterp(1,i,j)); %#ok<SAGROW>
    end
end % Calculate prop velocity (NPTs x NFs matrix)
Propdtf = gdivide(dxfI(PositionIndex),Vprop); % NPTs x NFs matrix containing incremental dt
Timf(repmat(ExtinctMaskf,1,NFs)) = MaskedTimf+Propdtf; % Increment times
Locf(ExtinctMaskf) = MaskedLocf+dxfI(PositionIndex);  % Increment locations
xPathf(ExtinctMaskf) = DeclareCell(xPathf(ExtinctMaskf),Locf(ExtinctMaskf),iINTf);    % Declare positions and times to the cells
tPathf(repmat(ExtinctMaskf,1,NFs)) = DeclareCell(tPathf(repmat(ExtinctMaskf,1,NFs)),Timf(repmat(ExtinctMaskf,1,NFs)),iINTf);
end

for iINTb = 2:max(PTLb)+1      % Integration of Vexact in backward direction
ExtinctMaskb = (iINTb <= PTLb+1);  % mask is 0 when path is extinct
MaskedLocb = Locb(ExtinctMaskb);
MaskedTimb = Timb(repmat(ExtinctMaskb,1,NFs));
[~,PositionIndex] = ismembertol(MaskedLocb,x,1e-10);    % Determine positions indices
VInterp = permute(Vexact(PositionIndex,:,:),[2 1 3]); % Setting velocity and times up for interpolation with interp1
TimbInterp = permute(MaskedTimb,[3 1 2]);
tInterp = timeflow'+dt/2;
for i=1:NPTs
    for j=1:NFs
Vprop(i,j) = interp1(tInterp,VInterp(:,i,j),TimbInterp(1,i,j));
    end
end
% Calculate prop velocity (NPTs x NFs matrix)
Propdtb = gdivide(dxbI(PositionIndex),Vprop); % NPTs x NFs matrix containing incremental dt
Timb(repmat(ExtinctMaskb,1,NFs)) = MaskedTimb+Propdtb; % Increment times
Locb(ExtinctMaskb) = MaskedLocb-dxbI(PositionIndex);  % Increment locations
xPathb(ExtinctMaskb) = DeclareCell(xPathb(ExtinctMaskb),Locb(ExtinctMaskb),iINTb);    % Declare positions and times to the cells
tPathb(repmat(ExtinctMaskb,1,NFs)) = DeclareCell(tPathb(repmat(ExtinctMaskb,1,NFs)),Timb(repmat(ExtinctMaskb,1,NFs)),iINTb);
end

% When end of simulationdomain is reached we get nan, delete these:
NaNmaskf = cellfun(@isnan,tPathf,'UniformOutput',false); % 1 if NaN
NaNmaskb = cellfun(@isnan,tPathb,'UniformOutput',false);
tPathfC = cellfun(@(x,y) x(~y),tPathf,NaNmaskf,'UniformOutput',false); 
tPathbC = cellfun(@(x,y) x(~y),tPathb,NaNmaskb,'UniformOutput',false);
tPathfLength = num2cell(cellfun(@length,tPathfC));
tPathbLength = num2cell(cellfun(@length,tPathbC));
xPathfC = cellfun(@(x,y) x(1:y),repmat(xPathf,1,NFs),tPathfLength,'UniformOutput',false);
xPathbC = cellfun(@(x,y) x(1:y),repmat(xPathb,1,NFs),tPathbLength,'UniformOutput',false);
% NPTs x NFs matrices of theoretical approximated MSOAf (calculated by
% dividing appdx by appdt
appdxf= cellfun(@(x) x(end),xPathfC)-cellfun(@(x) x(1),xPathfC); 
appdxb = cellfun(@(x) x(1),xPathbC)-cellfun(@(x) x(end),xPathbC);
appdtf = cellfun(@(t) t(end),tPathfC)-cellfun(@(t) t(1),tPathfC);
appdtb = cellfun(@(t) t(end),tPathbC)-cellfun(@(t) t(1),tPathbC);
appMSOAf = gdivide(appdxf,appdtf);
appMSOAb = gdivide(appdxb,appdtb);
end

% 8.4.3. Plot the activation table
uifigure = figure();
nr = size(TablePosandTime,1);
TablePosandTime = 1000.*TablePosandTime; % position in time in (mm) and (ms) respectively.
TablePosandTimeCell = num2cell(TablePosandTime);

% Out-comment
% global PositionAPmm;
% if isempty(TablePosandTime)
% PositionAPmm = nan;
% else
% PositionAPmm = TablePosandTime(1);
% end 

ColumnName = {'Activated' 'Synapse activated' 'nr.' 'pos. (mm)' 'time (ms)' 'MSOAf (m/s)' 'MSOAb (m/s)'};
appMSOAfName = cellstr([repmat('appMSOAf (',NFs,1) num2str(fExact') repmat('Hz)',NFs,1)])';
appMSOAbName = cellstr([repmat('appMSOAb (',NFs,1) num2str(fExact') repmat('Hz)',NFs,1)])';
ColumnName = horzcat(ColumnName,appMSOAfName,appMSOAbName); ColumnName(cellfun(@isempty,ColumnName)) = []; %#ok<AGROW>
RowName = vertcat(cellstr(TableName),cell(nr-1,1));
if NeuronActivated
nrs = cellfun(@horzcat,cellstr(num2str((1:nr)')),repmat({'.'},nr,1),'UniformOutput',false); % create a cell: {'1.'; '2.'; ...}
Activated = 'yes';
if NFs ~= 0
appMSOAfcell = num2cell(appMSOAf); appMSOAbcell = num2cell(appMSOAb);
else
appMSOAfcell = cell(nr,1); appMSOAbcell = cell(nr,1);
end
else Activated = 'no';
nrs = '-'; TablePosandTimeCell = {('-') ('-')}; MSOAfcell = {'-'}; MSOAbcell = {'-'}; appMSOAfcell =  repmat({'-'},NFs,1); 
appMSOAbcell = repmat({'-'},NFs,1);
end
if SynapseActivated
    Synapse = 'yes'; else Synapse = 'no';
end
Synapse = vertcat(Synapse,cell(nr-1,1)); %#ok<AGROW>
Activated = vertcat(Activated,cell(nr-1,1)); %#ok<AGROW>
if NFs ~= 0
DATA = horzcat(Activated,Synapse,nrs,TablePosandTimeCell,MSOAfcell,MSOAbcell,appMSOAfcell,appMSOAbcell);
else
DATA = horzcat(Activated,Synapse,nrs,TablePosandTimeCell,MSOAfcell,MSOAbcell);
end
ActivationTABLE = uitable(uifigure,'Data',DATA,'RowName',RowName,...
    'ColumnName',ColumnName);   % create the activation table
ActivationTABLE.Position(3) = ActivationTABLE.Extent(3);  % Adjust the outer rectangle to match the table
ActivationTABLE.Position(4) = ActivationTABLE.Extent(4);
set(uifigure,'Position',[500 500 ActivationTABLE.Position(3)+50 ActivationTABLE.Position(4)+50]); % The width of the figure is adjusted to match the table
end
end 
end
if SweepActI
    switch ActivationDetectionMode
        case 1, SweepActivatedPar = NeuronActivated; 
        case 2, SweepActivatedPar = TerminalActivated;
    end
Ielecs(SweepActivatedPar+1) = IIelecs;
if SearchMode % SearchMode = 1
IIelecs = (Ielecs(1)+Ielecs(2))/2;
PrecisionCheck=((Ielecs(2)-Ielecs(1))<10^(floor(log10(Ielecs(2)))-ElecsPrecision));
else % SearchMode = 0
if ~any(~Ielecs)
SearchMode = 1;
IIelecs = (Ielecs(1)+Ielecs(2))/2;
PrecisionCheck=((Ielecs(2)-Ielecs(1))<10^(floor(log10(Ielecs(2)))-ElecsPrecision));
else
IIelecs = IIelecs*10^(3-2*find(Ielecs));
end
end
    switch SENNmodel
        case 1 
            modelstr = 'HH';
        case 2
            modelstr = 'FH';
        case 3
            modelstr = 'CRRSS';
        case 4
            modelstr = 'SE';
        case 5
            modelstr = 'SRB';
    end
Checkpoint=['Sweep-' num2str(LLx0s) '-' num2str(LLy0s) '-' num2str(LLz0s) '-'...  
    num2str(LLTp) '-' num2str(LLancat) '-' num2str(LLmonbi) '(' modelstr '-Conf' num2str(conf) '):' num2str(Ielecs)];
if DISPLAY
fprintf([Checkpoint '\n']);
end
end
flg = ~SweepActI;
flg2 = 1;
if flg, break; end
end
SweepActImat(Inx0s,Iny0s,Inz0s,InTps,Iancat,Imb) = IIelecs;
if SweepActI
NeuronActivated = 0;
TerminalActivated = 0;
end
 if flg, break; end
end
 if flg, break; end   
end
 if flg, break; end
end
 if flg, break; end
end
 if flg, break; end
end
 if flg, break; end
end

           
if ~SweepActI
% 9. Plotting the recorders
if DISPLAY
disp('4. Calculating plots');
end
t = timeflow;

% 9.1. Membrane-voltage
if reference
VrIt = repmat(VrI(ceil(S.*PlotVt)),1,Nt+1);
VrIx = repmat(VrI,1,length(unique(ceil((Nt+1).*PlotVx))));
VrItot = repmat(VrI,1,Nt+1);
    if PlotTimeV
recorderV = recorderV-VrIt;
    end
    if PlotSpaceV
recorderVx = recorderVx-VrIx;
    end
    if VtotPlot
recorderVtot = recorderVtot-VrItot;
    end
end
if DISPLAY
disp('4.1. Calculating membrane voltage plot');
end
Sx = ceil(S.*PlotVt);
St = unique(ceil((Nt+1).*PlotVx));
Pos = x(Sx);
Tim = t(St);

if PlotTimeV
LegendTimeV = cell(length(PlotVt),1);      % Vector to display the legend on the plots
figure;
if ~PlotTimeOne
r = ceil(sqrt(length(PlotVt)));
else 
hold on
end
for i = 1:length(PlotVt)
    if ~PlotTimeOne
subplot(r,r,i);
hold on;
    end
plot(t,recorderV(i,:))  
    if ~PlotTimeOne
h=legend(['Membrane-voltage at ' num2str(Pos(i)) 'm'],'Location','northwest');
set(h,'Interpreter','latex','FontSize',16);
tit = title(['location' num2str(i)]);
set(tit,'Interpreter','latex','Fontsize',17);
xlab = xlabel('t (s)'); set(xlab,'Interpreter','latex','Fontsize',15);
ylab = ylabel('V (mV)'); set(ylab,'Interpreter','latex','Fontsize',15);
hold off;
    else
LegendTimeV(i) = {['Membrane-voltage at ' num2str(Pos(i)) 'm']};
    end
 end
if PlotTimeOne
h=legend(LegendTimeV,'Location','northwest');
set(h,'Interpreter','latex','FontSize',16);   
tit = title('Membrane voltage at constant time');
set(tit,'Interpreter','latex','Fontsize',17);
xlab = xlabel('t (s)'); set(xlab,'Interpreter','latex','Fontsize',15);
ylab = ylabel('V (mV)'); set(ylab,'Interpreter','latex','Fontsize',15);
hold off  
end
end

if PlotSpaceV 
LegendSpaceV = cell(length(unique(ceil((Nt+1).*PlotVx))),1);
figure;
if ~PlotSpaceOne
r = ceil(sqrt(length(unique(ceil((Nt+1).*PlotVx)))));
else
hold on;
end
for i = 1:length(unique(ceil((Nt+1).*PlotVx)))
    if ~PlotSpaceOne
subplot(r,r,i);
hold on;
    end
plot(x,recorderVx(:,i))

    if ~PlotSpaceOne
h=legend(['Membrane-voltage at ' num2str(Tim(i)) 's'],'Location','northwest');
set(h,'Interpreter','latex','FontSize',16); 
tit = title(['Time' num2str(i)]);
set(tit,'Interpreter','latex','Fontsize',17);
xlab = xlabel('x (m)'); ylab = ylabel('V (mV)');
set(xlab,'Interpreter','latex','Fontsize',15);
set(ylab,'Interpreter','latex','Fontsize',15);
hold off;
    else
LegendSpaceV(i) = {['Membrane-voltage at ' num2str(Tim(i)) 's']};
    end
end  
if PlotSpaceOne
h=legend(LegendSpaceV,'Location','northwest');
set(h,'Interpreter','latex','FontSize',16);
tit = title('Membrane voltage at constant position');
set(tit,'Interpreter','latex','Fontsize',17);
xlab = xlabel('x (m)'); ylab = ylabel('V (mV)');
set(xlab,'Interpreter','latex','Fontsize',15);
set(ylab,'Interpreter','latex','Fontsize',15);
hold off;
end
end

% Plotting total V
if VtotPlot 
    figure;
    % Pseudocolormap of the membrane-voltage matrix
    pcolor(1000*t,1000*x,recorderVtot);
    shading interp
    % axis equal;       % Uncomment for equal data unit lengths
    % caxis([-2 2]);    % Uncomment to specify color variation 
    h = colorbar;       % Include colorbar
    tit = title(h,'Membrane voltage [mV]');
    set(tit,'Interpreter','latex','Fontsize',17);
    xlab = xlabel('t (ms)'); ylab = ylabel('x (mm)');
    set(xlab,'Interpreter','latex','Fontsize',15);
    set(ylab,'Interpreter','latex','Fontsize',15);
    
    if VexactPlot
    LegendVexact = cell(3,1);
    colorSpec = colorSpectrum(NFs);
    for iFreq = 1:NFs
    for iPath=1:NPTs
    pHandle(iFreq) = plot(tPathf{iPath,iFreq},xPathf{iPath},'LineWidth',2.5,'color',colorSpec(iFreq,:));
    plot(tPathb{iPath,iFreq},xPathb{iPath},'LineWidth',2.5,'color',colorSpec(iFreq,:));
    xlim([t(1) t(end)]);ylim([x(1) x(end)]);
    end
    LegendVexact(iFreq) = {['Activation pulse at ' num2str(fExact(iFreq)) 'Hz']};
    end
    legend(pHandle,LegendVexact,'location','northwest');
    end
end


% 9.2. Gate-parameters
if DISPLAY
disp('4.2. Calculating gate-parameter plots');
end
t = timeflow+(dt/2).*ones(1,Nt+1);        % Timeflow of gate parameters is shifted over dt/2
Tim = t(St);

% 9.2.1. m-parameter
if DISPLAY
disp('4.2.1. m-gate parameter');
end
if sum(Model~=0) ~= 0  
if PlotTimem 
LegendTimem = cell(length(Plotmt),1);
figure;
if ~PlotTimeOne
r = ceil(sqrt(length(Plotmt)));
else
hold on 
end
for i = 1:length(Plotmt)
    if ~PlotTimeOne
subplot(r,r,i);
hold on;
    end
plot(t,recorderm(i,:))
    if ~PlotTimeOne
h=legend(['Gate parameter m at location ' num2str(Pos(i)) 'm'],'Location','northwest');
set(h,'Interpreter','latex','FontSize',16);  
tit = title(['location' num2str(i)]);
set(tit,'Interpreter','latex','Fontsize',17);
xlab = xlabel('t (s)'); ylab = ylabel('m');
set(xlab,'Interpreter','latex','Fontsize',15);
set(ylab,'Interpreter','latex','Fontsize',15);
hold off;
    else
LegendTimem(i) = {['Gate parameter m at location ' num2str(Pos(i)) 'm']};
    end
end
if PlotTimeOne
h=legend(LegendTimem,'Location','northwest');
set(h,'Interpreter','latex','FontSize',16);  
tit = title('Gate parameter at constant location');
set(tit,'Interpreter','latex','Fontsize',17);
xlab = xlabel('t (s)'); set(ylab,'Interpreter','latex','Fontsize',15);
ylab = ylabel('m'); set(xlab,'Interpreter','latex','Fontsize',15);
hold off; 
end
end

if PlotSpacem 
figure;
LegendSpacem = cell(length(unique(ceil((Nt+1).*Plotmx))),1);
if ~PlotSpaceOne
r = ceil(sqrt(length(unique(ceil((Nt+1).*Plotmx)))));
else
hold on
end
for i = 1:length(Plotmx)
    if ~PlotSpaceOne
subplot(r,r,i);
hold on;
    end
plot(x,recordermx(:,i))
    if ~PlotSpaceOne
h=legend(['m-gate parameter at ' num2str(Tim(i)) 's'],'Location','northwest');
set(h,'Interpreter','latex','FontSize',16);   
tit = title(['Time' num2str(i)]);
set(tit,'Interpreter','latex','Fontsize',17);
xlab = xlabel('x (m)'); ylab = ylabel('m');
set(xlab,'Interpreter','latex','Fontsize',15);
set(ylab,'Interpreter','latex','Fontsize',15);
hold off;
    else
LegendSpacem(i) = {['m-gate parameter at ' num2str(Tim(i)) 's']};
    end
end 
if PlotSpaceOne
h=legend(LegendSpacem,'Location','northwest');
set(h,'Interpreter','latex','FontSize',16); 
tit = title('m-gate parameter at constant time');
set(tit,'Interpreter','latex','Fontsize',17);
xlab = xlabel('x (m)'); ylab = ylabel('m');
set(xlab,'Interpreter','latex','Fontsize',15);
set(ylab,'Interpreter','latex','Fontsize',15);
hold off;
end
end

if mtotPlot 
% Pseudocolormap of the m-gate parameter matrix
    figure;
    pcolor(t,x,recordermtot);
    shading interp
    % axis equal;               % Uncomment for equal data unit lengths 
    % caxis([-2 2]);            % Uncomment to specify color variation
    h = colorbar;                % Include colorbar
    tit = title(h,'m-gate parameter');
    set(tit,'Interpreter','latex','Fontsize',17);
    xlab = xlabel('t (s)'); ylab = ylabel('x (m)');
    set(xlab,'Interpreter','latex','Fontsize',15);
    set(ylab,'Interpreter','latex','Fontsize',15);
end

% 9.2.2. n-parameter
if DISPLAY
disp('4.2.2. n-gate parameter');
end
if sum(Model~=3)~=0
if PlotTimen 
figure;
LegendTimen = cell(length(Plotnt),1);
if ~PlotTimeOne
r = ceil(sqrt(length(Plotnt)));
else
hold on;
end
for i = 1:length(Plotnt)
    if ~PlotTimeOne
subplot(r,r,i);
hold on;
    end
plot(t,recordern(i,:))
    if ~PlotTimeOne
h=legend(['Membrane-voltage at ' num2str(Pos(i)) 'm'],'Location','northwest');
set(h,'Interpreter','latex','FontSize',16); 
tit = title(['location ' num2str(i)]); set(tit,'Interpreter','latex','Fontsize',17);
xlab = xlabel('t (s)'); set(xlab,'Interpreter','latex','Fontsize',15);
ylab = ylabel('n'); set(ylab,'Interpreter','latex','Fontsize',15);
hold off;
    else
LegendTimen(i) = {['Membrane-voltage at ' num2str(Pos(i)) 'm']};
    end
end
if PlotTimeOne
h=legend(['Membrane-voltage at ' num2str(Pos(i)) 'm'],'Location','northwest');
set(h,'Interpreter','latex','FontSize',16);  
tit = title('Membrane voltage at constant location');
set(tit,'Interpreter','latex','Fontsize',17);
xlab = xlabel('t (s)'); ylab = ylabel('n');
set(xlab,'Interpreter','latex','Fontsize',15);
set(ylab,'Interpreter','latex','Fontsize',15);
hold off;    
end
end

if PlotSpacen 
figure;
LegendSpacen = cell(length(unique(ceil((Nt+1).*Plotnx))),1);
if ~PlotSpaceOne
r = ceil(sqrt(length(unique(ceil((Nt+1).*Plotnx)))));
else
hold on;
end
for i = 1:length(unique(ceil((Nt+1).*Plotnx)))
    if ~PlotSpaceOne
subplot(r,r,i);
hold on;
    end
plot(x,recordernx(:,i))
    if ~PlotSpaceOne
h=legend(['n-gate parameter at time ' num2str(Tim(i)) 's'],'Location','northwest');
set(h,'Interpreter','latex','FontSize',16);   
tit = title(['Time ' num2str(i)]); set(tit,'Interpreter','latex','Fontsize',17);
xlab = xlabel('x (m)'); set(xlab,'Interpreter','latex','Fontsize',15);
ylab = ylabel('m'); set(ylab,'Interpreter','latex','Fontsize',15);
hold off;
    else
LegendSpacen(i) = {['n-gate parameter at time ' num2str(Tim(i)) 's']};
    end
end  
if PlotSpaceOne
h=legend(LegendSpacen,'Location','northwest');
set(h,'Interpreter','latex','FontSize',16);
tit = title('n-gate parameter at constant time');
set(tit,'Interpreter','latex','Fontsize',17);
xlab = xlabel('x (m)'); set(xlab,'Interpreter','latex','Fontsize',15);
ylab = ylabel('m'); set(ylab,'Interpreter','latex','Fontsize',15);
hold off;    
end
end

if ntotPlot 
% Pseudocolormap of the n-gate parameter matrix
    figure;
    pcolor(t,x,recorderntot);
    shading interp
    % axis equal;                 % Uncomment for equal data unit lengths 
    % caxis([-2 2]);              % Uncomment to specify color variation
    h = colorbar;                 % Include colorbar
    tit = title(h,'n-gate parameter');
    set(tit,'Interpreter','latex','Fontsize',17);
    xlab = xlabel('t (s)'); ylab = ylabel('x (m)');
    set(xlab,'Interpreter','latex','Fontsize',15);
    set(ylab,'Interpreter','latex','Fontsize',15);
end
end

% 9.2.3. h-parameter
if DISPLAY
disp('4.2.3. h-gate parameter');
end
if PlotTimeh 
figure;
LegendTimeh = cell(length(Plotht),1);
if ~PlotTimeOne
r = ceil(sqrt(length(Plotht)));
else
hold on;
end
for i = 1:length(Plotht)
    if ~PlotTimeOne
subplot(r,r,i);
hold on;
    end
plot(t,recorderh(i,:))
    if ~PlotTimeOne
h=legend(['Gate parameter h at ' num2str(Pos(i)) 'm'],'Location','northwest');
set(h,'Interpreter','latex','FontSize',16);   
tit = title(['location ' num2str(i)]); 
set(tit,'Interpreter','latex','Fontsize',17);
xlab = xlabel('t (s)'); set(xlab,'Interpreter','latex','Fontsize',15);
ylab = ylabel('h'); set(ylab,'Interpreter','latex','Fontsize',15);
hold off;
    else
LegendTimeh(i) = {['Gate parameter h at ' num2str(Pos(i)) 'm']};
    end
end
if PlotTimeOne
h=legend(LegendTimeh,'Location','northwest');
set(h,'Interpreter','latex','FontSize',16);      
tit = title('h-gate parameter at constant location');
set(tit,'Interpreter','latex','Fontsize',17);
xlab = xlabel('t (s)'); set(xlab,'Interpreter','latex','Fontsize',15);
ylab = ylabel('h'); set(ylab,'Interpreter','latex','Fontsize',15);
hold off;
end
end

if PlotSpaceh 
figure;
LegendSpaceh = cell(unique(ceil(length(Nt+1).*Plothx)),1);
if ~PlotSpaceOne
r = ceil(sqrt(unique(ceil(length((Nt+1).*Plothx)))));
else
hold on;
end
for i = 1:unique(ceil(length((Nt+1).*Plothx)))
    if ~PlotSpaceOne
subplot(r,r,i);
hold on;
    end
plot(x,recorderhx(:,i))
    if ~PlotSpaceOne
h=legend(['h-gate parameter at time ' num2str(Tim(i)) 's'],'Location','northwest');
set(h,'Interpreter','latex','FontSize',16); 
title(['Time ' num2str(i)])
xlab = xlabel('x (m)'); set(xlab,'Interpreter','latex','Fontsize',15);
ylab = ylabel('m'); set(xlab,'Interpreter','latex','Fontsize',15);
hold off;
    else
LegendSpaceh(i) = {['h-gate parameter at time ' num2str(Tim(i)) 's']};
    end
end  
if PlotSpaceOne
h=legend(['h-gate parameter at time ' num2str(Tim(i)) 's'],'Location','northwest');
set(h,'Interpreter','latex','FontSize',16);     
tit = title('h-gate parameter at constant time');
set(tit,'Interpreter','latex','Fontsize',17);
xlab = xlabel('x (m)'); set(xlab,'Interpreter','latex','Fontsize',15);
ylab = ylabel('m'); set(ylab,'Interpreter','latex','Fontsize',15);
hold off;  
end
end

if htotPlot 
% Pseudocolormap of the h-gate parameter matrix
    figure;
    pcolor(t,x,recorderhtot);
    shading interp
    % axis equal;                 % Uncomment for equal data unit lengths
    % caxis([-2 2]);              % Uncomment to specify color variation
    h = colorbar;                  % Include colorbar
    tit = title(h,'h-gate parameter');
    set(tit,'Interpreter','latex','FontSize',17);
    xlab = xlabel('t (s)'); ylab = ylabel('x (m)');
    set(xlab,'Interpreter','latex','Fontsize',15);
    set(ylab,'Interpreter','latex','Fontsize',15);
end

% 9.2.4. p-parameter
if DISPLAY
disp('4.2.4. p-gate parameter');
end
if sum((Model==2)|(Model==5))~=0
if PlotTimep
LegendTimep = cell(length(Plotpt),1);
figure;
if ~PlotTimeOne
r = ceil(sqrt(length(Plotpt)));
end
for i = 1:length(Plotpt)
    if ~PlotTimeOne
subplot(r,r,i);
hold on;
    end
plot(t,recorderp(i,:))
    if ~PlotTimeOne
h=legend(['Gate parameter p at location ' num2str(Pos(i)) 'm'],'Location','northwest');
set(h,'Interpreter','latex','FontSize',16); 
tit = title(['location ' num2str(i)]); set(tit,'Interpreter','latex','Fontsize',17);
xlab = xlabel('t (s)'); set(xlab,'Interpreter','latex','Fontsize',15);
ylab = ylabel('p'); set(ylab,'Interpreter','latex','Fontsize',15);
hold off;
    else
LegendTimep(i) = {['Gate parameter p at location ' num2str(Pos(i)) 'm']};
    end
end
if PlotTimeOne
h=legend(LegendTimep,'Location','northwest');
set(h,'Interpreter','latex','FontSize',16);   
tit = title('h-gate parameter at constant position');
set(tit,'Interpreter','latex','Fontsize',17);
xlab = label('t (s)'); set(xlab,'Interpreter','latex','Fontsize',17);
ylab = ylabel('p'); set(ylab,'Interpreter','latex','Fontsize',17);
hold off;
end
end

if PlotSpacep 
figure;
LegendSpacep = cell(length(unique(ceil((Nt+1).*Plotpx))),1);
if ~PlotSpaceOne
r = ceil(sqrt(length(unique(ceil((Nt+1).*Plotpx)))));
end
for i = 1:length(Plotpx)
    if ~PlotSpaceOne
subplot(r,r,i);
hold on;
    end
plot(x,recorderpx(:,i))
    if ~PlotSpaceOne
h=legend(['p-gate parameter at ' num2str(Tim(i)) 's'],'Location','northwest');
set(h,'Interpreter','latex','FontSize',16);  
tit = title(['Time ' num2str(i)]); set(tit,'Interpreter','latex','Fontsize',17);
xlab = xlabel('x (m)'); set(xlab,'Interpreter','latex','Fontsize',17);
ylab = ylabel('m'); set(ylab,'Interpreter','latex','Fontsize',17);
hold off;
    else
LegendSpacep(i) = {['p-gate parameter at ' num2str(Tim(i)) 's']};
    end
end  
if PlotSpaceOne
h=legend(LegendSpacep,'Location','northwest');
set(h,'Interpreter','latex','FontSize',16);      
h=legend('p-gate parameter','Location','northwest');
set(h,'Interpreter','latex','FontSize',16);
tit = title('p-gate parameter at constant time'); 
set(tit,'Interpreter','latex','Fontsize',17);
xlab = xlabel('x (m)'); set(xlab,'Interpreter','latex','Fontsize',15);
ylab = ylabel('m'); set(ylab,'Interpreter','latex','Fontsize',15);
hold off;
end
end

if ptotPlot 
    % Pseudocolormap of the p-gate parameter matrix
    figure;
    pcolor(t,x,recorderptot);
    shading interp
%     axis equal;                 % Uncomment for equal data unit lengths 
%     caxis([-2 2]);              % Uncomment to specify color variation
    h = colorbar;               % Include colorbar
    tit = title(h,'p-gate parameter'); set(tit,'Interpreter','latex','FontSize',15);
    xlab = xlabel('t (s)'); ylab = ylabel('x (m)');
    set(xlab,'Interpreter','latex','Fontsize',15); 
    set(ylab,'Interpreter','latex','Fontsize',15);
end
end
end

if Plotf
if DISPLAY
disp('4.3. Plotting activating function');
end
figure;
hold on;   
plot(x,f)   
h=legend('Activating function f','Location','northwest');
set(h,'Interpreter','latex','FontSize',16);  
tit = title('Activating function'); set(tit,'Interpreter','latex','FontSize',17);
xlab = xlabel('x (m)'); set(xlab,'Interpreter','latex','FontSize',15);
ylab = ylabel('f (mV/ms)'); set(ylab,'Interpreter','latex','FontSize',15);
hold off;
end

% 9.3. Theoretical MSOA
if DISPLAY
disp('4.3. Plotting velocities');
end
% 9.3.1. Phase velocity
if DISPLAY
disp('4.3.1. Phase velocity');
end
if PhaseVPlot
if PlotVelocityOne      % Only one plot is made of the phase velocity for all frequencies
while length(t) > 100       % Resampling because plotting several color maps in one plot is computationally heavy
    t = t(1:2:end);
    Vphase = Vphase(:,1:2:end,:);
end
[T,X,F] = meshgrid(t,x,fPhaseV);
figure;
hold on;
slice(T,X,F,Vphase,[],[],fPhaseV);
h = title('Phase velocity');
set(h,'Interpreter','latex','Fontsize',17);
xlab = xlabel('t (s)'); set(xlab,'Interpreter','latex','Fontsize',15);
ylab = ylabel('x (m)'); set(ylab,'Interpreter','latex','Fontsize',15);
zlab = zlabel('f (Hz)'); set(zlab,'Interpreter','latex','Fontsize',15);
hold off;

else
    
 % Pseudocolormap of all frequencies of the phase velocity
 for i = 1:length(fPhaseV)
    figure;
    pcolor(t,x,Vphase(:,:,i));
    shading interp
%     axis equal;                 % Uncomment for equal data unit lengths 
%     caxis([-2 2]);              % Uncomment to specify color variation
    h = colorbar;               % Include colorbar
    tit = title(h,'Phase velocity');
    set(tit,'Interpreter','latex','Fontsize',15);
    g = title(['Phase velocity at ' num2str(fPhaseV(i)) 'Hz']);
    set(g,'Interpreter','latex','Fontsize',17);
    xlab=xlabel('t (s)');ylab=ylabel('x (m)'); 
    set(xlab,'Interpreter','latex','Fontsize',15); set(ylab,'Interpreter','latex','Fontsize',15);
 end
end
end
% 9.3.2. Group velocity
if DISPLAY
disp('4.3.2. Group velocity');
end
if GroupVPlot
if PlotVelocityOne
while length(t) > 100       % Resampling because plotting several color maps in one plot is computationally heavy
    t = t(1:2:end);
    Vgroup = Vgroup(:,1:2:end,:);
end 
[T,X,F] = meshgrid(t,x,fGroupV);
figure;
hold on;
slice(T,X,F,Vgroup,[],[],fGroupV);
h = title('Group velocity');
set(h,'Interpreter','latex','FontSize',17);
xlab = xlabel('t (s)'); set(xlab,'Interpreter','latex','Fontsize',15);
ylab = ylabel('x (m)'); set(ylab,'Interpreter','latex','Fontsize',15);
zlab = zlabel('f (Hz)'); set(zlab,'Interpreter','latex','Fontsize',15);
hold off;
else
    
 % Pseudocolormap of all frequencies of the group velocity
 for i = 1:length(fGroupV)
    figure;
    pcolor(t,x,Vgroup(:,:,i));
    shading interp
%     axis equal;                 % Uncomment for equal data unit lengths 
%     caxis([-2 2]);              % Uncomment to specify color variation
    h = colorbar;                 % Include colorbar
    tit = title(h,'Group velocity');
    set(tit,'Interpreter','latex','Fontsize',15);
    g = title(['Group velocity at ' num2str(fGroupV(i)) 'Hz']);
    set(g,'Interpreter','latex','Fontsize',17);
    xlab=xlabel('t (s)');ylab=ylabel('x (m)'); 
    set(xlab,'Interpreter','latex','Fontsize',15); set(ylab,'Interpreter','latex','Fontsize',15);
 end
end
end
% 9.3.3. Integrate group velocity (Vexact)
if DISPLAY
disp('4.3.2. Integrate group velocity');
end
if VexactPlot
 
 % Pseudocolormap of all frequencies of the group velocity
 for iFreq = 1:length(fExact)
    figure;
    hold on;
    pcolor(t,x,Vexact(:,:,iFreq));
    shading interp
%     axis equal;                 % Uncomment for equal data unit lengths 
%     caxis([-2 2]);              % Uncomment to specify color variation
    h = colorbar;                 % Include colorbar
    tit = title(h,'Group velocity');
    set(tit,'Interpreter','latex','Fontsize',15);
    g = title(['Group velocity at ' num2str(fExact(iFreq)) 'Hz']);
    set(g,'Interpreter','latex','Fontsize',17);
    xlab=xlabel('t (s)');ylab=ylabel('x (m)'); 
    set(xlab,'Interpreter','latex','Fontsize',15); set(ylab,'Interpreter','latex','Fontsize',15);
    for iPath=1:NPTs
    plot(tPathf{iPath,iFreq},xPathf{iPath},'LineWidth',2.5,'Color','k');
    xlim([t(1) t(end)]);ylim([x(1) x(end)]);
    plot(tPathb{iPath,iFreq},xPathb{iPath},'LineWidth',2.5,'Color','k');
    xlim([t(1) t(end)]);ylim([x(1) x(end)]);
    end
    hold off;
 end
end


else            % SweepActI = 1 -> post-processing Sweep parameters
    % 1) Rheobase = f(x,y,z,ancat)
Rheobase = permute(SweepActImat(:,:,:,end,:,1),[1 2 3 5 4]);  % Rheobase obtained by longest monophasic pulse duration
    % 2) Polarisation ratio = f(x,y,z,Tp,mb)
PolRatio = permute(SweepActImat(:,:,:,:,1,:)./SweepActImat(:,:,:,:,2,:),[1 2 3 4 6 5]);
    % 3) Minimum charge Qe = f(x,y,z,ancat)
Qe = permute(Tps(1)*SweepActImat(:,:,:,1,:,1),[1 2 3 5 4]);
    % 4) SD time constant SDtau = f(x,y,z,ancat)
SDtau = Qe./Rheobase;
    % 5) Bi/monophasic ratio Rbm = f(x,y,z,Tp,ancat)
Rbm = SweepActImat(:,:,:,:,:,2)./SweepActImat(:,:,:,:,:,1);
    % 6) Store results
	 if LLancat == 1
        if LLmonbi == 1
    SweepAct = SweepActImat(1,1,1,1,1,1);
        else if LLmonbi == 2
                SweepAct = SweepActImat(1,1,1,1,1,2);
            end
        end
    else if LLancat == 2
            if LLmonbi == 1
                SweepAct = SweepActImat(1,1,1,1,2,1);
            else if LLmonbi == 2
                    SweepAct = SweepActImat(1,1,1,1,2,2);
                end
            end
        end
    end
    switch SENNmodel
        case 1 
            modelstr = 'HH';
        case 2
            modelstr = 'FH';
        case 3
            modelstr = 'CRRSS';
        case 4
            modelstr = 'SE';
        case 5
            modelstr = 'SRB';
    end
save(['SweepSENNSR-' num2str(LLx0s) '-' num2str(LLy0s) '-' num2str(LLz0s) '-'...  
    num2str(LLTp,5) '-' num2str(LLancat) '-' num2str(LLmonbi) '(' modelstr '-Conf' num2str(LLconf) ').mat'],'SweepAct');
%save('x0sValuesSENN_FH_Conf1.mat','x0sValues');
%save('y0sValuesSENN_FH_Conf1.mat','y0sValues');
%save('z0sValuesSENN_FH_Conf1.mat','z0sValues');
%save('TpsValuesSENN_FH_Conf1.mat','TpsValues');

TotalTime = toc/60;                            % End stopwatch timer
%save('TT_SENN_FH_Conf1.mat','TotalTime');
if ~SweepActI
if DISPLAY
disp('5. Rendering plots and saving data');
end
else 
if DISPLAY
disp('2. Rendering plots and saving data');
end
end
end
end