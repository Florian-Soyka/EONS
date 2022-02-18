% Parameter 2 CSV-file.
% ------------------------------------------------------------------
% Input parameters are: 
% model (1-5, HH - SRB), D (mm), L (mm), max_step
%(us), order of solution, type (1=sine, 2=pulse), frequency
%(Hz), beat (Hz), duration (ms), lower limit (V/m), upper limit (V/m), 
% nAPs (-), APdur (s)
clear all; %#ok<CLALL>

% Define parameter names for .csv file format
ParNames = {'model','D','L','max_step','orderofsolution','type','frequency','beat','duration',...
    'lowerlimit','upperlimit','nAPs','APdur'};
saveName = 'BioEM2022_Abstract2_fullbatch';

% All combinations of the input parameter arrays are generated.
% The Val-parameters indicate how the parameters are defined:
% Val = 1 -> each element in the array indicates a value of the parameter
% Val = 2 -> Par(2) linearly distributed elements between Par(1) and Par(3) 
% Val = 3 -> Par(2) logarithmicly distributed elements between Par(1) and Par(3)
modelVal = 2;
DVal = 1;
LVal = 1;
max_stepVal = 1;
orderofsolutionVal = 1;
typeVal = 1;
frequencyVal = 1;
beatVal = 1;
durationVal = 1;
lowerlimitVal = 1;
upperlimitVal = 1;
nAPsVal = 2;
APdurVal = 1;
% ---COMBINE ALL VALUES IN PARVal-----------------------------------------
PARVal = [modelVal,DVal,LVal,max_stepVal,orderofsolutionVal,typeVal,frequencyVal,beatVal,durationVal,...
    lowerlimitVal,upperlimitVal,nAPsVal,APdurVal];
% ------------------------------------------------------------------------
model = [1,5,5]; %#ok<*NBRAK>
D = [0.02];
L = [20.5];
max_step = [25];
orderofsolution = [2];
type = [1];
frequency = [1,2,3,5,10,20,30,50,100,200,300,500,1000,2000,3000,...
    5000,10000,20000,30000,50000,100000];
beat = [0];
duration = 1000*min([1./frequency,0.025]);
lowerlimit = [0];
upperlimit = [0];
nAPs = [1,5,5];
APdur = 0.001;
% -------------COMBINE ALL PARAMETERS IN PAR cell--------------------------
PAR = {model;D;L;max_step;orderofsolution;type;frequency;beat;duration;lowerlimit;...
    upperlimit;nAPs;APdur};
% -------------------------------------------------------------------------
% CALCULATE EXPLICIT PARValues
PARValues = cell(length(PAR),1); NPARStep = zeros(length(PAR),1);
for i=1:length(PAR)
if PARVal(i) == 1
PARValues{i} = PAR{i};
else
if PAR{i}(1) == PAR{i}(3)
PARValues{i} = (PAR{i}(1));  
else
if PARVal(i) == 3
PARValues{i} =(PAR{i}(3)/PAR{i}(1)).^((0:1:(PAR{i}(2)-1))/(PAR{i}(2)-1))*PAR{i}(1);
elseif PARVal(i) == 2 
PARStep = (PAR{i}(3)-PAR{i}(1))/(PAR{i}(2)-1);
PARValues{i} = (PAR{i}(1):PARStep:PAR{i}(3));
end
end
end
NPARStep(i) = length(PARValues{i});
end

Config = zeros(prod(NPARStep),length(PAR));
digitVal = cumprod(NPARStep(end:-1:1));
digitVal = [digitVal(end-1:-1:1);1];
IndexList = cell(prod(NPARStep),1);
for i=0:prod(NPARStep)-1
    Digit=i; Index = ones(1,length(PAR));
for j=1:length(PAR)
while Digit>=digitVal(j)
Index(j) = Index(j)+1;
Digit = Digit-digitVal(j);
end
end
IndexList{i+1} = Index;
end

for i=1:size(Config,1)
Config(i,:) = cellfun(@(x,y) x(y),PARValues',num2cell(IndexList{i}));
end
Config = vertcat(ParNames,num2cell(Config)); 
xlswrite(saveName,Config);