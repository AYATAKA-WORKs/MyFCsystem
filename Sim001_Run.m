clear
%% Initalize
Sim000_Initialize

%% Simulate
tic % Start timer
out = sim(".\SLX\Integrated_Model.slx");
toc % Get simulation time