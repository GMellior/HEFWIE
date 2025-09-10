% Solves transitions from NICL to
% Tuition subsidies funded by everyone (TS)
% Tuition subsidies funded with graduate taxes (GT)

clear
clc

addpath('Functions')                                                   
subsvec          = [0.5,0.75,1,1]';                                        % Subsidy %, first 3 TS, then 1 GT
gi               = [0;0;0;1];                                              % 1 -> GT - 0 -> TS
Res              = cell(numel(subsvec,1),1);                               % Store some results (will generate .mat files anyways)
gridpoints       = [205,229,229,229];                                      % Grid points in b
   
for hhh=1:numel(subsvec)
    Res{hhh} = NICLtoTSGT(subsvec(hhh,1),gi(hhh,1),gridpoints(hhh));
end