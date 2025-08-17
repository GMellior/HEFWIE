% Solves steady states with student loans. 
% Baseline (NICL), ICL1 and ICL2

clear
clc

addpath('Functions')
Edpricescale     = 1;                                                      % Scale education costs
subsvec          = 0.47*ones(3,1);                                         % Subsidy % - fixed at baseline value
SLindicator      = [0:2]';                                                 % 0 -> ICL2 - 1 -> NICL - 2 -> ICL1
aminfedpctg      = 1;
input_vector     = [subsvec,SLindicator];
Res              = cell(numel(SLindicator),1);                             % Store some results (will generate .mat files anyways)

for hhh=1:numel(SLindicator)
    Res{hhh} = SSNICLandICLs(subsvec(hhh,1),0,Edpricescale,aminfedpctg,SLindicator(hhh));
end