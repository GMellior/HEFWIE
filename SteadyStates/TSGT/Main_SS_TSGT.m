% Solves steady states with no student loans. 
% Tuition subsidies funded by everyone (TS)
% Tuition subsidies funded with graduate taxes (GT)

clear
clc

addpath('Functions')
Edpricescale     = 1;                                                      % Scale education costs
subsvec          = [0.5,0.75,1,0.5,0.75,1]';                               % Subsidy %, first 3 GT, then 3 TS
gi               = [1;1;1;0;0;0];                                          % 1 -> GT - 0 -> TS
input_vector     = [subsvec,gi];
Res              = cell(size(input_vector,1),1);                           % Store some results (will generate .mat files anyways)
smallgrid        = [1,1,0,1,0,0];                                          % Set to 1 for fastest results
   
for hhh=1:size(input_vector,1)
    Res{hhh} = SSTSGT(subsvec(hhh,1),gi(hhh,1),Edpricescale,smallgrid(hhh));
end

