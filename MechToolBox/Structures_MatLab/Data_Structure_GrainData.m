function  graindatafile    = Data_Structure_GrainData(grain_volumes,rotations,grain_schmidfactors)

%
% Data_Structure_GrainData
%
%   Creates a MatLab structure file with data from MechMet simulations 
%  Used in MechMet

%
graindatafile = struct('grain_volumes',grain_volumes,'rotations',rotations,'schmidfactors',grain_schmidfactors);