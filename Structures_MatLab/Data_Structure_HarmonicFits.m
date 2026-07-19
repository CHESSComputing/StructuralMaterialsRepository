function  harmonicfitfile    = Data_Structure_HarmonicFits(field_data,field_modes_orth,beta_all_orth,error_all_orth)

%
% Data_Structure_HarmonicFits
%
%   Creates a MatLab structure file with data from fitting discrete harmonics to FE data. 
%  Used in MechMet

%
harmonicfitfile = struct('field_data',field_data,'field_modes_orth',field_modes_orth,'beta_all_orth',beta_all_orth,'error_all_orth',error_all_orth);
