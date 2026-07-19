function propertyfile = CrystalPropertyStructure(crystal_type_all,r_matrix_all,sss_all,sssr_all,covera_all)
%
% CrystalPropertyStructure
%
%   Creates a MatLab structure file with crystal properties
%

propertyfile = struct('crystal_type_all',crystal_type_all,'r_matrix_all',r_matrix_all,'sss_all',sss_all,'sssr_all',sssr_all,'covera_all',covera_all);

