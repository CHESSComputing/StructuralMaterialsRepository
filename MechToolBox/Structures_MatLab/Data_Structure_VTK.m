function  bandsdatafile    = Data_Structure_VTK(element_bands,element_numbers,datafile)
%
%   Data_Structure_VTK
%
% Structure file used in MechBands to write VTK file
%

%
bandsdatafile = struct('element_bands',element_bands,'element_numbers',element_numbers, 'datafile', datafile);