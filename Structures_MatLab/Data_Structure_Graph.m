function  graphdatafile    = Data_Structure_Graph(element_bands,element_numbers,datafile)

%
% Data_Structure_Graph
%
%   Creates a MatLab structure file with data in MechBands simulations 
%  

%
graphdatafile = struct('element_bands',element_bands,'element_numbers',element_numbers, 'datafile', datafile);