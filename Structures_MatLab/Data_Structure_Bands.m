function  bandarrays    = Data_Structure_Bands(nbands,bandsizes,bands,bigbands,bandorder,element_bands)
%
% Data_Structure_Bands
%
%   Creates a MatLab structure file with band data
%   Used in MechBands
%

%
bandarrays = struct('nbands',nbands,'bandsizes',bandsizes,'bands',bands,...
    'bigbands',bigbands,'bandorder',bandorder,'element_bands',element_bands);