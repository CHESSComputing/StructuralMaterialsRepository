function microstructurefile = MicrostructureStructure(phases,grains,orientations,rotations)
%
% MicrostructureStructure - Create microstructure structure from grain phase and orientation data.
%                         - used in MechMet and Mechmonics


microstructurefile = struct('phases',phases,'grains',grains,'orientations',orientations,'rotations',rotations);
