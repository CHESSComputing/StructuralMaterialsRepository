function mesh = MeshStructureGrain(crd, con, npcrossref,npicrossref,surfnp,solid4surf,surfcentroid,surfnormal,surfarea,eqv, varargin)
% MeshStructureGrain - Create mesh structure from mesh data related to a grain.
% Used for grain-by-grain calculations in GrainbyGrainMesh_Modes.
%   
%
%-------------------- Defaults and Keyword Arguments
%
optcell = {...
    'ElementType',  '', ...
    'Symmetries',   []  ...
       };
%
options = OptArgs(optcell, varargin);
%
%-------------------- Main Code
%
if (nargin == 0)
  crd = [];
  con = [];
  eqv = [];
elseif (nargin == 2)
  eqv = [];
end
%
mesh = struct('crd', crd, 'con', con,'npcrossref', npcrossref,'npicrossref', npicrossref,...
    'surfnp',surfnp,'solid4surf',solid4surf,'surfcentroid',surfcentroid,'surfnormal',surfnormal,'surfarea',surfarea,...
    'eqv', eqv);
%
%  Process optional arguments.
%
if nargin > 3
  %
  if options.ElementType
    mesh.etype = ElementTypeStruct(options.ElementType);
  end
  %
  if ~isempty(options.Symmetries)
    mesh.symmetries = options.Symmetries;
  end
  %
end
