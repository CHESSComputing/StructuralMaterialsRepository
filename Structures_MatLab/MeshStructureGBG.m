function mesh = MeshStructureGBG(crd, con, nepercrd, nepercon, old2new, grn4np, nepergrn4np, eqv, varargin)
% MeshStructure - Create mesh structure from mesh data.
%   
%   USAGE:
%
%   mesh = MeshStructure(crd, con, nepercrd,nepercon, nepergrn4np, eqv)
%   mesh = MeshStructure(crd, con, eqv, <options>)
%
%   INPUT:
%
%   crd is e x n
%       the array of nodal point locations
%   con is d x m  (integer)
%       the mesh connectivity
%   grn4np is a np x 6 (integer)
%       grains for each nodal point
%   eqv is 2 x k, (integer, optional) 
%       the equivalence array
%
%   <options> is a sequence of parameter, value pairs
%             Available options are listed below with default values
%             shown in brackets.
%
%   'ElementType'    string
%                    See ElementTypeStruct for available types.
%   'Symmetries'     4 x n
%                    an array of quaternions for crystal symmetries
%
%   OUTPUT:
%
%   mesh is a MeshStructure, 
%        the basic MeshStructure consists of three fields,
%        the nodal point coordinates (.crd), the connectivity
%        (.con) and the nodal point equivalence array (.eqv).
%        
%   NOTES:
% 
%   *  With no arguments, this function returns and empty mesh
%      structure.  With only two arguments, it sets the equivalence 
%      array to be empty.
%
%
MyName = mfilename;
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
  grn4np = [];
  eqv = [];
elseif (nargin == 2)
  eqv = [];
end
%
mesh = struct('crd', crd, 'con', con, 'nepercrd', nepercrd, ...
    'nepercon', nepercon, 'old2new', old2new, 'grn4np', grn4np, 'nepergrn4np', nepergrn4np,...
 'eqv', eqv);
%
%  Process optional arguments.
%
if nargin > 4
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
