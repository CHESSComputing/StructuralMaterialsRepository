function  shafac_structure    = Data_Structure_Shafac(nnpe,nqpt,wtqp,sfac,dndxi,dndet,dndze)
%
%  Data_Structure_Shafac
%
%  Structure file used in ShapeFunctionArrays
%

%
shafac_structure = struct('nnpe',nnpe,'nqpt',nqpt,'wtqp',wtqp,'sfac',sfac,'dndxi',dndxi,'dndet',dndet,'dndze',dndze);