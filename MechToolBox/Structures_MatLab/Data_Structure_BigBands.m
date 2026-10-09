function  bigbandgraphfile = Data_Structure_BigBands(graph_nodenames,graph_elements,graph_links_unique,graph_link_ele_unique,band_global_link)
%
% Data_Structure_BigBands
%
%   Creates a MatLab structure file with data related to graphs associated
%   with connected sets (two-way) 
%   Used in MechBands

%
bigbandgraphfile = struct('graph_nodenames',graph_nodenames,'graph_elements',graph_elements,...
    'graph_links_unique',graph_links_unique,'graph_link_ele_unique',graph_link_ele_unique,...
     'band_global_link',band_global_link);