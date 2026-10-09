function  bc_structure   = Data_Structure_BCs(bc_code,bc_value, ess_or_nat_code, force_or_disp, moment_or_angle, inertiavalues, neutralaxisvalues)
%
% Data_Structure_BCs
%
%   Creates a MatLab structure file with data for application of the
%   boundary conditions (natural and essential)
%   Used in MechMet and MechModuli
%


%
bc_structure = struct('bc_code',bc_code,'bc_value',bc_value, 'ess_or_nat_code', ess_or_nat_code, 'force_or_disp', force_or_disp, ...
    'moment_or_angle', moment_or_angle, 'inertiavalues', inertiavalues, 'neutralaxisvalues', neutralaxisvalues);