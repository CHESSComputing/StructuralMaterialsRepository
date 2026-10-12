function bc_structure =   BC_Definition(bc_filename,meshfile,inertiavalues,neutralaxisvalues,areavalues)
%
% read in boundary condition data from file with basename: bc_filename, and
%     extension: bcs (eg bc_filename.bcs).
% define and return (in bc_structure) the bccode and bcvalue arrays.  
%     Note that the surface force arrays needed in MechModuli are computed in separate function files.
%
% This function is called by MechMet and by MechModuli.
%   MechMet and MechModuli differ in the type of boundary conditions applied.  
%   MechMet applies only essential BCs to impose displacements for axial extension or bending rotations. 
%   MechModuli applies only natural BCs to impose tractions for tensile forces or bending moments.
%    Both MechMech and MechModuli apply essential boundary conditions to elminate rigid body modes.
%   
%  Boundary conditions are defined surface by surface.
%     Surfaces are the 6 surfaces given in Neper output under fasets
%     Edges and vertices are identified in this function.
%     BCs for any surfaces not defined here are taken to be natural BCs
%       with zero tractions (general fem default condition for displacement-based elasticity codes)
%     In accordance with the comment above:
%       when using MechMet, prescribe only essential bcs on the surfaces, or  
%       when using MechModuli, prescribe only natural bcs on the surfaces.
%       For both codes, either rigid body option (see below) can be used for removing rigid body modes 
%     
%  Information in bc_filename.bcs:

%   To eliminate rigid body modes, essential body conditions are set on up
%   to 4 nodes.  
%    For pure extension, the point of min x,y and z if fixed and rigid body
%    rotations are eliminated by setting u_z=0 at xmax, u_x=0 at ymax, and
%    u_y=0 at zmax.  
%    For problems with bending imposed, not all the rotation constraints
%    are needed.  If there is bending about the x axis, then the essential
%    bc at zmax is dropped. If there is bending about the y axis, then the essential
%    bc at xmax is dropped. If there is bending about the z axis, then the essential
%    bc at ymax is dropped. 
%
%   For each surface have bc data: the following input:    (One set for each surface with applied bc ) 
%     First line:   surfindentifier: integers between 1 and 6 according as designated using fasets in neper:  x0,x1,y0,y1,z0,z1) 
%     Second line:  essential/natural code:  0 for natural or 1 for essential  for each degree of freedom:  bccode  integer
%     Third line:   normal traction or displacement perpendicular to the surface -reals
%     Fourth line:  moment or rotation about each of the three coordinate direction  -reals
%       However, of the three possible directions, only the two in-plane diretions should by non-zero; the out-of-plane direction should be zero


edges4surf = [3 4 5 6;3 4 5 6; 1 2 5 6; 1 2 5 6; 1 2 3 4; 1 2 3 4];
surfs4edge = [ 3 5; 2 5; 4 5; 1 5; 2 3; 2 4; 1 4; 1 3];


% form list of surface nodes from nps arrays for all elements on a surface
nplist_s1 = unique(meshfile.nps_1(:,[2 3 4 5 6 7]));
nplist_s2 = unique(meshfile.nps_2(:,[2 3 4 5 6 7]));
nplist_s3 = unique(meshfile.nps_3(:,[2 3 4 5 6 7]));
nplist_s4 = unique(meshfile.nps_4(:,[2 3 4 5 6 7]));
nplist_s5 = unique(meshfile.nps_5(:,[2 3 4 5 6 7]));
nplist_s6 = unique(meshfile.nps_6(:,[2 3 4 5 6 7]));

nplist_s_size = [ size(nplist_s1,1) size(nplist_s2,1) size(nplist_s3,1) size(nplist_s4,1) size(nplist_s5,1) size(nplist_s6,1)];

% form list of nodes on edges from intersections of surface np arrays
nplist_e1 = intersect(nplist_s3,nplist_s5);
nplist_e2 = intersect(nplist_s2,nplist_s5);
nplist_e3 = intersect(nplist_s4,nplist_s5);
nplist_e4 = intersect(nplist_s1,nplist_s5);
nplist_e5 = intersect(nplist_s2,nplist_s3);
nplist_e6 = intersect(nplist_s2,nplist_s4);
nplist_e7 = intersect(nplist_s1,nplist_s4);
nplist_e8 = intersect(nplist_s1,nplist_s3);
nplist_e9 = intersect(nplist_s3,nplist_s6);
nplist_e10 = intersect(nplist_s2,nplist_s6);
nplist_e11 = intersect(nplist_s4,nplist_s6);
nplist_e12 = intersect(nplist_s1,nplist_s6);

nplist_e_size = [ size(nplist_e1,1) size(nplist_e2,1) size(nplist_e3,1) size(nplist_e4,1) size(nplist_e5,1) size(nplist_e6,1) ...;
    size(nplist_e7,1) size(nplist_e8,1) size(nplist_e9,1) size(nplist_e10,1) size(nplist_e11,1) size(nplist_e12,1)];

% find vertex nodes from intersection of edges
nplist_v1 = intersect(nplist_e4,nplist_e1);
nplist_v2 = intersect(nplist_e1,nplist_e2);
nplist_v3 = intersect(nplist_e3,nplist_e4);
nplist_v4 = intersect(nplist_e2,nplist_e3);
nplist_v5 = intersect(nplist_e12,nplist_e9);
nplist_v6 = intersect(nplist_e9,nplist_e10);
nplist_v7 = intersect(nplist_e11,nplist_e12);
nplist_v8 = intersect(nplist_e10,nplist_e11);

%  read the .bcs file for all bc data

bcfile = [bc_filename,'.bcs'];
fid = fopen(bcfile,'r');

surf_bc_specified = zeros(6,1);
ess_or_nat_code = zeros(6,1);
force_or_disp = zeros(6,1);
moment_or_angle = zeros(6,3);
    
while ~feof(fid)

% read in set of bc data for a surface
    sid = fscanf(fid, '%d', 1);
    eorn = fscanf(fid, '%d', 1); 
    fd = fscanf(fid, '%f', 1); 
    mv = fscanf(fid, '%f %f %f', 3);    

% put these values into arrays for all six surfaces
    surf_bc_specified(sid) = 1;
    ess_or_nat_code(sid) = eorn;
    force_or_disp(sid,:) = fd(:);
    moment_or_angle(sid,:) = mv(:);

end

% apply boundary conditions specified in BC file to each node the designated surfaces

    bc_code = [];
    bc_value = [];
    numnps_total = 0;

    codes=[1 0 0; 1 0 0; 0 1 0; 0 1 0; 0 0 1; 0 0 1];

    rot_xaxis = 0;
    rot_yaxis = 0;
    rot_zaxis = 0;

for isurf = 1:1:6

    numnps = nplist_s_size(isurf);                            
    rotang1 = moment_or_angle(isurf,1);
    rotang2 = moment_or_angle(isurf,2);
    moi1 = inertiavalues(isurf,1);
    moi2 = inertiavalues(isurf,1);
    nuax1 = neutralaxisvalues(isurf,1);
    nuax2 = neutralaxisvalues(isurf,2);

    rotangles = moment_or_angle(isurf,:);
    mois = inertiavalues(isurf,:);
    nuaxes = neutralaxisvalues(isurf,:);

    switch isurf
       case 1
            nplist = nplist_s1;
       case 2 
            nplist = nplist_s2;  
       case 3
            nplist = nplist_s3;
       case 4
            nplist = nplist_s4; 
       case 5
            nplist = nplist_s5;
       case 6
            nplist = nplist_s6;
    end

    if(ess_or_nat_code(isurf)==1)
    numnps_total = numnps_total+numnps;
        for j = 1:1:numnps
                npj = nplist(j);
                xj = meshfile.crd(1,npj);
                yj = meshfile.crd(2,npj);
                zj = meshfile.crd(3,npj);
                bccods = [npj; codes(isurf,1); codes(isurf,2); codes(isurf,3)];
                bc_code = [bc_code bccods];
                    switch isurf
                        case 1
                            npvalb = [(zj-nuaxes(3))*rotangles(2)+(yj-nuaxes(2))*rotangles(3) 0 0];
                            npvale = [force_or_disp(1) 0 0];
                            if(rotangles(2) ~= 0.0)
                                rot_yaxis = 1;
                            end
                            if(rotangles(3) ~= 0.0)
                                rot_zaxis = 1;
                            end
                        case 2 
                            npvalb = [(zj-nuaxes(3))*rotangles(2)+(yj-nuaxes(2))*rotangles(3) 0 0];
                            npvale = [force_or_disp(2) 0 0];
                            if(rotangles(2) ~= 0.0)
                                rot_yaxis = 1;
                            end
                            if(rotangles(3) ~= 0.0)
                                rot_zaxis = 1;
                            end
                        case 3
                            npvalb = [0 (xj-nuaxes(1))*rotangles(3)+(zj-nuaxes(3))*rotangles(1) 0];
                            npvale = [0 force_or_disp(3) 0];
                            if(rotangles(3) ~= 0.0)
                                rot_zaxis = 1;
                            end
                            if(rotangles(1) ~= 0.0)
                                rot_xaxis = 1;
                            end
                        case 4
                            npvalb = [0 (xj-nuaxes(1))*rotangles(3)+(zj-nuaxes(3))*rotangles(1) 0];
                            npvale = [0 force_or_disp(4) 0];
                            if(rotangles(3) ~= 0.0)
                                rot_zaxis = 1;
                            end
                            if(rotangles(1) ~= 0.0)
                                rot_xaxis = 1;
                            end
                        case 5
                            npvalb = [0 0 (yj-nuaxes(2))*rotangles(1)+(xj-nuaxes(1))*rotangles(2)];
                            npvale = [0 0 force_or_disp(5)];
                            if(rotangles(1) ~= 0.0)
                                rot_xaxis = 1;
                            end
                            if(rotangles(2) ~= 0.0)
                                rot_zaxis = 1;
                            end
                        case 6
                            npvalb = [0 0 (yj-nuaxes(2))*rotangles(1)+(xj-nuaxes(1))*rotangles(2)];
                            npvale = [0 0 force_or_disp(6)];
                            if(rotangles(1) ~= 0.0)
                                rot_xaxis = 1;
                            end
                            if(rotangles(2) ~= 0.0)
                                rot_zaxis = 1;
                            end
                    end

                bcvals = [npj; npvale(1)+npvalb(1); npvale(2)+npvalb(2); npvale(3)+npvalb(3)];
                bc_value = [bc_value bcvals];
        end
    end

end

bc_code= bc_code';
bc_value = bc_value';


% add bc assignments as needed to remove rigid body modes


    bcods = [nplist_v1(1) 1 1 1];
    bcvals = [nplist_v1(1) 0 0 0];
    numnps_total = numnps_total+1;
    bc_code = cat(1,bc_code,bcods);
    bc_value = cat(1,bc_value,bcvals);


    if(rot_yaxis==0)
    bcods = [nplist_v2(1) 0 0 1];
    bcvals = [nplist_v2(1) 0 0 0];
    numnps_total = numnps_total+1;
    bc_code = cat(1,bc_code,bcods);
    bc_value = cat(1,bc_value,bcvals);
    end

    if(rot_zaxis==0)
    bcods = [nplist_v3(1) 1 0 0];
    bcvals = [nplist_v3(1) 0 0 0];
    numnps_total = numnps_total+1;
    bc_code = cat(1,bc_code,bcods);
    bc_value = cat(1,bc_value,bcvals);
    end

    if(rot_xaxis==0)
    bcods = [nplist_v5(1) 0 1 0];
    bcvals = [nplist_v5(1) 0 0 0];
    numnps_total = numnps_total+1;
    bc_code = cat(1,bc_code,bcods);
    bc_value = cat(1,bc_value,bcvals);
    end

[uniquenodes,ia,ic] = unique(bc_code(:,1));

numuniquebc = size(uniquenodes,1);
bc_code_unique = zeros(numuniquebc,4);
bc_value_unique = zeros(numuniquebc,4);

for iunique = 1:1:numuniquebc
    bc_code_unique(iunique,1) = uniquenodes(iunique);
    bc_value_unique(iunique,1) = uniquenodes(iunique);
end

for itotal = 1:1:numnps_total
    bc_code_unique(ic(itotal),[2 3 4]) = bc_code_unique(ic(itotal),[2 3 4])+bc_code(itotal,[2 3 4]);
    bc_value_unique(ic(itotal),[2 3 4]) = bc_value_unique(ic(itotal),[2 3 4])+bc_value(itotal,[2 3 4]);
end

for iunique = 1:1:numuniquebc
    if(bc_code_unique(iunique,2)>1)
        bc_code_unique(iunique,2)=1;
    end
    if(bc_code_unique(iunique,3)>1)
        bc_code_unique(iunique,3)=1;
    end
    if(bc_code_unique(iunique,4)>1)
        bc_code_unique(iunique,4)=1;
    end
end


%  write bccode and bcvalue arrays to bcstructure

  bc_structure = Data_Structure_BCs(bc_code_unique,bc_value_unique, ess_or_nat_code, force_or_disp, moment_or_angle, inertiavalues, neutralaxisvalues);
%
status = fclose(fid);
