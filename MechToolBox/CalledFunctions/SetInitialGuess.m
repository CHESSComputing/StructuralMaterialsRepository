function [u] = SetInitialGuess(numnp, x,y,z,lenxyz,minxyz,loadcode)
%
% SetInitialGuess  -- defines a guess for the displacements base on the
% loading type
%
%  called by MechMet
%
%  
%  the code is defined as follows:
%   loadcode = 1  --  x-tension
%   loadcode = 2  --  y-tension
%   loadcode = 3  --  z-tension
%   loadcode = 4  --  x-y biaxial tension
%   loadcode = 5  --  y-z biaxial tension
%   loadcode = 6  --  z-x biaxial tension    
%
%   Does a linear interpolation across the mesh
%
%   OUTPUT:
%
%     displacement vector, u
%

u = zeros(3*numnp,1);
lenx = lenxyz(1);
leny = lenxyz(2);
lenz = lenxyz(3);
minx = minxyz(1);
miny = minxyz(2);
minz = minxyz(3);
%
for i=1:numnp
    iz=3*i;
    iy=iz-1;
    ix=iz-2;
    switch loadcode
        case 1  
    u(ix)= 0.001*(x(i)-minx)/lenx;
    u(iy)= -0.00025*(y(i)-miny)/leny;
    u(iz)= -0.00025*(z(i)-minz)/lenz;
        case 2  
    u(ix)= -0.00025*(x(i)-minx)/lenx;
    u(iy)= 0.001*(y(i)-miny)/leny;
    u(iz)= -0.0025*(z(i)-minz)/lenz;
        case 3  
    u(ix)= -0.00025*(x(i)-minx)/lenx;
    u(iy)= -0.00025*(y(i)-miny)/leny;
    u(iz)=  0.001*(z(i)-minz)/lenz;
        case 4  
    u(ix)= 0.001*(x(i)-minx)/lenx;
    u(iy)= 0.001*(y(i)-miny)/leny;
    u(iz)= -0.0005*(z(i)-minz)/lenz;
        case 5  
    u(ix)= -0.0005*(x(i)-minx)/lenx;
    u(iy)=  0.001*(y(i)-miny)/leny;
    u(iz)=  0.001*(z(i)-minz)/lenz;
        case 6  
    u(ix)=  0.001*(x(i)-minx)/lenx;
    u(iy)= -0.0005*(y(i)-miny)/leny;
    u(iz)=  0.001*(z(i)-minz)/lenz;
    end
            
end

%
