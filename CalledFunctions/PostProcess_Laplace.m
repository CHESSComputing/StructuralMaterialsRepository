function [grad_u]=PostProcess_Laplace(iele,nnpe,np,x,y,z,nqptv,sfac,dndxi,dndet,dndze,u)

% called by MechMonics

% Calculates spatial gradients of the field variable u  



detj = zeros(1,nqptv);
dndx = zeros(nnpe,nqptv);
dndy = zeros(nnpe,nqptv);
dndz = zeros(nnpe,nqptv);

grad_u = zeros(3,nqptv);

for i=1:1:nqptv  
    

tjac=zeros(3,3);
u_el=zeros(nnpe,1);

for j=1:1:nnpe
     
  j1=np(iele,j) ;  

  tjac(1,1)=tjac(1,1)+dndxi(j,i)*x(j1);                           
  tjac(1,2)=tjac(1,2)+dndet(j,i)*x(j1);	  
  tjac(1,3)=tjac(1,3)+dndze(j,i)*x(j1);                          
  tjac(2,1)=tjac(2,1)+dndxi(j,i)*y(j1);                        
  tjac(2,2)=tjac(2,2)+dndet(j,i)*y(j1);                         
  tjac(2,3)=tjac(2,3)+dndze(j,i)*y(j1); 	  
  tjac(3,1)=tjac(3,1)+dndxi(j,i)*z(j1);                          
  tjac(3,2)=tjac(3,2)+dndet(j,i)*z(j1);                          
  tjac(3,3)=tjac(3,3)+dndze(j,i)*z(j1); 
 
  % Form the elemental displacement vector
  u_el(j)=u(j1);
  
end

detj(i)=det(tjac);
tjaci = inv(tjac);

for  j=1:1:nnpe
  dndx(j,i)=tjaci(1,1)*dndxi(j,1)+tjaci(2,1)*dndet(j,1)+tjaci(3,1)*dndze(j,1);            
  dndy(j,i)=tjaci(1,2)*dndxi(j,1)+tjaci(2,2)*dndet(j,1)+tjaci(3,2)*dndze(j,1);       
  dndz(j,i)=tjaci(1,3)*dndxi(j,1)+tjaci(2,3)*dndet(j,1)+tjaci(3,3)*dndze(j,1);
end


end

% Form B-matrix

bigNprime = BigNprimeMat(nnpe,nqptv,dndx,dndy,dndz);

% Gradients

for i=1:1:nqptv  
    
grad_u(:,i) =bigNprime(:,:,i)*u_el(:);


end


