function schmid_tensors = CalcSchmidTensors(rotations,phase_schmid_tensors,phaseofgrain,maxss)
%
%  called by MechMet
%
%  Calculates the Schmid tensor for each grain in global coordinates
% 
%  converts to vector (6x1) form for tensor
%
  
  
  [mo,no,numgrains] = size(rotations);
  myrot = zeros(3,3);
  schmid_tensors = zeros(maxss,6,numgrains);
  
  for ngr = 1:numgrains

 myrot = rotations(:,:,ngr);
 
 phasenum = phaseofgrain(ngr);
 
 schmid_tensor(:,:,:) = phase_schmid_tensors(:,:,:,phasenum);

 sx_temp  = zeros(maxss,3,3);
 sti = zeros(3,3);
 sti_rotated = zeros(3,3);

 for iss = 1:1:maxss
   sti(:,:) = schmid_tensor(iss,:,:);
   sti_rotated = myrot*sti*myrot';
   sx_temp(iss,:,:) = sti_rotated(:,:);
 end

schmid_tensors(:,1,ngr) = sx_temp(:,1,1);
schmid_tensors(:,2,ngr) = sx_temp(:,2,2);
schmid_tensors(:,3,ngr) = sx_temp(:,3,3);
schmid_tensors(:,4,ngr) = 2.0*sx_temp(:,2,3);
schmid_tensors(:,5,ngr) = 2.0*sx_temp(:,1,3);
schmid_tensors(:,6,ngr) = 2.0*sx_temp(:,1,2);
    
  end
  
end

