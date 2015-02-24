% === MyLake model, version 1.1, 16.03.04 ===
% by Tom Andersen & Tuomo Saloranta, NIVA 2004
%
% Module for calculating new temperature profile in sediments
% Code checked by TSA, 16.03.2004

function [Tzy_sed_updated] = sedimentheat_v11(Tzy_sed, K_sed, dt)

% Grid resolution in sediment is 0.2 m in the 0-2 layer and 0.5 m in the 2-10 m layer 
% (totally 10 + 16 = 26 layers)
dz_s=[NaN; 0.2*ones(10-1,1); 0.35; 0.5*ones(16-1,1); NaN]; %length = 27

alpha_sed =  dt./(dz_s.^2);
[Ns,Nw]=size(Tzy_sed);

Tzy_sed_updated=zeros(Ns,Nw);

for i=1:Nw %do for every water layer in contact with sediment
    
   az = K_sed * alpha_sed(1:end-1);
   bz = K_sed * alpha_sed(2:end);
   az(1)=0; %top sediment boundary condition
   bz(1)=0; %top sediment boundary condition
   bz(end)=0; %deep sediment boundary condition
   
   Gi = [-bz (1 + az + bz) -az];
   Fi = spdiags(Gi,-1:1,Ns,Ns)';
   Tzy_sed_updated(:,i) = Fi \ Tzy_sed(:,i);
end  