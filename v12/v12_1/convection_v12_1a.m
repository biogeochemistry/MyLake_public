% === MyLake model, version 1.2, 15.03.05 ===
% by Tom Andersen & Tuomo Saloranta, NIVA 2004

% VERSION 1.2.1a, based on convection_v12 (with three modified lines of code, marked with NEW!!!)

% Convection module
% Code checked by TSA, 07.03.05
% Last modified by TSA, 17.07.07

function [Tz,Cz,Sz,Pz,Chlz,PPz,DOPz,DOCz] = ...
    convection_v12_1a(Tz_in,Cz_in,Sz_in,Pz_in,Chlz_in,PPz_in,DOPz_in,DOCz_in,Tprof_prev,Vz,Cw,f_par,lambdaz_wtot_avg,zz,swa_b0,tracer_switch,springautumn);

% Inputs (with extension "_in") and Outputs:
%       Tz   : Temperature profile
%       Cz   : Tracer profile
%       Sz   : Suspended inorg. matter profile
%       Pz   : Dissolved inorg. P profile
%       Chlz : Chlorophyll a profile
%       PPz  : Phosphorus bound to inorganic particles profile
%       DOPz  : Dissolved organic phosphorus profile
%       DOCz  : Particulate inorganic phosphorus profile

% Inputs:
%       Tprof_prev   : Temperature profile from previous timestep
%       etc.

% These variables are still global and not transferred by functions
global ies80;

Trhomax=3.98; %temperature of maximum water density (deg C)
Nz=length(zz); %total number of layers in the water column
dz=zz(2)-zz(1); %model grid step

% Convective mixing adjustment
% Mix successive layers until stable density profile
rho = polyval(ies80,max(0,Tz_in(:)))+min(Tz_in(:),0);	% Density (kg/m3)
d_rho=[diff(rho); 1]; %d_rho = how much a layer is lighter than layer below; last cell in "d_rho" is always positive (sediment)

while any(d_rho < 0),
    blnUnstb_layers=(d_rho <= 0); %1=layer is heavier or equal than layer below, 0=layer is lighter than layer below
    A_Unstb=find(diff([0; blnUnstb_layers])==1); %layer index(es) where unstable/neutral column(s) start(s)
    B_Unstb=find(diff([0; blnUnstb_layers])==-1)-1;%layer index(es) where unstable/neutral column(s) end(s)

    for n = 1:length(A_Unstb)
        j = [A_Unstb(n):B_Unstb(n)+1];
        Tmix = sum(Tz_in(j) .* Vz(j)) / sum(Vz(j));
        Tz_in(j) = Tmix * ones(size(Tz_in(j)));

        if (tracer_switch==1)
            Cmix = sum(Cz_in(j) .* Vz(j)) / sum(Vz(j));
            Cz_in(j) = Cmix * ones(size(Cz_in(j)));
        end

        Smix = sum(Sz_in(j) .* Vz(j)) / sum(Vz(j));
        Sz_in(j) = Smix * ones(size(Sz_in(j)));

        Pmix = sum(Pz_in(j) .* Vz(j)) / sum(Vz(j));
        Pz_in(j) = Pmix * ones(size(Pz_in(j)));

        Chlmix = sum(Chlz_in(j) .* Vz(j)) / sum(Vz(j));
        Chlz_in(j) = Chlmix * ones(size(Chlz_in(j)));

        PPmix = sum(PPz_in(j) .* Vz(j)) / sum(Vz(j));
        PPz_in(j) = PPmix * ones(size(PPz_in(j)));

        DOPmix = sum(DOPz_in(j) .* Vz(j)) / sum(Vz(j));
        DOPz_in(j) = DOPmix * ones(size(DOPz_in(j)));

        DOCmix = sum(DOCz_in(j) .* Vz(j)) / sum(Vz(j));
        DOCz_in(j) = DOCmix * ones(size(DOCz_in(j)));
    end

    rho = polyval(ies80,max(0,Tz_in(:))) + min(Tz_in(:),0);
    d_rho=[diff(rho); 1];
end;

if (springautumn==1)
    % Spring/autumn turnover
    % don't allow temperature jumps over temperature of maximum density
    if( ((Tprof_prev(1)>Trhomax)&(Tz_in(1)<Trhomax))|((Tprof_prev(1)<Trhomax)&(Tz_in(1)>Trhomax)) ) %NEW!!! (Tprof, ">" instead of ">=")!

        jumpinx=find( ((Tprof_prev>Trhomax)&(Tz_in<Trhomax))|((Tprof_prev<Trhomax)&(Tz_in>Trhomax)) ); %NEW!!!!!
        if (sum(jumpinx==1)==0);disp('NOTE: Non-surface jumps over temperature of maximum density');end %NEW!!!!!

        intSign=sign(Trhomax-Tz_in(1)); %plus in autumn turnover, minus in spring turnover
        XE_turn=cumsum((Tz_in-Trhomax).*Vz*Cw*intSign); %always starts negative
        Dummy=find(XE_turn>0);
        if(isempty(Dummy)==1)
            Tz_in(:)=Trhomax;
            if (intSign==1)
                Tz_in(1)=Tz_in(1)+intSign*XE_turn(end)/(Vz(1)*Cw); %put overshoot on top layer
            else
                Tz_in=Tz_in + (-diff( intSign*XE_turn(end) * (f_par * exp([0; -lambdaz_wtot_avg] .* [zz; zz(end)+dz]) + ...
                    (1-f_par) * exp([0; -swa_b0*ones(Nz,1)] .* [zz; zz(end)+dz])) ) ./(Vz*Cw));
                %distribute overshoot as shortwave energy
            end
        else
            Tz_in(1:Dummy(1)-1)=Trhomax;
            Tz_in(Dummy(1))=Trhomax+intSign*XE_turn(Dummy(1))/(Vz(Dummy(1))*Cw);
        end
    end
end %springautumn

Tz=Tz_in;
Cz=Cz_in;
Sz=Sz_in;
Pz=Pz_in;
PPz=PPz_in;
Chlz=Chlz_in;
DOPz=DOPz_in;
DOCz=DOCz_in;