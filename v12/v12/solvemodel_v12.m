% === MyLake model, version 1.2, 15.03.05 ===
% by Tom Andersen & Tuomo Saloranta, NIVA 2005
%
% Main module
% Code checked by TSA, xx.03.2005
% Last modified by TSA, 15.08.2006 (temperature profile sent to convection.m, Thermocline-bug corrected)


function [zz,Az,Vz,tt,Qst,Kzt,Tzt,Czt,Szt,Pzt,Chlzt,PPzt,DOPzt,DOCzt,Qzt_sed,lambdazt,...
         P3zt_sed,P3zt_sed_sc,His,DoF,DoM,MixStat,Wt] = ...
         solvemodel_v12(M_start,M_stop,Initfile,Initsheet,Inputfile,Inputsheet,Parafile,Parasheet,varargin);

warning off MATLAB:fzero:UndeterminedSyntax %suppressing a warning message

% Inputs (to function)
%       M_start     : Model start date [year, month, day]
%       M_stop      : Model stop date [year, month, day]
%               + Input filenames and sheetnames

% Inputs (received from input module):
%		tt		: Solution time domain (day)
%       In_Z    : Depths read from initial profiles file (m)
%       In_Az   : Areas read from initial profiles file (m2)
%       In_Tz   : Initial temperature profile read from initial profiles file (deg C)
%       In_Cz   : Initial tracer profile read from initial profiles file (-)
%       In_Sz   : Initial sedimenting tracer (or suspended inorganic matter) profile read from initial profiles file (kg m-3)
%       In_TPz  : Initial total P profile read from initial profiles file (incl. DOP & Chla) (mg m-3)
%       In_DOPz  : Initial dissolved organic P profile read from initial profiles file (mg m-3)
%       In_Chlz : Initial chlorophyll a profile read from initial profiles file (mg m-3)
%       In_DOCz  : Initial DOC profile read from initial profiles file (mg m-3)
%       In_TPz_sed  : Initial total P profile in the sediment compartments read from initial profiles file (mg m-3)
%       In_Chlz_sed : Initial chlorophyll a profile in the sediment compartments read from initial profiles file (mg m-3)
%       In_FIM      : Initial profile of volume fraction of inorganic matter in the sediment solids (dry weight basis)
%       Ice0            : Initial conditions, ice and snow thicknesses (m) (Ice, Snow)
%		Wt		        : Weather data
%       Inflow          : Inflow data
%       Phys_par        : Main 23 parameters that are more or less fixed
%       Phys_par_range  : Minimum and maximum values for Phys_par (23 * 2)
%       Phys_par_names  : Names for Phys_par
%       Bio_par         : Main 15 parameters that are more or less site specific
%       Bio_par_range   : Minimum and maximum values for Bio_par (15 * 2)
%       Bio_par_names   : Names for Bio_par

% Outputs (other than Inputs from input module):
%		Qst : Estimated surface heat fluxes ([sw, lw, sl] * tt) (W m-2)
%		Kzt	: Predicted vertical diffusion coefficient (tt * zz) (m2 d-1)
%		Tzt	: Predicted temperature profile (tt * zz) (deg C)
%		Czt	: Predicted passive tracer profile (tt * zz) (-)
%		Szt	: Predicted passive sedimenting tracer (or suspended inorganic matter) profile (tt * zz) (kg m-3)=(g L-1)
%		Pzt	: Predicted dissolved inorganic phosphorus profile (tt * zz) (mg m-3)
%		Chlzt	    : Predicted chlorophyll a profile (tt * zz) (mg m-3)
%		PPzt	    : Predicted particulate inorganic phosphorus profile (tt * zz) (mg m-3)
%		DOPzt	    : Predicted dissolved organic phosphorus profile (tt * zz) (mg m-3)
%		DOCzt	    : Predicted dissolved organic carbon (DOC) profile (tt * zz) (mg m-3)
%		Qz_sed      : Predicted  sediment-water heat flux (tt * zz) (W m-2, normalised to lake surface area)
%       lambdazt    : Predicted average total light attenuation coefficient down to depth z (tt * zz) (m-1)
%       P3zt_sed    : Predicted P conc. in sediment for P (mg m-3), PP(mg kg-1 dry w.) and Chl (mg kg-1 dry w.) (tt * zz * 3)
%       P3zt_sed_sc : Predicted P source from sediment for P, PP and Chl (mg m-3 day-1) (tt * zz * 3)
%       His         : Ice information matrix ([Hi Hs Hsi Tice Tair rho_snow IceIndicator] * tt)
%       DoF, DoM    : Days of freezing and melting (model timestep number)
%       MixStat     : Temporary variables used in model testing, see code (N * tt)

% These variables are still global and not transferred by functions
global ies80;

tic
disp(['Running MyLake from ' datestr(datenum(M_start)) ' to ' datestr(datenum(M_stop)) ' ...']);

% ===Switches===
snow_compaction_switch=1;       %snow compaction: 0=no, 1=yes
river_inflow_switch=1;          %river inflow: 0=no, 1=yes
sediment_heatflux_switch=0;     %heatflux from sediments: 0=no, 1=yes
selfshading_switch=1;           %light attenuation by chlorophyll a: 0=no, 1=yes
tracer_switch=1;               %simulate tracers:  0=no, 1=yes
% ==============

dt=1.0; %model time step = 1 day (DO NOT CHANGE!)

if (nargin>8) %if optional command line parameter input is used 
  disp('Bypassing input files...Running with input data & parameters given on command line');
 [In_Z,In_Az,tt,In_Tz,In_Cz,In_Sz,In_TPz,In_DOPz,In_Chlz,In_DOCz,In_TPz_sed,In_Chlz_sed,In_FIM,Ice0,Wt,Inflw,...
         Phys_par,Phys_par_range,Phys_par_names,Bio_par,Bio_par_range,Bio_par_names]...
    = deal(varargin{:});
else
 %Read input data
  [In_Z,In_Az,tt,In_Tz,In_Cz,In_Sz,In_TPz,In_DOPz,In_Chlz,In_DOCz,In_TPz_sed,In_Chlz_sed,In_FIM,Ice0,Wt,Inflw,...
         Phys_par,Phys_par_range,Phys_par_names,Bio_par,Bio_par_range,Bio_par_names]...
    = modelinputs_v12(M_start,M_stop,Initfile,Initsheet,Inputfile,Inputsheet,Parafile,Parasheet,dt);
end
            
%load H:\MyLake\air_sea\albedot1.mat; %load albedot1 table, in order to save execution time
load albedot1.mat; %load albedot1 table, in order to save execution time

% Unpack the more fixed parameter values from input array "Phys_par"
dz = Phys_par(1); %grid stepsize (m)

   zm = In_Z(end); %max depth
   zz = [0:dz:zm-dz]'; %solution depth domain

Kz_K1 = Phys_par(2); % open water diffusion parameter (-)
Kz_K1_ice = Phys_par(3); % under ice diffusion parameter (-)
Kz_N0 = Phys_par(4); % min. stability frequency (s-2)
C_shelter = Phys_par(5); % wind shelter parameter (-)
lat = Phys_par(6); %latitude (decimal degrees)
lon = Phys_par(7); %longitude (decimal degrees)
alb_melt_ice = Phys_par(8);   %albedo of melting ice (-)
alb_melt_snow = Phys_par(9); %albedo of melting snow (-)
PAR_sat = Phys_par(10);         %PAR saturation level for phytoplankton growth (mol(quanta) m-2 s-1) 
f_par = Phys_par(11);           %Fraction of PAR in incoming solar radiation (-)
beta_chl = Phys_par(12);        %Optical cross_section of chlorophyll (m2 mg-1)
lambda_i = Phys_par(13);       %PAR light attenuation coefficient for ice (m-1)
lambda_s = Phys_par(14);       %PAR light attenuation coefficient for snow (m-1)
F_sed_sld = Phys_par(15);      %volume fraction of solids in sediment (= 1-porosity)
I_scV = Phys_par(16); %scaling factor for inflow volume (-)
I_scT = Phys_par(17); %scaling coefficient for inflow temperature (-) 
I_scC = Phys_par(18); %scaling factor for inflow concentration of C (-)
I_scS = Phys_par(19); %scaling factor for inflow concentration of S (-)
I_scTP = Phys_par(20); %scaling factor for inflow concentration of total P (-)
I_scDOP = Phys_par(21); %scaling factor for inflow concentration of diss. organic P (-)
I_scChl = Phys_par(22); %scaling factor for inflow concentration of Chl a (-)
I_scDOC = Phys_par(23); %scaling factor for inflow concentration of DOC  (-)


% Unpack the more site specific parameter values from input array "Bio_par"

swa_b0 = Bio_par(1); % non-PAR light atteneuation coeff. (m-1)
swa_b1 = Bio_par(2); %  PAR light atteneuation coeff. (m-1)
S_res_epi = Bio_par(3);      %Particle resuspension mass transfer coefficient, epilimnion (m day-1, dry)
S_res_hypo = Bio_par(4);     %Particle resuspension mass transfer coefficient, hypolimnion (m day-1, dry)
H_sed = Bio_par(5);          %height of active sediment layer (m, wet mass)
Psat_L = Bio_par(6);           %Half saturation parameter for Langmuir isotherm
Fmax_L = Bio_par(7);           %Scaling parameter for Langmuir isotherm
w_s = Bio_par(8);            %settling velocity for S (m day-1)
w_chl = Bio_par(9);          %settling velocity for Chl a (m day-1)
Y_cp = Bio_par(10);            %yield coefficient (chlorophyll to carbon) * (carbon to phosphorus) ratio (-)
m_twty = Bio_par(11);          %loss rate (1/day) at 20 deg C
g_twty = Bio_par(12);          %specific growth rate (1/day) at 20 deg C
k_twty = Bio_par(13);          %specific Chl a to P transformation rate (1/day) at 20 deg C
dop_twty = Bio_par(14);        %specific DOP to P transformation rate (day-1) at 20 deg C
P_half = Bio_par(15);          %Half saturation growth P level (mg/m3)


% Other variables/parameters not read from the input file:

Nz=length(zz); %total number of layers in the water column
N_sed=26; %total number of layers in the sediment column

theta_m = exp(0.1*log(2));    %loss and growth rate parameter base, ~1.072
e_par = 240800;               %Average energy of PAR photons (J mol-1)

% diffusion parameterisation exponents
Kz_b1 = 0.43;
Kz_b1_ice  = 0.43; 

% ice & snow parameter values
rho_fw=1000;        %density of freshwater (kg m-3)
rho_ice=910;        %ice (incl. snow ice) density (kg m-3)
rho_new_snow=250;   %new-snow density (kg m-3)
max_rho_snow=450;   %maximum snow density (kg m-3)
L_ice=333500;       %latent heat of freezing (J kg-1)
K_ice=2.1;          %ice heat conduction coefficient (W m-1 K-1)
C1=7.0;             %snow compaction coefficient #1
C2=21.0;            %snow compaction coefficient #2

Tf=0;               %water freezing point temperature (deg C)

F_OM=1e+6*0.012;    %mass fraction [mg kg-1] of P of dry organic matter (assuming 50% of C, and Redfield ratio)

K_sed=0.035;      %thermal diffusivity of the sediments (m2 day-1) 
rho_sed=2500;      %bulk density of the inorganic solids in sediments (kg m-3)
rho_org=1000;      %bulk density of the organic solids in sediments (kg m-3)
cp_sed=1000;       %specific heat capasity of the sediments (J kg-1 K-1)
%=======
%================NEW testing NEW testing NEW testing NEW testing NEW testing 
ksw=1e-3; %(m/d) Porewater_water mass transfer coefficient; NEW testing 3.8.05
Fmax_L_sed=1*Fmax_L; %Fmax for sediment; NEW, testing 3.8.05
Fstable=655; % (mg/kg) Inactive P conc. in inorg. particles; NEW, testing 3.8.05            

% Allocate and initialise output data matrices
Qst = zeros(3,length(tt));
Kzt = zeros(Nz,length(tt));
Tzt = zeros(Nz,length(tt));
Czt = zeros(Nz,length(tt));
Szt = zeros(Nz,length(tt));
Pzt = zeros(Nz,length(tt));
Chlzt = zeros(Nz,length(tt));
PPzt = zeros(Nz,length(tt));
DOPzt = zeros(Nz,length(tt));
DOCzt = zeros(Nz,length(tt));
Qzt_sed = zeros(Nz,length(tt));
lambdazt = zeros(Nz,length(tt));
P3zt_sed = zeros(Nz,length(tt),4); %3-D
P3zt_sed_sc = zeros(Nz,length(tt),3); %3-D
His = zeros(7,length(tt)); 
MixStat = zeros(20,length(tt)); 

% Initial profiles

  Az = interp1(In_Z,In_Az,zz);
  Vz = dz * (Az + [Az(2:end); 0]) / 2;

   T0 = interp1(In_Z,In_Tz,zz+dz/2); % Initial temperature distribution (deg C)
   C0 = interp1(In_Z,In_Cz,zz+dz/2); % Initial passive tracer distribution (-)
   S0 = interp1(In_Z,In_Sz,zz+dz/2); % Initial passive sedimenting tracer (or suspended inorganic matter) distribution (kg m-3)
   TP0 = interp1(In_Z,In_TPz,zz+dz/2);	% Initial total P distribution (incl. DOP & Chla) (mg m-3)
   DOP0 = interp1(In_Z,In_DOPz,zz+dz/2);	% Initial dissolved organic P distribution (mg m-3)
   Chl0 = interp1(In_Z,In_Chlz,zz+dz/2);	% Initial chlorophyll a distribution (mg m-3)
   DOC0 = interp1(In_Z,In_DOCz,zz+dz/2);	% Initial DOC distribution (mg m-3)
   TP0_sed = interp1(In_Z,In_TPz_sed,zz+dz/2); % Initial total P distribution in bulk wet sediment ((mg m-3); particles + porewater)
   Chl0_sed = interp1(In_Z,In_Chlz_sed,zz+dz/2); % Initial chlorophyll a distribution in bulk wet sediment (mg m-3)
   FIM0 = interp1(In_Z,In_FIM,zz+dz/2);     % Initial sediment solids volume fraction of inorganic matter (-)

   VolFrac=1./(1+(1-F_sed_sld)./(F_sed_sld*FIM0)); %volume fraction: inorg sed. / (inorg.sed + pore water)

   if any((TP0-DOP0-Chl0./Y_cp-S0*Fstable)<0)
       error('Sum of initial DOP, stably particle bound P, and P contained in Chl a cannot be larger than TP')
   end
   
   if any((TP0_sed-DOP0-Chl0_sed./Y_cp-VolFrac*rho_sed*Fstable)<0)
       error('Sum of initial DOP stably, particle bound P, and P contained in Chl_sed a cannot be larger than TP_sed')
   end
      
   if (any(FIM0<0)|any(FIM0>1))
       error('Initial fraction of inorganic matter in sediments must be between 0 and 1')
   end

   if (any(ksw>(H_sed*(1-F_sed_sld))))
       error('Parameter ksw is larger than the volume (thickness) of porewater')
   end  %OBS! Ideally should also be that the daily diffused porewater should not be larger 
        %than the corresponding water layer volume, but this seems very unlike in practise

Tz = T0;
Cz = C0;
Sz = S0; % (kg m-3)
Chlz = Chl0;  % (mg m-3)
DOPz = DOP0;  % (mg m-3)
[Pz, trash] =Ppart(S0./rho_sed,TP0-DOP0-(Chl0./Y_cp),Psat_L,Fmax_L,rho_sed,Fstable);  % (mg m-3)
PPz = TP0-DOP0-(Chl0./Y_cp)-Pz; % (mg m-3) 
DOCz = DOC0;   % (mg m-3)
F_IM = FIM0; %initial VOLUME fraction of inorganic particles of total dry sediment solids

%== P-partitioning in sediments==
%Pdz_store: %diss. inorg. P in sediment pore water (mg m-3)
%Psz_store: %P conc. in inorganic sediment particles (mg kg-1 dry w.)

[Pdz_store, Psz_store]=Ppart(VolFrac,TP0_sed-(Chl0_sed./Y_cp)-DOP0,Psat_L,Fmax_L_sed,rho_sed,Fstable);

%Chlsz_store: %Chla conc. in organic sediment particles (mg kg-1 dry w.)
 Chlsz_store = Chl0_sed./(rho_org*F_sed_sld*(1-F_IM)); %(mg kg-1 dry w.)

 
% assume linear initial temperature profile in sediment (4 deg C at the bottom)
clear Tzy_sed
 for j=1:Nz
  Tzy_sed(:,j) = interp1([0.2 10], [Tz(j) 4], [0.2:0.2:2 2.5:0.5:10])';
 end

S_resusp=S_res_hypo*ones(Nz,1); %hypolimnion resuspension assumed on the first time step

rho_snow=rho_new_snow;   %initial snow density (kg m-3)
Tice=NaN;                %ice surface temperature (initial value, deg C)
XE_melt=0;               %energy flux that is left from last ice melting (initial value, W m-2)
XE_surf=0;               %energy flux from water to ice (initial value,  J m-2 per day)

%Initialisation of ice & snow variables
Hi=Ice0(1);               %total ice thickness (initial value, m)
WEQs=(rho_snow/rho_fw)*Ice0(2); %snow water equivalent  (initial value, m)
Hsi=0;                %snow ice thickness (initial value = 0 m)

if ((Hi<=0)&(WEQs>0))
    error('Mismatch in initial ice and snow thicknesses')    
end

if (Hi<=0)
    IceIndicator=0;     %IceIndicator==0 means no ice cover
else
    IceIndicator=1;
end

pp=1; %initial indexes for ice freezing/melting date arrays
qq=1;
DoF=[]; %initialize
DoM=[]; %initialize

% >>>>>> Start of the time loop >>>>>>
Resuspension_counter=zeros(Nz,1); %kg
Sedimentation_counter=zeros(Nz,1); %kg
SS_decr=0; %kg
for i = 1:length(tt)

    Pz_prev=Pz; PPz_prev=PPz; DOPz_prev=DOPz; Chlz_prev=Chlz; %for monitoring purposes
    
 % Surface heat fluxes (W m-2), wind stress (N m-2) & daylight fraction (-), based on Air-Sea Toolbox
 [Qsw,Qlw,Qsl,tau,DayFrac,DayFracHeating] = heatflux_v12(tt(i),Wt(i,1),Wt(i,2),Wt(i,3),Wt(i,4),Wt(i,5),Wt(i,6),Tz(1), ...
     lat,lon,WEQs,Hi,alb_melt_ice,alb_melt_snow,albedot1);     %Qlw and Qsl are functions of Tz(1)            

 % Calculate total mean PAR and non-PAR light extinction coefficient in water (background + due to Chl a)
 lambdaz_wtot_avg=zeros(Nz,1);
 lambdaz_NP_wtot_avg=zeros(Nz,1);
 
 if (selfshading_switch==1)
  lambdaz_wtot=swa_b1 * ones(Nz,1) + beta_chl*Chlz; %at layer z
  lambdaz_NP_wtot=swa_b0 * ones(Nz,1) + beta_chl*Chlz; %at layer z
  for j=1:Nz
     lambdaz_wtot_avg(j)=mean(swa_b1 * ones(j,1) + beta_chl*Chlz(1:j)); %average down to layer z 
     lambdaz_NP_wtot_avg(j)=mean(swa_b0 * ones(j,1) + beta_chl*Chlz(1:j)); %average down to layer z     
  end
 else %constant with depth
  lambdaz_wtot=swa_b1 * ones(Nz,1); 
  lambdaz_wtot_avg=swa_b1 * ones(Nz,1); 
  lambdaz_NP_wtot=swa_b0 * ones(Nz,1); 
  lambdaz_NP_wtot_avg=swa_b0 * ones(Nz,1); 
 end %if selfshading...
    
%Photosynthetically Active Radiation                                
H_sw_z=NaN*zeros(Nz,1);

if(IceIndicator==0)
  IceSnowAttCoeff=1; %no extra light attenuation due to snow and ice
  PAR_z=((3/2) / (e_par * DayFrac)) * f_par * Qsw  * exp(-lambdaz_wtot_avg .* zz);   
  %Irradiance at noon (mol m-2 s-1) at levels zz 
else    %extra light attenuation due to ice and snow
  IceSnowAttCoeff=exp(-lambda_i * Hi) * exp(-lambda_s * (rho_fw/rho_snow)*WEQs);  
  PAR_z=((3/2) / (e_par * DayFrac)) * IceSnowAttCoeff * f_par *...
        Qsw  * exp(-lambdaz_wtot_avg .* zz); 
end

U_sw_z=PAR_z./PAR_sat; %scaled irradiance at levels zz                                                             
inx_u=find(U_sw_z<=1); %undersaturated
inx_s=find(U_sw_z>1);  %saturated

H_sw_z(inx_u)=(2/3)*U_sw_z(inx_u);  %undersaturated
 
  dum_a=sqrt(U_sw_z);
  dum_b=sqrt(U_sw_z-1);
H_sw_z(inx_s)=(2/3)*U_sw_z(inx_s) + log((dum_a(inx_s) + dum_b(inx_s))./(dum_a(inx_s) ...  %saturated
            - dum_b(inx_s))) - (2/3)*(U_sw_z(inx_s)+2).*(dum_b(inx_s)./dum_a(inx_s));
        
% Update biological rates at layers z
Growth_bioz=g_twty*theta_m.^(Tz-20) .* (Pz./(P_half+Pz)) .* (DayFrac./(dz*lambdaz_wtot)) .* diff([-H_sw_z; 0]); 
Loss_bioz=m_twty*theta_m.^(Tz-20);
R_bioz = min((Growth_bioz-Loss_bioz),(Y_cp*Pz./(Chlz*dt))); %growth rate is limited by available phosphorus 

   
   Tprof_prev=Tz; %temperature profile at previous time step (for modified convection.m)
   rho = polyval(ies80,max(0,Tz(:))) + min(Tz(:),0);  % Density (kg/m3)                                                      

   % Sediment vertical heat flux, Q_sed 
   % (averaged over the whole top area of the layer, although actually coming only from the "sides")   
    if (sediment_heatflux_switch==1)
     % update top sediment temperatures
     dz_sf = 0.2; %fixed distance between the two topmost sediment layers (m)
     Tzy_sed(1,:) = Tz';    
     Tzy_sed_upd = sedimentheat_v11(Tzy_sed, K_sed, dt);
     Tzy_sed=Tzy_sed_upd;
     Qz_sed=K_sed*rho_sed*cp_sed*(1/dz_sf)*(-diff([Az; 0])./Az) .* (Tzy_sed(2,:)'-Tzy_sed(1,:)'); %(J day-1 m-2)
        %positive heat flux => from sediment to water
    else
     Qz_sed = zeros(Nz,1);
    end
     
   Cw = 4.18e+6;	% Volumetric heat capacity of water (J K-1 m-3)
     
   %Heat sources/sinks:  
   %Total attenuation coefficient profile, two-band extinction, PAR & non-PAR
   Par_Attn=exp([0; -lambdaz_wtot_avg] .* [zz; zz(end)+dz]);
   NonPar_Attn=exp([0; -lambdaz_NP_wtot_avg] .* [zz; zz(end)+dz]);
   Attn_z=(f_par * (Par_Attn(1:end-1) - ([Az(2:end);0]./Az).*Par_Attn(2:end)) + ...
       (1-f_par) * (NonPar_Attn(1:end-1) - ([Az(2:end);0]./Az).*NonPar_Attn(2:end))); 

   
   if(IceIndicator==0)
     % 1) Vertical heating profile for open water periods (during daytime heating)
    Qz = (Qsw + XE_melt) * Attn_z; %(W m-2)
    Qz(1) = Qz(1) + DayFracHeating*(Qlw + Qsl); %surface layer heating    
    XE_melt=0; %Reset    
    dT = Az .* ((60*60*24*dt) * Qz + DayFracHeating*Qz_sed) ./ (Cw * Vz); %Heat source (K day-1) (daytime heating, ice melt, sediment);
   else
     % Vertical heating profile for ice-covered periods (both day- and nighttime)
    Qz = Qsw * IceSnowAttCoeff * Attn_z; %(W/m2)
    dT = Az .* ((60*60*24*dt) * Qz + Qz_sed) ./ (Cw * Vz); %Heat source (K day-1) (solar rad., sediment);          
   end
   
   Tz = Tz + dT;        %Temperature change after daytime surface heatfluxes (or whole day in ice covered period)

% Convective mixing adjustment (mix successive layers until stable density profile)
% and 
% Spring/autumn turnover (don't allow temperature jumps over temperature of maximum density)
[Tz,Cz,Sz,Pz,Chlz,PPz,DOPz,DOCz] = convection_v12_1a(Tz,Cz,Sz,Pz,Chlz,PPz,DOPz,DOCz,Tprof_prev,Vz,Cw,f_par,lambdaz_wtot_avg,zz,swa_b0,tracer_switch,1);

   if(IceIndicator==0)
     % 2) Vertical heating profile for open water periods (during nighttime heating)
    [Qsw,Qlw_2,Qsl_2,tau,DayFrac,DayFracHeating] = heatflux_v12(tt(i),Wt(i,1),Wt(i,2),Wt(i,3),Wt(i,4),Wt(i,5),Wt(i,6),Tz(1), ...
     lat,lon,WEQs,Hi,alb_melt_ice,alb_melt_snow,albedot1); %Qlw and Qsl are functions of Tz(1)            
    Qz(1) = (1-DayFracHeating)*(Qlw_2 + Qsl_2); %surface layer heating
    Qz(2:end)=0; %No other heating below surface layer   
    dT = Az .* ((60*60*24*dt) * Qz + (1-DayFracHeating)*Qz_sed) ./ (Cw * Vz); %Heat source (K day-1) (longwave & turbulent fluxes);
   
    Tz = Tz + dT;         %Temperature change after nighttime surface heatfluxes

   % Convective mixing adjustment (mix successive layers until stable density profile)  
   % and 
   % Spring/autumn turnover (don't allow temperature jumps over temperature of maximum density)
   [Tz,Cz,Sz,Pz,Chlz,PPz,DOPz,DOCz] = convection_v12_1a(Tz,Cz,Sz,Pz,Chlz,PPz,DOPz,DOCz,Tprof_prev,Vz,Cw,f_par,lambdaz_wtot_avg,zz,swa_b0,tracer_switch,1);

   Qlw = DayFracHeating*Qlw + (1-DayFracHeating)*Qlw_2; %total amounts, only for output purposes
   Qsl = DayFracHeating*Qsl + (1-DayFracHeating)*Qsl_2; %total amounts, only for output purposes
  end
      
   % Vertical turbulent diffusion
   g   = 9.81;							% Gravity acceleration (m s-2)
   rho = polyval(ies80,max(0,Tz(:))) + min(Tz(:),0);  % Water density (kg m-3)                                                      
    % Note: in equations of rho it is assumed that every supercooled degree lowers density by 
    % 1 kg m-3 due to frazil ice formation (probably no practical meaning, but included for "safety")
    
   N2  = g * (diff(log(rho)) ./ diff(zz));	% Brunt-Vaisala frequency (s-2) for level (zz+1)              
   if (IceIndicator==0)
     Kz  = Kz_K1 * max(Kz_N0, N2).^(-Kz_b1);	% Vertical diffusion coeff. in ice free season (m2 day-1)
                                                % for level (zz+1)              
   else 
     Kz  = Kz_K1_ice * max(Kz_N0, N2).^(-Kz_b1_ice); % Vertical diffusion coeff. under ice cover (m2 day-1)
                                                     % for level (zz+1)              
   end
   
   Fi = tridiag_DIF_v11([NaN; Kz],Vz,Az,dz,dt); %Tridiagonal matrix for general diffusion

   Tz = Fi \ (Tz);        %Solving new temperature profile (diffusion, sources/sinks already added to Tz above)

 % Convective mixing adjustment (mix successive layers until stable density profile)  
 % and 
 % Spring/autumn turnover (don't allow temperature jumps over temperature of maximum density)
 [Tz,Cz,Sz,Pz,Chlz,PPz,DOPz,DOCz] = convection_v12_1a(Tz,Cz,Sz,Pz,Chlz,PPz,DOPz,DOCz,Tprof_prev,Vz,Cw,f_par,lambdaz_wtot_avg,zz,swa_b0,tracer_switch,1);
       
   %Tracers
  if (tracer_switch==1)
   Cz = Fi \ Cz;             %Solving new passive tracer profile (diffusion)
  end 
  
   dDOP =  dop_twty * DOPz .* theta_m.^(Tz-20);  %Mineralisation to P
   DOPz = Fi \ (DOPz - dDOP);         %Solving new dissolved inorganic P profile (diffusion)
  
   % Suspended solids, particulate inorganic P 
   Fi_ad = tridiag_HAD_v11([NaN; Kz],w_s,Vz,Az,dz,dt); %Tridiagonal matrix for advection and diffusion 
       
   dSz_inorg = rho_sed*S_resusp.*F_IM.*(-diff([Az; 0])./Vz);  % Dry inorganic particle resuspension source from sediment (kg m-3 day-1)
   Sz = Fi_ad \ (Sz + dSz_inorg);           %Solving new suspended solids profile (advection + diffusion)    

   dPP = dSz_inorg.*Psz_store;  % PP resuspension source from sediment((kg m-3 day-1)*(mg kg-1) = mg m-3 day-1) 
   PPz = Fi_ad \ (PPz + dPP);     %Solving new suspended particulate inorganic P profile (advection + diffusion)
   
   %Chlorophyll a
   dSz_org = rho_org*S_resusp.*(1-F_IM).*(-diff([Az; 0])./Vz);  %Dry organic particle resuspension source from sediment (kg m-3 day-1)
   dChl_res = dSz_org.*Chlsz_store;  %Chl a resuspension source from sediment resusp. ((kg m-3 day-1)*(mg kg-1) = mg m-3 day-1);  
   dChl_growth = Chlz .* R_bioz; %Chl a growth source 
   dChl = dChl_growth +  dChl_res; % Total Chl a source
   Fi_ad = tridiag_HAD_v11([NaN; Kz],w_chl,Vz,Az,dz,dt); %Tridiagonal matrix for advection and diffusion    
   Chlz = Fi_ad \ (Chlz + dChl);  %Solving new phytoplankton profile (advection + diffusion) (always larger than background level)
   
   %Dissolved inorganic phosphorus
   dP = dDOP - dChl_growth./ Y_cp; %DOP source, P sink = Chla growth source
   Pz = Fi \ (Pz + dP); %Solving new dissolved inorganic P profile (diffusion)
       

   %Dissolved organic carbon
   dDOC = 0; %no sources/sinks yet implemented
   DOCz = Fi \ (DOCz + dDOC); %Solving new dissolved inorganic P profile (diffusion)
   
   %Sediment-water exchange (DOP source missing, neglect it????)
   %-porewater to water

   PwwFrac=ksw*(-diff([Az; 0]))./Vz; %fraction between resuspended porewater and water layer volumes
   %PwwFrac=(((1-F_sed_sld)/F_sed_sld)*S_resusp.*(-diff([Az; 0]))./Vz); %fraction between resuspended porewater and water layer volumes
   EquP1 = (1-PwwFrac).*Pz + PwwFrac.*Pdz_store; %Mixture of porewater and water 
   dPW_up = EquP1-Pz; %"source/sink" for output purposes
    
   %-water to porewater 
   PwwFrac=ksw./((1-F_sed_sld)*H_sed); %NEW testing 3.8.05; fraction between resuspended (incoming) water and sediment layer volumes
   %PwwFrac=S_resusp./(F_sed_sld*H_sed); %fraction between resuspended (incoming) water and sediment layer volumes
   EquP2 = PwwFrac.*Pz + (1-PwwFrac).*Pdz_store; %Mixture of porewater and water 
   dPW_down = EquP2-Pdz_store; %"source/sink" for output purposes
   
   %-update concentrations
   Pz = EquP1;
   Pdz_store=EquP2;


    %Calculate the thickness ratio of newly settled net sedimentation and mix these
    %two to get new sediment P concentrations in sediment (taking into account particle resuspension) 
    delPP_inorg=NaN*ones(Nz,1); %initialize
    delC_inorg=NaN*ones(Nz,1); %initialize
    delC_org=NaN*ones(Nz,1); %initialize
    
    delA=diff([Az; 0]); %Area difference for layer i (OBS: negative)
    meanA=0.5*(Az+[Az(2:end); 0]);

    
    %sedimentation is calculated from "Funnelling-NonFunnelling" difference
    %(corrected 03.10.05)
    delPP_inorg(1)=(0 - PPz(1)*delA(1)./meanA(1))./(dz/(dt*w_s) + 1);
    delC_inorg(1)=(0 - Sz(1)*delA(1)./meanA(1))./(dz/(dt*w_s) + 1);
    delC_org(1)=(0 - Chlz(1)*delA(1)./meanA(1))./(dz/(dt*w_chl) + 1);
    
    for ii=2:Nz
     delPP_inorg(ii)=(delPP_inorg(ii-1) - PPz(ii)*delA(ii)./meanA(ii))./(dz/(dt*w_s) + 1); %(mg m-3)
     delC_inorg(ii)=(delC_inorg(ii-1) - Sz(ii)*delA(ii)./meanA(ii))./(dz/(dt*w_s) + 1); %(kg m-3)
     delC_org(ii)=(delC_org(ii-1) - Chlz(ii)*delA(ii)./meanA(ii))./(dz/(dt*w_chl) + 1); %(mg m-3)
    end
    
    H_netsed_inorg=max(0, (Vz./(-diff([Az; 0]))).*delC_inorg./rho_sed - F_IM.*S_resusp); %inorganic(m day-1, dry), always positive
    H_netsed_org=max(0, (Vz./(-diff([Az; 0]))).*delC_org./(F_OM*Y_cp*rho_org) - (1-F_IM).*S_resusp); %organic (m day-1, dry), always positive
    
%     H_netsed_inorg=max(0, w_s*Sz./rho_sed - F_IM.*S_resusp); %inorganic(m day-1, dry), always positive
%     H_netsed_org=max(0, w_chl*Chlz./(F_OM*rho_org) - (1-F_IM).*S_resusp); %organic (m day-1, dry), always positive
    H_totsed=H_netsed_org + H_netsed_inorg;  %total (m day-1), always positive
  
    F_IM_NewSed=F_IM;    
    inx=find(H_totsed>0);
    F_IM_NewSed(inx)=H_netsed_inorg(inx)./H_totsed(inx); %volume fraction of inorganic matter in net settling sediment
    
    NewSedFrac = min(1, H_totsed./(F_sed_sld*H_sed)); %Fraction of newly fallen net sediment of total active sediment depth, never above 1
    NewSedFrac_inorg = min(1, H_netsed_inorg./(F_IM.*F_sed_sld*H_sed)); %Fraction of newly fallen net inorganic sediment of total active sediment depth, never above 1
    NewSedFrac_org = min(1, H_netsed_org./((1-F_IM).*F_sed_sld*H_sed)); %Fraction of newly fallen net organic sediment of total active sediment depth, never above 1
    
    Psz_store_prev=Psz_store; %For monitoring purposes
    Chlsz_store_prev=Chlsz_store; %For monitoring purposes
    
    %Psz_store: %P conc. in inorganic sediment particles (mg kg-1 dry w.)
    Psz_store = (1-NewSedFrac_inorg).*Psz_store + NewSedFrac_inorg.*PPz./Sz; %(mg kg-1)  

    %Update counters
    Sedimentation_counter = Sedimentation_counter + Vz.*(delC_inorg + delC_org./(F_OM*Y_cp)); %Inorg.+Org. (kg)
    Resuspension_counter = Resuspension_counter + Vz.*(dSz_inorg + dSz_org); %Inorg.+Org. (kg) 

    
    %Chlsz_store: %Chl a conc. in sediment particles (mg kg-1 dry w.)
    Chlsz_store = (1-NewSedFrac_org).*Chlsz_store + NewSedFrac_org.*F_OM*Y_cp; %(mg kg-1)     
    %Subtract degradation to P in pore water
    Chlz_seddeg = k_twty * Chlsz_store .* theta_m.^(Tz-20);
    Chlsz_store = Chlsz_store - Chlz_seddeg;
    Pdz_store=Pdz_store + Chlz_seddeg .* (rho_org*F_sed_sld*(1-F_IM))./Y_cp;
   
    %== P-partitioning in sediments==   
    VolFrac=1./(1+(1-F_sed_sld)./(F_sed_sld*F_IM)); %volume fraction: inorg sed. / (inorg.sed + pore water)
    TIP_sed =rho_sed*VolFrac.*Psz_store + (1-VolFrac).*Pdz_store; %total inorganic P in sediments (mg m-3) 
    [Pdz_store, Psz_store]=Ppart(VolFrac,TIP_sed,Psat_L,Fmax_L_sed,rho_sed,Fstable);
    %calculate new VOLUME fraction of inorganic particles of total dry sediment
    F_IM=min(1,((k_twty *(1-F_IM).*theta_m.^(Tz-20)) + F_IM)).*(1-NewSedFrac) + F_IM_NewSed.*NewSedFrac; 
    
    
  % Inflow calculation
% Inflw(:,1) Inflow volume (m3 day-1)
% Inflw(:,2) Inflow temperature (deg C)
% Inflw(:,3) Inflow tracer concentration (-)
% Inflw(:,4) Inflow sedimenting tracer (or suspended inorganic matter) concentration (kg m-3)
% Inflw(:,5) Inflow total phosphorus (TP) concentration  (incl. DOP & Chla) (mg m-3)
% Inflw(:,6) Inflow dissolved organic phosphorus (DOP) concentration (mg m-3)
% Inflw(:,7) Inflow chlorophyll a concentration (mg m-3)
% Inflw(:,8) Inflow DOC concentration (mg m-3)

  if (river_inflow_switch==1)
   Iflw = I_scV * Inflw(i,1); % (scaled) inflow rate
   Iflw_T = I_scT + Inflw(i,2); %(adjusted) inflow temperature
      if (Iflw_T<Tf) %negative temperatures changed to Tf
         Iflw_T=Tf;
      end
   Iflw_C = I_scC * Inflw(i,3); %(scaled) inflow C concentration
   Iflw_S = I_scS * Inflw(i,4); %(scaled) inflow S concentration
   Iflw_TP = I_scTP * Inflw(i,5); %(scaled) inflow TP concentration (incl. DOP & Chla)
   Iflw_DOP = I_scDOP * Inflw(i,6); %(scaled) inflow DOP concentration
   Iflw_Chl = I_scChl * Inflw(i,7); %(scaled) inflow Chl a concentration
   Iflw_DOC = I_scDOC * Inflw(i,8); %(scaled) inflow DOC concentration
 
       
%       if any((1-(Iflw_DOP+Iflw_Chl./Y_cp)./Iflw_TP-(Iflw_S*Fstable)./Iflw_TP)<0.1); %If Rbio<0.1 (bioavailability factor)
%           Iflw_S_dum = (1 - (Iflw_DOP+Iflw_Chl./Y_cp)./Iflw_TP - 0.1).*(Iflw_TP./Fstable); %Truncate Iflw_S!
%           SS_decr=SS_decr+(Iflw_S-Iflw_S_dum)*Iflw;
%           Iflw_S=Iflw_S_dum;
%       end

      if any((Iflw_TP-Iflw_DOP-Iflw_Chl./Y_cp-Iflw_S*Fstable)<0)  
         error('Sum of DOP, inactive PP, and P contained in Chl a in inflow cannot be larger than TP')
      end
      
      
   if(Iflw>0)
       if (isnan(Iflw_T))
        lvlD=0;
        Iflw_T=Tz(1);
       else
        rho = polyval(ies80,max(0,Tz(:)))+min(Tz(:),0);	% Density (kg/m3)
        rho_Iflw=polyval(ies80,max(0,Iflw_T))+min(Iflw_T,0);
        lvlG=find(rho>=rho_Iflw);
        if (isempty(lvlG))
            lvlG=length(rho);
        end
        lvlD=zz(lvlG(1)); %level zz above which inflow is put
       end %if isnan...
       
       
     %Changes in properties due to inflow 
    dummy=IOflow_v11(dz, zz, Vz, Tz, lvlD, Iflw, Iflw_T);
    Tz=dummy; %Temperature

     if (tracer_switch==1)    
      dummy=IOflow_v11(dz, zz, Vz, Cz, lvlD, Iflw, Iflw_C);
      Cz=dummy; %Passive tracer
     end
     
    dummy=IOflow_v11(dz, zz, Vz, Sz, lvlD, Iflw, Iflw_S);
    Sz=dummy; %Sedimenting tracer
     
    dummy=IOflow_v11(dz, zz, Vz, DOPz, lvlD, Iflw, Iflw_DOP);
    DOPz=dummy; %Particulate organic P
     
    TIPz=Pz + PPz; % Total inorg. phosphorus (excl. Chla and DOP) in the water column (mg m-3)
    dummy=IOflow_v11(dz, zz, Vz, TIPz, lvlD, Iflw, Iflw_TP-(Iflw_Chl./Y_cp)-Iflw_DOP);
    TIPz=dummy; %Total inorg. phosphorus (excl. Chla and DOP) 
    

    %== P-partitioning in water==
    [Pz, trash]=Ppart(Sz./rho_sed,TIPz,Psat_L,Fmax_L,rho_sed,Fstable);
    PPz=TIPz-Pz;
    
    dummy=IOflow_v11(dz, zz, Vz, Chlz, lvlD, Iflw, Iflw_Chl);
    Chlz=dummy; %Chlorophyll a
    
    dummy=IOflow_v11(dz, zz, Vz, DOCz, lvlD, Iflw, Iflw_DOC);
    DOCz=dummy; %DOC

   else
    lvlD=NaN;
   end %if(Iflw>0)
   
  else
    Iflw=0;
    lvlD=NaN;     
  end  %if (river_inflow_switch==1)

   % Convective mixing adjustment (mix successive layers until stable density profile)  
   % and 
   % Spring/autumn turnover (don't allow temperature jumps over temperature of maximum density)
   [Tz,Cz,Sz,Pz,Chlz,PPz,DOPz,DOCz] = convection_v12_1a(Tz,Cz,Sz,Pz,Chlz,PPz,DOPz,DOCz,Tprof_prev,Vz,Cw,f_par,lambdaz_wtot_avg,zz,swa_b0,tracer_switch,1);
                              
   if (IceIndicator==0)
    
    TKE=C_shelter*Az(1)*sqrt(tau^3/rho(1))*(24*60*60*dt); %Turbulent kinetic energy (J day-1) over the whole lake
     
    %Wind mixing
        WmixIndicator=1;
            Bef_wind=sum(diff(rho)==0); %just a watch variable
        while (WmixIndicator==1)
         d_rho=diff(rho);
         inx=find(d_rho>0);
          if (isempty(inx)==0); %if water column not already fully mixed
           zb=inx(1);   
           MLD=dz*zb; %mixed layer depth
           dD=d_rho(zb); %density difference
           Zg=sum( Az(1:zb+1) .* zz(1:zb+1) ) / sum(Az(1:zb+1)); %Depth of center of mass of mixed layer       
           V_weight=Vz(zb+1)*sum(Vz(1:zb))/(Vz(zb+1)+sum(Vz(1:zb)));
           POE=(dD*g*V_weight*(MLD + dz/2 - Zg));
           KP_ratio=TKE/POE;
            if (KP_ratio>=1)
             Tmix=sum( Vz(1:zb+1).*Tz(1:zb+1) ) / sum(Vz(1:zb+1));
             Tz(1:zb+1)=Tmix;

              if (tracer_switch==1)
               Cmix=sum( Vz(1:zb+1).*Cz(1:zb+1) ) / sum(Vz(1:zb+1));
               Cz(1:zb+1)=Cmix;
              end
              
             Smix=sum( Vz(1:zb+1).*Sz(1:zb+1) ) / sum(Vz(1:zb+1));
             Sz(1:zb+1)=Smix;
              
             Pmix=sum( Vz(1:zb+1).*Pz(1:zb+1) ) / sum(Vz(1:zb+1));
             Pz(1:zb+1)=Pmix;

             Chlmix=sum( Vz(1:zb+1).*Chlz(1:zb+1) ) / sum(Vz(1:zb+1));
             Chlz(1:zb+1)=Chlmix;
             
             PPmix=sum( Vz(1:zb+1).*PPz(1:zb+1) ) / sum(Vz(1:zb+1));
             PPz(1:zb+1)=PPmix;
             
             DOPmix=sum( Vz(1:zb+1).*DOPz(1:zb+1) ) / sum(Vz(1:zb+1));
             DOPz(1:zb+1)=DOPmix;

             DOCmix=sum( Vz(1:zb+1).*DOCz(1:zb+1) ) / sum(Vz(1:zb+1));
             DOCz(1:zb+1)=DOCmix;

             rho = polyval(ies80,max(0,Tz(:))) + min(Tz(:),0);
             TKE=TKE-POE;
            else %if KP_ratio < 1, then mix with the remaining TKE part of the underlying layer 
             Tmix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*Tz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
             Tz(1:zb)=Tmix;
             Tz(zb+1)=KP_ratio*Tmix + (1-KP_ratio)*Tz(zb+1);

              if (tracer_switch==1)
               Cmix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*Cz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
               Cz(1:zb)=Cmix;
               Cz(zb+1)=KP_ratio*Cmix + (1-KP_ratio)*Cz(zb+1);
              end 
              
             Smix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*Sz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
             Sz(1:zb)=Smix;
             Sz(zb+1)=KP_ratio*Smix + (1-KP_ratio)*Sz(zb+1);
                             
             Pmix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*Pz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
             Pz(1:zb)=Pmix;
             Pz(zb+1)=KP_ratio*Pmix + (1-KP_ratio)*Pz(zb+1);

             Chlmix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*Chlz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
             Chlz(1:zb)=Chlmix;
             Chlz(zb+1)=KP_ratio*Chlmix + (1-KP_ratio)*Chlz(zb+1);
             
             PPmix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*PPz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
             PPz(1:zb)=PPmix;
             PPz(zb+1)=KP_ratio*PPmix + (1-KP_ratio)*PPz(zb+1);

             DOPmix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*DOPz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
             DOPz(1:zb)=DOPmix;
             DOPz(zb+1)=KP_ratio*DOPmix + (1-KP_ratio)*DOPz(zb+1);

             DOCmix=sum( [Vz(1:zb); KP_ratio*Vz(zb+1)].*DOCz(1:zb+1) ) / sum([Vz(1:zb); KP_ratio*Vz(zb+1)]);
             DOCz(1:zb)=DOCmix;
             DOCz(zb+1)=KP_ratio*DOCmix + (1-KP_ratio)*DOCz(zb+1);

             rho = polyval(ies80,max(0,Tz(:))) + min(Tz(:),0);
             TKE=0;       
             WmixIndicator=0;
            end %if (KP_ratio>=1)
           else
           WmixIndicator=0;
           end %if water column (not) already mixed               
        end %while

            Aft_wind=sum(diff(rho)==0); %just a watch variable
            
   else % ice cover module               
        XE_surf=(Tz(1)-Tf) * Cw * dz; %Daily heat accumulation into the first water layer (J m-2)
        Tz(1)=Tf; %Ensure that temperature of the first water layer is kept at freezing point
        TKE=0; %No energy for wind mixing under ice
      
      if (Wt(i,3)<Tf) %if air temperature is below freezing          
         %Calculate ice surface temperature (Tice)
         if(WEQs==0) %if no snow
          alfa=1/(10*Hi);
          dHsi=0;
         else
          K_snow=2.22362*(rho_snow/1000)^1.885; %Yen (1981)
          alfa=(K_ice/K_snow)*(((rho_fw/rho_snow)*WEQs)/Hi);
          %Slush/snow ice formation (directly to ice)
          dHsi=max([0, Hi*(rho_ice/rho_fw-1)+WEQs]);         
          Hsi=Hsi+dHsi;
         end
       Tice=(alfa*Tf+Wt(i,3))/(1+alfa);

        %Ice growth by Stefan's law       
       Hi_new=sqrt((Hi+dHsi)^2+(2*K_ice/(rho_ice*L_ice))*(24*60*60)*(Tf-Tice));
        %snow fall
       dWEQnews=0.001*Wt(i,7); %mm->m
       dWEQs=dWEQnews-dHsi*(rho_ice/rho_fw); % new precipitation minus snow-to-snowice in snow water equivalent
       dHsi=0; %reset new snow ice formation       
      else %if air temperature is NOT below freezing
        Tice=Tf;    %ice surface at freezing point
        dWEQnews=0; %No new snow
        if (WEQs>0)
        %snow melting in water equivalents
        dWEQs=-max([0, (60*60*24)*(((1-IceSnowAttCoeff)*Qsw)+Qlw+Qsl)/(rho_fw*L_ice)]);
                if ((WEQs+dWEQs)<0) %if more than all snow melts...
                Hi_new=Hi+(WEQs+dWEQs)*(rho_fw/rho_ice); %...take the excess melting from ice thickness
                else
                Hi_new=Hi; %ice does not melt until snow is melted away
                end
        else    
        %total ice melting
        dWEQs=0;
        Hi_new=Hi-max([0, (60*60*24)*(((1-IceSnowAttCoeff)*Qsw)+Qlw+Qsl)/(rho_ice*L_ice)]);
        %snow ice part melting
        Hsi=Hsi-max([0, (60*60*24)*(((1-IceSnowAttCoeff)*Qsw)+Qlw+Qsl)/(rho_ice*L_ice)]);
            if (Hsi<=0)
                Hsi=0;
            end
        end %if there is snow or not
      end %if air temperature is or isn't below freezing
      
      
      %Update ice and snow thicknesses
      Hi=Hi_new-(XE_surf/(rho_ice*L_ice)); %new ice thickness (minus melting due to heat flux from water)
      XE_surf=0; %reset energy flux from water to ice (J m-2 per day)
      WEQs=WEQs+dWEQs; %new snow water equivalent   
      
      if(Hi<Hsi)
          Hsi=max(0,Hi);    %to ensure that snow ice thickness does not exceed ice thickness 
                            %(if e.g. much ice melting much from bottom)
      end
            
     
      if(WEQs<=0)
       WEQs=0; %excess melt energy already transferred to ice above
       rho_snow=rho_new_snow;
      else
       %Update snow density as weighed average of old and new snow densities
       rho_snow=rho_snow*(WEQs-dWEQnews)/WEQs + rho_new_snow*dWEQnews/WEQs;
         if (snow_compaction_switch==1)
         %snow compaction
           if (Wt(i,3)<Tf) %if air temperature is below freezing        
           rhos=1e-3*rho_snow; %from kg/m3 to g/cm3
           delta_rhos=24*rhos*C1*(0.5*WEQs)*exp(-C2*rhos)*exp(-0.08*(Tf-0.5*(Tice+Wt(i,3))));
           rho_snow=min([rho_snow+1e+3*delta_rhos, max_rho_snow]);  %from g/cm3 back to kg/m3 
           else
           rho_snow=max_rho_snow;    
           end
         end
      end
      
      if(Hi<=0)
      IceIndicator=0;
      disp(['Ice-off, ' datestr(datenum(M_start)+i-1)])
      XE_melt=(-Hi-(WEQs*rho_fw/rho_ice))*rho_ice*L_ice/(24*60*60); 
             %(W m-2) snow part is in case ice has melted from bottom leaving some snow on top (reducing XE_melt)
      Hi=0;
      WEQs=0;
      Tice=NaN;
      DoM(pp)=i;
      pp=pp+1;
      end
      
   end %of ice cover module           
   
  %== P-partitioning in water==
  TIPz=Pz + PPz; % Total inorg. phosphorus (excl. Chla and DOP) in the water column (mg m-3)
  [Pz, trash]=Ppart(Sz./rho_sed,TIPz,Psat_L,Fmax_L,rho_sed,Fstable);
  PPz=TIPz-Pz;
  
        %Initial freezing
            Supercooled=find(Tz<Tf);
            if (isempty(Supercooled)==0)
              IceIndicator=1;
              InitIceEnergy=sum((Tf-Tz(Supercooled)).*Vz(Supercooled)*Cw);
              Hi=Hi+(InitIceEnergy/(rho_ice*L_ice))/Az(1);
              Tz(Supercooled)=Tf;
              Tz(1)=Tf; %set temperature of the first layer to freezing point
              DoF(qq)=i;
              disp(['Ice-on, ' datestr(datenum(M_start)+i-1)])
              qq=qq+1;
            end   
           
        % Calculate pycnocline depth
        pycno_thres=0.1;  %treshold density gradient value (kg m-3 m-1)
        rho = polyval(ies80,max(0,Tz(:))) + min(Tz(:),0);
        dRdz = [NaN; abs(diff(rho))];
        di=find((dRdz<(pycno_thres*dz)) | isnan(dRdz));
        %dRdz(di)=NaN; 
        %TCz = nansum(zz .* dRdz) ./ nansum(dRdz);
        dRdz(di)=0; %modified for MATLAB version 7
        TCz = sum(zz .* dRdz) ./ sum(dRdz);
        
        %vector with S_res_epi above, and S_res_hypo below the pycnocline
        inx=find(zz <= TCz);
        S_resusp(inx)=S_res_epi;
        inx=find(zz > TCz);
        S_resusp(inx)=S_res_hypo;
        
        if (IceIndicator==1)
        S_resusp(:)=S_res_hypo;  %only hypolimnetic type of resuspension allowed under ice
        end
        
        if( isnan(TCz) & (IceIndicator==0) )  
         S_resusp(:)=S_res_epi;   %if no pycnocline and open water, resuspension allowed from top to bottom   
        end
            
            
   % Output matrices     
   Qst(:,i) = [Qsw Qlw Qsl]';
   Kzt(:,i) = [0;Kz];
   Tzt(:,i) = Tz;
   Czt(:,i) = Cz;
   Szt(:,i) = Sz;
   Pzt(:,i) = Pz;
   Chlzt(:,i) = Chlz;
   PPzt(:,i) = PPz;
   DOPzt(:,i) = DOPz;
   DOCzt(:,i) = DOCz;
   Qzt_sed(:,i) = Qz_sed./(60*60*24*dt); %(J m-2 day-1) -> (W m-2)
   lambdazt(:,i) = lambdaz_wtot_avg;

   P3zt_sed(:,i,1) = Pdz_store; %diss. P conc. in sediment pore water (mg m-3)
   P3zt_sed(:,i,2) = Psz_store; %P conc. in inorganic sediment particles (mg kg-1 dry w.)
   P3zt_sed(:,i,3) = Chlsz_store; %Chl conc. in organic sediment particles (mg kg-1 dry w.)
   P3zt_sed(:,i,4) = F_IM; %VOLUME fraction of inorganic particles of total dry sediment
   P3zt_sed(:,i,5) = Sedimentation_counter; %H_netsed_inorg; %Sedimentation (m/day) of inorganic particles of total dry sediment
   P3zt_sed(:,i,6) = Resuspension_counter; %H_netsed_org; %Sedimentation (m/day) of organic particles of total dry sediment
   P3zt_sed(:,i,7) = NewSedFrac; %(monitoring variables)
   
   P3zt_sed_sc(:,i,1) = dPW_up; %(mg m-3 day-1) change in Pz due to exchange with pore water
   P3zt_sed_sc(:,i,2) = dPP; %(mg m-3 day-1)
   P3zt_sed_sc(:,i,3) = dChl_res; %(mg m-3 day-1)
   
   His(1,i) = Hi;
   His(2,i) = (rho_fw/rho_snow)*WEQs;
   His(3,i) = Hsi;
   His(4,i) = Tice;
   His(5,i) = Wt(i,3);
   His(6,i) = rho_snow;
   His(7,i) = IceIndicator; 
   
   MixStat(1,i) = Iflw_S;
   MixStat(2,i) = Iflw_TP;
   MixStat(3,i) = lambdaz_wtot(2);%Iflw_DOC;
        MixStat(4,i) = mean(Growth_bioz);
        MixStat(5,i) = mean(Loss_bioz);
        MixStat(6,i) = Iflw;
          if (IceIndicator == 1)
           MixStat(7:11,i) = NaN;
          else
           dum=interp1(zz,Pz,[0:0.1:4]);
           MixStat(7,i) = mean(dum); %diss-P conc. 0-4m in ice-free period
              
           dum=interp1(zz,Chlz,[0:0.1:4]);
           MixStat(8,i) = mean(dum); %Chla conc. 0-4m in ice-free period
           
           dum=interp1(zz,PPz,[0:0.1:4]);
           MixStat(9,i) = mean(dum); %particulate inorg. P conc. 0-4m in ice-free period

           dum=interp1(zz,DOPz,[0:0.1:4]);
           MixStat(10,i) = mean(dum); %dissolved organic P conc. 0-4m in ice-free period
           
           dum=interp1(zz,Sz,[0:0.1:4]);
           MixStat(11,i) = mean(dum); %particulate matter conc. 0-4m in ice-free period
          end

        MixStat(12,i) = TCz; %pycnocline depth

        MixStat(13,i) = 1e-6*Iflw*Iflw_TP; %total P inflow (kg day-1)
        if (Iflw>Vz(1))
           disp('Large inflow!!')    
        end
        MixStat(14,i) = 1e-6*Iflw*(Pz(1)+PPz(1)+DOPz(1)+Chlz(1)); %total P outflow (kg day-1)
        MixStat(15,i) = sum(1e-6*Vz.*(delPP_inorg + delC_org)); %total P sink due to sedimentation (kg day-1)
        MixStat(16,i) = sum(1e-6*(dPP+dPW_up).*Vz); %Internal P loading (kg day-1, excluding Chla)
        MixStat(17,i) = sum(1e-6*dChl_res.*Vz); %Internal Chla loading (kg day-1)
        MixStat(18,i)= sum(1e-6*Vz.*((Pz+PPz+DOPz+Chlz) - TP0)); %Net P change kg
        MixStat(19,i)= sum(1e-6*((dPP+dPW_up-delPP_inorg+dChl_res-delC_org).*Vz - (1-F_sed_sld)*H_sed*(-diff([Az; 0])).*dPW_down)); %Net P flux from sediment kg
        MixStat(20,i) = 1e-6*Iflw*(Iflw_TP-Iflw_Chl./Y_cp-Iflw_DOP-Fstable*Iflw_S); %total algae-available P inflow (kg day-1)
end; %for i

runtime=toc;

disp(['Total model runtime: ' int2str(floor(runtime/60)) ' min ' int2str(round(mod(runtime,60))) ' s']);
% disp(['Total matter sedimented: ' num2str(round(nansum(Sedimentation_counter))) ' kg']);
% disp(['Total matter resuspended: ' num2str(round(nansum(Resuspension_counter))) ' kg']);
% disp(['Net dry matter sedimented (mean): '  num2str(1000*round(nansum(Sedimentation_counter-Resuspension_counter))./(rho_sed*Az(1))) ' mm']); %m -> mm

% >>>>>> End of the time loop >>>>>>

% Below are the two functions for calculating tridiagonal matrix Fi for solving the 
% 1) diffusion equation (tridiag_DIF_v11), and 
% 2) advection-diffusion equation (tridiag_HAD_v11) by fully implicit hybrid exponential numerical scheme, 
% based on Dhamotharan et al. 1981, 
%'Unsteady one-dimensional settling of suspended sediments', Water Resources Research 17(4), 1125-1132
% code checked by TSA, 16.03.2004


%Inputs:
% Kz    diffusion coefficient at layer interfaces (plus surface) N (N,1)
% U     vertical settling velocity (scalar)
% Vz    layer volumes (N,1)
% Az    layer interface areas (N,1)
% dz    grid size
% dt    time step

%Output:
% Fi    tridiagonal matrix for solving new profile Cz

   az = (dt/dz) * [0; Kz] .* (Az ./ Vz);
   bz = (dt/dz) * [Kz; 0] .* ([Az(2:end); 0] ./ Vz);
   Gi = [-bz (1 + az + bz) -az];

%=== DIFFUSIVE EQUATION ===
function[Fi]=tridiag_DIF_v11(Kz,Vz,Az,dz,dt)

Nz=length(Vz); %number of grid points/layers

% Linearized heat conservation equation matrix (diffusion only)
   az = (dt/dz) * Kz .* (Az ./ Vz);                                        %coefficient for i-1
   cz = (dt/dz) * [Kz(2:end); NaN] .* ([Az(2:end); NaN] ./ Vz);            %coefficient for i+1
   bz = 1 + az + cz;                                                       %coefficient for i+1
%Boundary conditions, surface

  az(1) = 0;
  %cz(1) remains unchanged 
  bz(1)= 1 + az(1) + cz(1);
   

%Boundary conditions, bottom

  %az(end) remains unchanged 
  cz(end) = 0;
  bz(end) = 1 + az(end) + cz(end);
  
   Gi = [-cz bz -az];
   Fi = spdiags(Gi,-1:1,Nz,Nz)';
%end of function


%=== ADVECTIVE-DIFFUSIVE EQUATION ===
function[Fi]=tridiag_HAD_v11(Kz,U,Vz,Az,dz,dt)

if (U<0)
    error('only positive (downward) velocities allowed')
end

if (U==0)
    U=eps; %set Vz next to nothing (=2.2204e-016) in order to avoid division by zero
end

Nz=length(Vz); %number of grid points/layers

theta=U*(dt/dz);

az = theta.*(1 + (1./(exp( (U*Vz)./(Kz.*Az) ) - 1)));                   %coefficient for i-1
cz = theta./(exp( (U*Vz)./([Kz(2:end); NaN].*[Az(2:end); NaN]) ) - 1);  %coefficient for i+1
bz = 1 + az + cz;                                                       %coefficient for i

%Boundary conditions, surface

az(1) = 0;
%cz(1) remains unchanged 
bz(1) = 1 + theta + cz(1);

%Boundary conditions, bottom

%az(end) remains unchanged 
cz(end) = 0;
bz(end) = 1 + az(end);

Gi = [-cz bz -az];
Fi = spdiags(Gi,-1:1,Nz,Nz)';
%end of function


% Below are three functions for calculating the partitioning between
% dissolved and inorganic particle bound phosphorus. 
% based on N.N, 
%'Title', Paper 1(4), 11xx-11xx
% code checked by TSA, x.x.200x
%================================

% function [Pdiss, Pfpart]=Ppart(C,TIP,C_alpha)
% %C:     Suspended inorganic matter (mg/kg)
% %TIP:   Total inorganic phosphorus (mg/m3)
% 
% %opts=optimset('TolX',0.01); %Set P partitioning iteration treshold (for dissolved P conc.)
% %C_alpha=75;     %P isotherm factor
% C_beta=0.33;    %P isotherm exponent
% N=length(TIP);
% 
% if (C_beta==0.33) %If isotherm exponent b=1/3 an analytical solution can be used instead of zero-finding method (speeding up the run)
%     Pdiss=Ppart_analyt(C_alpha*C, TIP);  
% else
%    for w=1:N
%     Pdiss(w) = fzero(@ekviP, [0 10000], [], C, TIP, C_alpha, C_beta); %(mg/m3)
%    end
% end
% Pfpart = (TIP - Pdiss)./C; %inorg. P conc. in sediment particles(mg/kg dry w.) 
% %end of function


function [Pdiss, Pfpart]=Ppart(vf,TIP,Psat,Fmax,rho_sed,Fstable)
% Function for calculating the partitioning between
% dissolved and inorganic particle bound phosphorus. 
% Based on Langmuir isotherm approach 
%vf:    volume fraction of suspended inorganic matter (m3 m-3); S/rho_sed OR (1-porosity)
%TIP:   Total inorganic phosphorus (mg m-3)
%Psat, mg m-3 - Langmuir half-saturation parameter
%Fmax, mg kg-1  - Langmuir scaling parameter
%rho_sed, kg m-3 - Density of dry inorganic sediment mass
%Fstable, mg kg-1 - Inactive P conc. in inorg. particles

N=length(TIP);
Pdiss=NaN*ones(N,1);

for w=1:N
 a = vf(w)-1;
 b = TIP(w) + (vf(w)-1)*Psat - vf(w)*rho_sed*(Fmax+Fstable);
 c = Psat*TIP(w) - vf(w)*rho_sed*Fstable*Psat ;
Pdiss(w) = max(real(roots([a b c])));
end

Pfpart = (TIP - (1-vf).*Pdiss)./(rho_sed*vf); %inorg. P conc. in sediment particles(mg kg-1 dry w.) 
%end of function
       

function [tobezero]=ekviP(Pd,C,Ptot,A,B)
%calculating the phosphorus equilibrium between particles and dissolved
%phase (Webster and Grace, 2001) 
tobezero=C*A*Pd.^B + Pd - Ptot;
%end of function


function [Pdiss]=Ppart_analyt(aa,bb)
%calculating the phosphorus equilibrium between particles and dissolved
%phase (Webster and Grace, 2001) when the isotherm exponent is equal to 1/3
%aa: C_alpha * Suspended inorganic matter (mg/kg)
%bb: Total inorganic phosphorus (mg/m3)
Pdiss = -aa.*(1/6*(108*bb+12*(12*aa.^3+81*bb.^2).^(1/2)).^(1/3)-2*aa./(108.*bb+12.*(12.*aa.^3+81.*bb.^2).^(1/2)).^(1/3))+bb;
%end of function


