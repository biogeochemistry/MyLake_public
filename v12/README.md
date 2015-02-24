MyLake see MyLake_v12_1_note.pdf
======

Revisions made to MyLake v.1.2.1 
21.08.07 by Tuomo Saloranta 
The basis for the revisions implemented in code v.1.2.1 and described here were mainly to improve 
the robustness and performance of the MyLake v.1.2 code. In addition, the photodegradation of DOC 
as well as two phytoplankton groups were added to the code. Some of the revisions were based on 
issues that arose during applications and testing ofMyLake v.1.2. These revisions are implemented in 
two new code files (convection_v12_1a.m, solvemodel_v12_1b.m) and a new parameter file. 
NB! In order to implement the new parameter file correctly, the vector index “40” must be changed to 
“48” in code file “modelinputs_v12.m” (rows 53-55). In addition, to avoid annoying "division by 
zero" warnings under polar night conditions, onecan add the Dayfrac==0 condition on the ifstatement on row 76 in heatflux_v12.m, like this: 
if((isnan(GR))|(Dayfrac==0)) % if global radiation input is missing, or Dayfrac=0, calculate total 
transmissivity from formulaes 

1. Problems encountered with v.1.2 
Three model performance issues were noted during testsimulation of some 350 lakes in a batch with 
climatological mean monthly forcing and assumedly no river input (CRU/Euregi simulations; see the 
corresponding research note), as well as with simulations in Lake Vansjø: 
1)  The lake ice cover on/off date sometimes “flips” and produced e.g. up to 40 on/off dates 
during one year for a lake with the minimum mean monthly air temperature slightly above 
zero. As there was no threshold limit for ice cover, these “flips” usually occurred in the 
beginning of the freezing period. This issue was resolved by introducing a frazil ice process as 
explained below (section 2.1). 
2)  There were some negative values in chlorophyll and dissolved P for some lakes. In the code 
there was a restriction that phytoplankton cannot grow more than there is dissolved P left. 
However, there was a convection process between the calculation of the growth rates and their 
application which mixed the profiles sometimes so that negative values were generated. These 
negative values caused also some numerical instability. This issue was resolved by rearranging 
some lines of the code as explained below (section 2.2). 
3)  The new spring/autumn turnover algorithm in the revised convection file (used since june 
2006, see the development log) produced some abnormal peaks in temperature. These were 
probably related to complications caused by the river inflow. This issue was resolved by 
taking the older version of the convection file into use, and modifying it slightly (section 2.2).
2. Solutions 
2.1 Frazil ice 
In the initial phase of formation of ice cover the ice crystals are suspended in the water column (socalled “frazil ice”). After the mass of this frazil ice grows the ice crystals eventually float to the 
surface and form a slushy layer which freezes to form the initial ice cover. The more wind and 
turbulence in the water column there is the more frazil ice the water column can sustain before it forms 
the initial ice cover. In the MyLake v.1.2.1 the frazil ice is handled as explained in the following. 
When supercooled water is encountered this is turnedinto ice (by turning the sensible heat deficit in 
the supercooled water into latent heat of freezing). Inthe previous version 1.2 of the model code this 
formed the initial ice cover. However, now this ice is defined as frazil ice and allowed to accumulated 
up to a threshold thickness (now set to 3 cm, but this could later be defined e.g.as function of the wind 
mixing energy) before it is turned into the initial solid ice cover in the model. It is first after the 
formation of this solid ice cover that the ice related processes are enabled in the model (new 
calculations for heat fluxes, no wind mixing, new albedo, etc.). Thus the model functions in “open 
water mode” under frazil ice conditions. The “thickness” of frazil ice (i.e. the volume of frazil ice 
divided by the lake surface area) increases whenever new supercooled water is encountered, and 
decreases whenever the water columnreceives heat to melt the frazil ice. Since no mixing depth of the 
frazil ice is modelled, the possible excess heat (in case all the frazil ice is melted) is distributed 
similarly as the original heating profile. 
2.2 Code rearrangements 
In order to solve the problem encountered with the spring/autumn turnover code in the revised 
convection file, this file was replaced by the original convection code, where only the surface 
temperature “jumps” over the temperature of the maximum density are controlled (instead of all such 
jumps in the water column, as in the revised version). Three lines are modified in the original code 
(see convection_v12_1a.m). Also in solvemodel_v12_1b.m, the spring/autumn turnover is only 
allowed after the two first occasions of convection_v12_1a.m, after surface heat flux calculations. The 
previous temperature profile variable is updated after each occasion of convection_v12_1a.m. The 
new modifications to the original convection filedid not change significantly the 1990-1999 TotP, 
Chl, or suspended solids mean values in Lake Vansjø-Storefjorden (old: 17.467; 6.502; 0.00344; new: 
17.468; 6.502; 0.00344). If purely non-surface jumps occur, the model code now displays a note about 
this on the command line. 
In order to solve the problem encountered with the negative chemistry values, the calculation of the 
growth rate in solvemodel_v12_1b.m was moved to after the solving of the temperature profile and the 
consequent convection routine. The PAR light profile is here recalculated too. These new code 
modifications did not change significantly the 1990-1999 TotP, Chl, or suspended solids mean values 
in Lake Vansjø-Storefjorden (old: 17.467; 6.502; 0.00344; new: 17.382; 6.474; 0.00343). 
3. Other new features 
3.1 Photodegradation of DOC 
The DOC in the water column was in the previous version 1.2 of the model code modelled as a passive 
dissolved tracer. In the MyLake v.1.2.1 a photodegradation process due to solar radiation is added and 
the rate of degradation ?CDOC [mg m
-3
day
-1
] is formulated in the following way: 
PAR PAR DOC q DOC DOC
I E y C C · · · · = ? ) / 1 ( ß
where CDOCis the DOC concentration [mg m
-3
] (taken from the previous time step in the model), yqis 
the degradation quantum yield [mg DOC / mol quanta], ßDOC is the optical cross-section of DOC [m
2
mg
-1
DOC], EPARis the energy of photons [J/mol quanta], and IPARis the PAR irradiation [J m
-2
day
-1
]. 
Thus two new parameters yqand ßDOC are introduced in the model code. 
3.2 Two phytoplankton groups 
A new phytoplankton group, in terms of a new chlorophyll (PChl) variable, is introduced in the model 
code. In practise this means that a second set of values is given for six phytoplankton-related 
parameters: the settling velocity wChl,light saturation level of photosynthesis I’, specific loss rate at 20 
°Cm(20), the maximal attainable specific growth rate at 20 °C µ’(20), the half-saturation parameter P’
and the optical cross section of chlorophyll ßC. Now two different phytoplankton growth and loss rates 
are calculated, but otherwise the solving of the equations remain the same as before in the model code. 
The resuspension of chlorophyll from sediments is assumed to be divided equally (50-50) between the 
two phytoplankton groups. The variable for dissolve tracer (Cz) in the model code is now adopted for 
the second phytoplankton (chlorophyll) group. 
3.3 Miscellaneous 
The model parameter file was extended to cover the eight new DOC and second phytoplankton group 
parameters. NB! In order to implement this extended filecorrectly, the vector index “40” must be 
changed to “48” in code file “modelinputs_v12.m” (rows 53-55). 
As the old variable for dissolve tracer (Cz) is now used for the second phytoplankton (chlorophyll) 
group, the “tracer_switch” variable has become redundantin the code and must be set equal to 1. (This 
variable is retained in the convection code for version compatibility). 
A small bug was corrected in the calculation of the attenuated solar radiation (variable “Attn_z”). This 
bug caused the sum of “Attn_z” to be slightly above one (should be equal to one). 
(12/11/2007) 
A bug, which caused a model error with latitudes northof the polar circle (ca. 66 N), was fixed. 