%% mg/dL (mg%) --> M
val = 150;

mw = 46.068; % molar mass of 
mg2g = 1e-3;
dL2L = 10;

val*mg2g/mw*dL2L

%% ug/uL --> M
val = 100;

mw = 46.068; % molar mass of 
ug2g = 1e-6;
uL2L = 1e6;

val*ug2g/mw*uL2L

%% v/v% --> M

val = 1.8;
density = 0.79; %g/mL
ml2L = 1000;

val / 100 * density / mw * ml2L

%% v/v% --> mg/dL

val = 0.018
mL2dL = 100;
mg2g =1000;

val*mL2dL*density*mg2g