function Check_Shim_configuration(varargin)
%-------------------------------------------------------------------------
% Compute the solution of the optimized shimming configuration
%
% Dipole model of PM based on paper by R. Engel-Herbert and T. Hesjedal:
%     R. Engel-Herbert and T. Hesjedal
%     Calculation of the magnetic stray field of a uniaxial magnetic domain
%     J. Appl. Phys., Vol. 97, 2005, 074504
%
% The dipole model also include:
%   - the reaction field due to the other PMs
%   - the nonlinear JH behaviour of PM
%   (B = Mu0*(H+J), with J = f(H))
%
% Author: O. Bottauscio (first version: 2025)
%-------------------------------------------------------------------------
global MAGNETI_STATO
global POINTS
Code_Version='2.3';

fprintf('Code: Check_Shim_configuration, Version: %s\n',Code_Version);
DirRun=pwd;   %Directory di run
fprintf('Running directory: %s\n',DirRun);

warning('off','all');    %ALL Warning message are suppressed
fprintf('ALL Warnings suppressed\n');

% constants
Mu0=4e-7*pi;

% Input from prompt line
n_input=size(varargin,2);
if n_input >= 1
  input_data = varargin{1};
else
  fprintf('Input file missing\n')
end

load(input_data);

MAGNETI_STATO = struct;

[Bmin,Bmax,Bmean,DEV1_ppm,DEV2_ppm] = util.variability(POINTS.B5);
fprintf('Original configuration (before shimming steps) - \n');
fprintf('Bmin,Bmax,Bmean [mT]: %f %f %f\n',Bmin*1000,Bmax*1000,Bmean*1000);
fprintf('DEV1 [ppm]: %f\n',DEV1_ppm);
fprintf('DEV2 [ppm]: %f\n',DEV2_ppm);

B5_tot = deviation_in_sphere3B(Mu0);
[Bmin,Bmax,Bmean,DEV1_ppm,DEV2_ppm] = util.variability(B5_tot);
fprintf('Ending situation - \n');
fprintf('Bmin,Bmax,Bmean [mT]: %.2f %.2f %.2f\n',Bmin*1000,Bmax*1000,Bmean*1000);
fprintf('DEV1 [ppm]: %.1f\n',DEV1_ppm);
fprintf('DEV2 [ppm]: %.1f\n',DEV2_ppm);


end

function B5_tot = deviation_in_sphere3B(Mu0)
global MAGNETI_GEO
global TEMPERATURE
global MAGNETI_STATO
global MATERIALI
global SIMUL_DATA
global POINTS

%-------------------------------------------------------------------
%  Assegnazione parametri magnetici per PM di configurazione shimming
%-------------------------------------------------------------------
MAGNETI_STATO.Temperature_magneti=ones(MAGNETI_GEO.Ndipoli_tot,1)*TEMPERATURE.Tactual;

for n=1:MATERIALI.Nlist
  mat_number=MATERIALI.Lista(n,1);
  Jresidua_relativo(mat_number,1)=1.0;   %For the DET simulations the JH curves are not rescaled
end
[MAGNETI_STATO] = dipoles_sub.assegna_valori_dipoles(MAGNETI_GEO,MATERIALI,MAGNETI_STATO,Jresidua_relativo);
recompute=true;
first=true;
% Initialization
L=0;
U=0;
P=0;
Tnoto2=0;
output_JHmag=true;
[~,~,~,~,BFIELD,MAGNETI_STATO] = dipoles_sub.output_dipoles(first,SIMUL_DATA.reaction,SIMUL_DATA.NL,...
    recompute,MAGNETI_GEO,MAGNETI_STATO,MATERIALI,POINTS,L,U,P,Tnoto2,Mu0,output_JHmag);
B5 = util.estrai_output(SIMUL_DATA,BFIELD);

B5_tot=POINTS.B5+B5;
return
end