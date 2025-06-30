function Array_Dipoles_Shim(varargin)
%-------------------------------------------------------------------------
% Optimization of the shimming for an array of permanent magnets (PMs)
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
global MAGNETI_GEO
global SHIMMING_OPT
global SHIMMING_GEO
global MAGNETI_STATO
global MATERIALI
global TEMPERATURE
global SIMUL_DATA
global POINTS
Code_Version='1.2';

fprintf('Code: Array_Dipoles_Shim, Version: %s\n',Code_Version);
DirRun=pwd;   %Directory di run
fprintf('Running directory: %s\n',DirRun);

% constants
Mu0=4e-7*pi;

% Input from prompt line
n_input=size(varargin,2);
if n_input >= 1
  input_data = varargin{1};
else
  fprintf('Input file missing\n')
end
if ~contains(input_data,'.toml')
  input_data=append(input_data,'.toml');
end

% Opening input data file
if isfile(input_data) == 0
  fprintf('File: %s does not exist\n',input_data)
  return
end

% Reading toml file
toml_data = toml.read(input_data);

input_file_description = toml_data.title;
input_file_version = toml_data.version;

% Units
[scale] = read_labels.readUnit(toml_data);

% Read shimming geometrica data
SHIMMING_GEO = struct;
[ierr,SHIMMING_GEO] = input_shimming_geo(toml_data,SHIMMING_GEO,scale);
if ierr > 0
    return
end

MAGNETI_GEO = struct;

MAGNETI_STATO = struct;

% Read material data
MATERIALI = struct;
[ierr,MATERIALI] = read_labels.readMATERIAL(toml_data,MATERIALI);
if ierr > 0
    return
end

% Read magnet temperature status
TEMPERATURE = struct;
[ierr,TEMPERATURE] = read_labels.readTEMPERATURE(toml_data,TEMPERATURE);
if ierr > 0
    return
end

% Read measurements in DSV
POINTS = struct;
[ierr,POINTS] = input_DSVMeas(toml_data,POINTS);
if ierr > 0
    return
end

% Input simulations data
SIMUL_DATA = struct;
[ierr,SIMUL_DATA] = read_labels.readSIMUL(toml_data,SIMUL_DATA);
if ierr > 0
    return
end

if strcmpi(MATERIALI.JHcurve,'NO') && SIMUL_DATA.NL
  fprintf('Nonlinear simulation impossible without JH curves\n')
  return
end

% Read data for optimization
SHIMMING_OPT = struct;
[ierr,SHIMMING_OPT,SHIMMING_GEO,MAGNETI_GEO_PREVIOUS] = input_shimming_opt(toml_data,SHIMMING_OPT,SHIMMING_GEO,POINTS,scale);
if ierr==1
    return
end

SHIMMING_OPT.Mu0=Mu0;

OPT=true;
if ~OPT
  x_input(1,1:SHIMMING_OPT.Ndipoli_tot_shim_max)=1;
  x_input(1,SHIMMING_OPT.Ndipoli_tot_shim_max+1:2*SHIMMING_OPT.Ndipoli_tot_shim_max)=0;
  deviation_tot_ppm = deviation_in_sphere3(x_input);
else
  [Bmin,Bmax,Bmean,DEV1_ppm,DEV2_ppm] = util.variability(SHIMMING_OPT.B5);
  fprintf('Starting situation - \n');
  fprintf('Bmin,Bmax,Bmean [mT]: %f %f %f\n',Bmin*1000,Bmax*1000,Bmean*1000);
  fprintf('DEV1 [ppm]: %f\n',DEV1_ppm);
  fprintf('DEV2 [ppm]: %f\n',DEV2_ppm);

  if strcmpi(SHIMMING_OPT.DevOpt,'DEV1')
    SHIMMING_OPT.original_deviation_ppm=DEV1_ppm;
  elseif strcmpi(SHIMMING_OPT.DevOpt,'DEV2')
    SHIMMING_OPT.original_deviation_ppm=DEV2_ppm;
  end

  nvars=3*SHIMMING_OPT.Ndipoli_tot_shim_max;
  MaxGenerations1=SHIMMING_OPT.MaxGenerations1;
  MaxStallGenerations=SHIMMING_OPT.MaxStallGenerations;
  Tolerance=SHIMMING_OPT.Tolerance;
  
  fprintf('----OPTIMIZATION----\n');
  fprintf('N. of parameters: %d\n',nvars);
  MaxGenerations=MaxGenerations1*nvars;
  options = optimoptions('ga','MaxGenerations',MaxGenerations,...
                              'MaxStallGenerations',MaxStallGenerations,...
                              'FunctionTolerance',Tolerance);
  nnn=nvars/3; %tre blocchi di variabili
% Block 1: PM yes/no
% Block 2: PM position in sector (1:max_PM_per_sector)
% Block 3: PM rotation (-Variazione_angolo_PM_shim:Variazione_angolo_PM_shim)

  intcon=1:nnn*2;   %I primi 2/3 degli input sono interi
  x_input_Min(1,1:nnn)=0;
  x_input_Max(1,1:nnn)=1;
  x_input_Min(1,nnn+1:2*nnn)=1;
  x_input_Max(1,nnn+1:2*nnn)=SHIMMING_GEO.Max_PM_per_settore;
  x_input_Min(1,2*nnn+1:nvars)=-SHIMMING_OPT.Variazione_angolo_PM_shim;
  x_input_Max(1,2*nnn+1:nvars)=SHIMMING_OPT.Variazione_angolo_PM_shim;
  f = @(x_input) deviation_in_sphere3(x_input);
%  
  [x_input,evaluation_opt,exitflag] = ga(f,nvars,[],[],[],[],x_input_Min,x_input_Max,[],intcon,options);

  fprintf('Optimal value of Deviation after shimming (ppm): %f\n',evaluation_opt);
  fprintf('Original value of Deviation before shimming (ppm): %f\n',SHIMMING_OPT.original_deviation_ppm);
  fprintf('Exit flag: %d\n',exitflag);

  mag_yes_no=transpose(x_input(1,1:nnn));
  pos_in_sector=transpose(x_input(1,nnn+1:2*nnn));
  dir_mag_shim=transpose(x_input(1,2*nnn+1:nvars));

% Reconstruct position of PM added in shimming and adjust final field in DSV
  configuration_at_end(mag_yes_no,pos_in_sector,dir_mag_shim);

  if SHIMMING_OPT.Run == 1
      mkdir(SHIMMING_OPT.ShimmingDir);
  end
  filemat=append(SHIMMING_OPT.ShimmingDir,'\',SHIMMING_OPT.ending_state);
  if MAGNETI_GEO_PREVIOUS.exist
     MAGNETI_GEO = assemblea_magneti(MAGNETI_GEO_PREVIOUS,MAGNETI_GEO);
     Delta_magnets=MAGNETI_GEO.Ndipoli_tot-MAGNETI_GEO_PREVIOUS.Ndipoli_tot;
  else
     Delta_magnets=MAGNETI_GEO.Ndipoli_tot;
  end
  fprintf('Added magnets: %d\n',Delta_magnets);
  MatrixShimming=SHIMMING_GEO.MatrixShimming;
  FieldValues=SHIMMING_OPT.B5;
  Run=SHIMMING_OPT.Run;
  Deviation_after_optimization=evaluation_opt;
  save(filemat,'MAGNETI_GEO','MatrixShimming','FieldValues','Run','Delta_magnets','Deviation_after_optimization');
  save(filemat,'Code_Version','input_file_description','input_file_version','-append');
%
  [Bmin,Bmax,Bmean,DEV1_ppm,DEV2_ppm] = util.variability(SHIMMING_OPT.B5);
  fprintf('Ending situation - \n');
  fprintf('Bmin,Bmax,Bmean [mT]: %f %f %f\n',Bmin*1000,Bmax*1000,Bmean*1000);
  fprintf('DEV1 [ppm]: %f\n',DEV1_ppm);
  fprintf('DEV2 [ppm]: %f\n',DEV2_ppm);
  Bmean_mT=Bmean*1000;
  save(filemat,'DEV1_ppm','DEV2_ppm','Bmean_mT','-append');

  fprintf('Results of shimming saved in file: %s\n',filemat);
end


end


function MAGNETI_GEO = assemblea_magneti(MAGNETI_GEO_PREVIOUS,MAGNETI_GEO)
nold = MAGNETI_GEO_PREVIOUS.Ndipoli_tot;
new = MAGNETI_GEO.Ndipoli_tot;
MAGNETI_GEO_PREVIOUS.xdip(nold+1:nold+new,1)=MAGNETI_GEO.xdip(1:new,1);
MAGNETI_GEO_PREVIOUS.ydip(nold+1:nold+new,1)=MAGNETI_GEO.ydip(1:new,1);
MAGNETI_GEO_PREVIOUS.zdip(nold+1:nold+new,1)=MAGNETI_GEO.zdip(1:new,1);
MAGNETI_GEO_PREVIOUS.angle(nold+1:nold+new,1)=MAGNETI_GEO.angle(1:new,1);
MAGNETI_GEO_PREVIOUS.cubo_dimU2(nold+1:nold+new,1)=MAGNETI_GEO.cubo_dimU2(1:new,1);
MAGNETI_GEO_PREVIOUS.cubo_dimV2(nold+1:nold+new,1)=MAGNETI_GEO.cubo_dimV2(1:new,1);
MAGNETI_GEO_PREVIOUS.cubo_dimW2(nold+1:nold+new,1)=MAGNETI_GEO.cubo_dimW2(1:new,1);
MAGNETI_GEO_PREVIOUS.mat_number(nold+1:nold+new,1)=MAGNETI_GEO.mat_number(1:new,1);

MAGNETI_GEO=MAGNETI_GEO_PREVIOUS;
MAGNETI_GEO.Ndipoli_tot=nold+new;
MAGNETI_GEO=rmfield(MAGNETI_GEO,'exist');
return
end



function configuration_at_end(mag_yes_no,pos_in_sector,dir_mag_shim)
global MAGNETI_GEO
global SHIMMING_OPT
global MAGNETI_STATO
global MATERIALI
global TEMPERATURE
global SIMUL_DATA
global POINTS

build_dipole_shim_end(mag_yes_no,pos_in_sector,dir_mag_shim);
%
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
    recompute,MAGNETI_GEO,MAGNETI_STATO,MATERIALI,POINTS,L,U,P,Tnoto2,SHIMMING_OPT.Mu0,output_JHmag);
B5 = util.estrai_output(SIMUL_DATA,BFIELD);
SHIMMING_OPT.B5=SHIMMING_OPT.B5+B5;
return
end


function [ierr,POINTS] = input_DSVMeas(toml_data,POINTS)
ierr=0;
if ~isfield(toml_data.shimming,'measurement')
    fprintf('Field [shimming.pm.measurement] not present\n');
    ierr=1;
    return
end
if ~isfield(toml_data.shimming.measurement,'filemis')
    fprintf('Field [shimming.pm.measurement.filemis] not present\n');
    ierr=1;
    return
end
if ~isfield(toml_data.shimming.measurement,'component')
    fprintf('Field [shimming.pm.measurement.component] not present\n');
    ierr=1;
    return
end
nomemis=toml_data.shimming.measurement.filemis;
aus=load(nomemis);
POINTS.Npoint=size(aus.acqPos,1);
POINTS.xp=aus.acqPos(:,1)*1e-3;
POINTS.yp=aus.acqPos(:,2)*1e-3;
POINTS.zp=aus.acqPos(:,3)*1e-3;
if strcmpi(toml_data.shimming.measurement.component,'x')
  POINTS.B5=abs(aus.realField(:,1))/1000.;
elseif strcmpi(toml_data.shimming.measurement.component,'y')
  POINTS.B5=abs(aus.realField(:,2))/1000.;
elseif strcmpi(toml_data.shimming.measurement.component,'z')
  POINTS.B5=abs(aus.realField(:,3))/1000.;
elseif strcmpi(toml_data.shimming.measurement.component,'max')
  POINTS.B5=sqrt(aus.realField(:,1)^2+aus.realField(:,2)^2+aus.realField(:,3)^2)/1000.;
end
return
end


function [ierr,SHIMMING_OPT,SHIMMING_GEO,MAGNETI_GEO_PREVIOUS] = input_shimming_opt(toml_data,SHIMMING_OPT,SHIMMING_GEO,POINTS,scale)
ierr=0;
if ~isfield(toml_data.shimming,'opt_rules')
    fprintf('Field [shimming.opt_rules] not present\n');
    ierr=1;
    return
end
if ~isfield(toml_data.shimming.opt_rules,'PM_add_per_sector')
    fprintf('Field [shimming.opt_rules.PM_add_per_sector] not present\n');
    ierr=1;
    return
end
if ~isfield(toml_data.shimming.opt_rules,'PM_angleMax')
    fprintf('Field [shimming.opt_rules.PM_angleMax] not present\n');
    ierr=1;
    return
end
if ~isfield(toml_data.shimming.opt_rules,'PM_size')
    fprintf('Field [shimming.opt_rules.PM_size] not present\n');
    ierr=1;
    return
end
if ~isfield(toml_data.shimming.opt_rules,'PM_material')
    fprintf('Field [shimming.opt_rules.PM_material] not present\n');
    ierr=1;
    return
end
if ~isfield(toml_data.shimming.opt_rules,'Zaxis_symmetry')
    fprintf('Field [shimming.opt_rules.Zaxis_symmetry] not present\n');
    ierr=1;
    return
end
if ~isfield(toml_data.shimming.opt_rules,'ring_to_be_used')
    fprintf('Field [shimming.opt_rules.ring_to_be_used] not present\n');
    ierr=1;
    return
end
if ~isfield(toml_data.shimming.opt_rules,'ring_to_be_used')
    fprintf('Field [shimming.opt_rules.ring_to_be_used] not present\n');
    ierr=1;
    return
end

SHIMMING_OPT.PM_add_per_sector=toml_data.shimming.opt_rules.PM_add_per_sector;
SHIMMING_OPT.Variazione_angolo_PM_shim=toml_data.shimming.opt_rules.PM_angleMax;
SHIMMING_OPT.cubo_dimU2_shim_ref=toml_data.shimming.opt_rules.PM_size(1)*scale;
SHIMMING_OPT.cubo_dimV2_shim_ref=toml_data.shimming.opt_rules.PM_size(2)*scale;
SHIMMING_OPT.cubo_dimW2_shim_ref=toml_data.shimming.opt_rules.PM_size(3)*scale;
SHIMMING_OPT.PM_material=toml_data.shimming.opt_rules.PM_material;
SHIMMING_OPT.Zaxis_symmetry=toml_data.shimming.opt_rules.Zaxis_symmetry;

SHIMMING_OPT.ring_to_be_used=toml_data.shimming.opt_rules.ring_to_be_used;
% Verify consistency with used rings
if ~all(SHIMMING_GEO.ring_available(SHIMMING_OPT.ring_to_be_used))
    ierr=1;
    return
end
if SHIMMING_OPT.Zaxis_symmetry
    aaa=SHIMMING_OPT.ring_to_be_used+SHIMMING_GEO.Nring_tot/2;
    if ~all(SHIMMING_GEO.ring_available(aaa))
        ierr=1;
        return
    end
end
Nrused=length(SHIMMING_OPT.ring_to_be_used);
SHIMMING_OPT.Ndipoli_tot_shim_max=SHIMMING_OPT.PM_add_per_sector*SHIMMING_GEO.Nsector_Shim*Nrused;

if isfield(toml_data.shimming.opt_rules,'dev_to_be_optimized')
  SHIMMING_OPT.DevOpt=toml_data.shimming.opt_rules.dev_to_be_optimized;
else
  SHIMMING_OPT.DevOpt='DEV1';
end

if isfield(toml_data.shimming.opt_rules,'shimming_Dir')
  SHIMMING_OPT.ShimmingDir=toml_data.shimming.opt_rules.shimming_Dir;
else
  SHIMMING_OPT.ShimmingDir='ShimmingDir';
end

if isfield(toml_data.shimming.opt_rules,'ending_state')
  SHIMMING_OPT.ending_state=toml_data.shimming.opt_rules.ending_state;
else
  SHIMMING_OPT.ending_state='shimming_output';
end

if isfield(toml_data.shimming.opt_rules,'restart')
  SHIMMING_OPT.restart=toml_data.shimming.opt_rules.restart;
else
  SHIMMING_OPT.restart=false;
end

if SHIMMING_OPT.restart
  if ~isfield(toml_data.shimming.opt_rules,'starting_state')
    fprintf('Field [shimming.opt_rules.starting_state] not present\n');
    ierr=1;
    return
  end
  SHIMMING_OPT.starting_state=toml_data.shimming.opt_rules.starting_state;
  filemat=append(SHIMMING_OPT.ShimmingDir,'\',SHIMMING_OPT.starting_state);
  aus=load(filemat);
  SHIMMING_GEO.MatrixShimming=aus.MatrixShimming;
  MAGNETI_GEO_PREVIOUS=aus.MAGNETI_GEO;
  MAGNETI_GEO_PREVIOUS.exist=true;
  SHIMMING_OPT.B5=aus.FieldValues;
  SHIMMING_OPT.Run=aus.Run+1;
  clear aus
else
  MAGNETI_GEO_PREVIOUS.exist=false;
  SHIMMING_OPT.B5=POINTS.B5;
  SHIMMING_OPT.Run=1;
end

if isfield(toml_data.shimming.opt_rules,'MaxGenerations1')
  SHIMMING_OPT.MaxGenerations1=toml_data.shimming.opt_rules.MaxGenerations1;
else
  SHIMMING_OPT.MaxGenerations1=100;
end
if isfield(toml_data.shimming.opt_rules,'MaxStallGenerations')
  SHIMMING_OPT.MaxStallGenerations=toml_data.shimming.opt_rules.MaxStallGenerations;
else
  SHIMMING_OPT.MaxStallGenerations=50;
end
if isfield(toml_data.shimming.opt_rules,'Tolerance')
  SHIMMING_OPT.Tolerance=toml_data.shimming.opt_rules.Tolerance;
else
  SHIMMING_OPT.Tolerance=1e-4;
end
return
end



function [ierr,SHIMMING_GEO] = input_shimming_geo(toml_data,SHIMMING_GEO,scale)
ierr=0;
if ~isfield(toml_data,'shimming')
    fprintf('Field [shimming] not present\n');
    ierr=1;
    return
end
if ~isfield(toml_data.shimming,'pm')
    fprintf('Field [shimming.pm] not present\n');
    ierr=1;
    return
end
if ~isfield(toml_data.shimming.pm,'raggio')
    fprintf('Field [shimming.pm.raggio] not present\n');
    ierr=1;
    return
end
if ~isfield(toml_data.shimming.pm,'Nsector')
    fprintf('Field [shimming.pm.Nsector] not present\n');
    ierr=1;
    return
end
if ~isfield(toml_data.shimming.pm,'sector_angular_width')
    fprintf('Field [shimming.pm.sector_angular_width] not present\n');
    ierr=1;
    return
end
if ~isfield(toml_data.shimming.pm,'max_PM_per_sector')
    fprintf('Field [shimming.pm.max_PM_per_sector] not present\n');
    ierr=1;
    return
end
if ~isfield(toml_data.shimming.pm,'ring_available')
    fprintf('Field [shimming.pm.ring_available] not present\n');
    ierr=1;
    return
end
if ~isfield(toml_data.shimming.pm,'ring_Zposition')
    fprintf('Field [shimming.pm.ring_Zposition] not present\n');
    ierr=1;
    return
end
SHIMMING_GEO.Raggio_shim=toml_data.shimming.pm.raggio;
SHIMMING_GEO.Nsector_Shim=toml_data.shimming.pm.Nsector;
SHIMMING_GEO.Ampiezza_angolare_settore=toml_data.shimming.pm.sector_angular_width;
SHIMMING_GEO.Max_PM_per_settore=toml_data.shimming.pm.max_PM_per_sector;
SHIMMING_GEO.ring_available=cell2mat(toml_data.shimming.pm.ring_available);
SHIMMING_GEO.ring_Zposition=toml_data.shimming.pm.ring_Zposition*scale;
SHIMMING_GEO.Nring_tot=length(SHIMMING_GEO.ring_Zposition);

Dalfa_shim=SHIMMING_GEO.Ampiezza_angolare_settore/SHIMMING_GEO.Max_PM_per_settore;   %Angolo in gradi di separazione tra magneti del settore
Dang_shim=360./SHIMMING_GEO.Nsector_Shim;
aus=0:1:SHIMMING_GEO.Nsector_Shim-1;
Pos_Angular_shim=Dang_shim*aus;
% Structure which keep trace of position filled with magnets (the value corresponds to the Run number)
SHIMMING_GEO.MatrixShimming=zeros(SHIMMING_GEO.Nring_tot,SHIMMING_GEO.Nsector_Shim,SHIMMING_GEO.Max_PM_per_settore);  %Matrice che conterrà posizioni occupate

SHIMMING_GEO.angoloslot=zeros(SHIMMING_GEO.Nsector_Shim,SHIMMING_GEO.Max_PM_per_settore);  %Angular position of each slot of the sectors
for Ns=1:SHIMMING_GEO.Nsector_Shim
  for Nm=1:SHIMMING_GEO.Max_PM_per_settore
    mm=mod(SHIMMING_GEO.Max_PM_per_settore,2);
    if mm==0
      Nsub=SHIMMING_GEO.Max_PM_per_settore/2+0.5;
      SHIMMING_GEO.angoloslot(Ns,Nm)=Pos_Angular_shim(Ns)+(Nm-Nsub)*Dalfa_shim;
    else
      Nsub=fix(SHIMMING_GEO.Max_PM_per_settore/2)+1;
      SHIMMING_GEO.angoloslot(Ns,Nm)=Pos_Angular_shim(Ns)+(Nm-Nsub)*Dalfa_shim;
    end
  end
end
return
end


function build_dipole_shim(mag_yes_no,pos_in_sector,dir_mag_shim)

global MAGNETI_GEO
global SHIMMING_OPT
global SHIMMING_GEO
%----------------------------------------------------------------------
% Input dati magneti shimming
%----------------------------------------------------------------------
Nring_used=length(SHIMMING_OPT.ring_to_be_used);   %Number of used rings in this opt run
Nmax=SHIMMING_OPT.Ndipoli_tot_shim_max;
if SHIMMING_OPT.Zaxis_symmetry
  Nmax=Nmax*2;
end
MAGNETI_GEO.xdip=zeros(Nmax,1);
MAGNETI_GEO.ydip=zeros(Nmax,1);
MAGNETI_GEO.zdip=zeros(Nmax,1);
MAGNETI_GEO.angle=zeros(Nmax,1);
MAGNETI_GEO.cubo_dimU2=zeros(Nmax,1);
MAGNETI_GEO.cubo_dimV2=zeros(Nmax,1);
MAGNETI_GEO.cubo_dimW2=zeros(Nmax,1);
MAGNETI_GEO.mat_number=zeros(Nmax,1);
MAGNETI_GEO.Ndipoli_tot=0;
ii=0;
for Nu=1:Nring_used
  Nr=SHIMMING_OPT.ring_to_be_used(Nu);   %Number of used Ring
  for Ns=1:SHIMMING_GEO.Nsector_Shim
    for Nm=1:SHIMMING_OPT.PM_add_per_sector
      ii=ii+1;
      if ii > Nmax
          fprintf('ERRORE: conteggio ii\n')
          return
      end
      if mag_yes_no(ii,1)==1
        ipos=pos_in_sector(ii,1);
        angolo=SHIMMING_GEO.angoloslot(Ns,ipos);
        if SHIMMING_GEO.MatrixShimming(Nr,Ns,ipos) == 0
          MAGNETI_GEO.Ndipoli_tot=MAGNETI_GEO.Ndipoli_tot+1;
          MAGNETI_GEO.xdip(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_GEO.Raggio_shim*cos(angolo/180*pi);
          MAGNETI_GEO.ydip(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_GEO.Raggio_shim*sin(angolo/180*pi);
          MAGNETI_GEO.zdip(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_GEO.ring_Zposition(Nr);
          MAGNETI_GEO.angle(MAGNETI_GEO.Ndipoli_tot,1)=2*angolo+dir_mag_shim(ii,1);
          MAGNETI_GEO.cubo_dimU2(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.cubo_dimU2_shim_ref;
          MAGNETI_GEO.cubo_dimV2(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.cubo_dimV2_shim_ref;
          MAGNETI_GEO.cubo_dimW2(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.cubo_dimW2_shim_ref;
          MAGNETI_GEO.mat_number(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.PM_material;
        end
        if SHIMMING_OPT.Zaxis_symmetry
          Nr2=Nr+SHIMMING_GEO.Nring_tot/2;
          if SHIMMING_GEO.MatrixShimming(Nr2,Ns,ipos) == 0
            MAGNETI_GEO.Ndipoli_tot=MAGNETI_GEO.Ndipoli_tot+1;
            MAGNETI_GEO.xdip(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_GEO.Raggio_shim*cos(angolo/180*pi);
            MAGNETI_GEO.ydip(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_GEO.Raggio_shim*sin(angolo/180*pi);
            MAGNETI_GEO.zdip(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_GEO.ring_Zposition(Nr2);
            MAGNETI_GEO.angle(MAGNETI_GEO.Ndipoli_tot,1)=2*angolo+dir_mag_shim(ii,1);
            MAGNETI_GEO.cubo_dimU2(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.cubo_dimU2_shim_ref;
            MAGNETI_GEO.cubo_dimV2(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.cubo_dimV2_shim_ref;
            MAGNETI_GEO.cubo_dimW2(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.cubo_dimW2_shim_ref;
            MAGNETI_GEO.mat_number(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.PM_material;
          end
        end
      end
    end
  end
end
return
end


function build_dipole_shim_end(mag_yes_no,pos_in_sector,dir_mag_shim)

global MAGNETI_GEO
global SHIMMING_OPT
global SHIMMING_GEO

MAGNETI_GEO.Ndipoli_tot=0;
MAGNETI_GEO.xdip(1:end)=[];
MAGNETI_GEO.ydip(1:end)=[];
MAGNETI_GEO.zdip(1:end)=[];
MAGNETI_GEO.angle(1:end)=[];
MAGNETI_GEO.cubo_dimU2(1:end)=[];
MAGNETI_GEO.cubo_dimV2(1:end)=[];
MAGNETI_GEO.cubo_dimW2(1:end)=[];
MAGNETI_GEO.mat_number(1:end)=[];
Nring_used=length(SHIMMING_OPT.ring_to_be_used);   %Number of used rings in this opt run
ii=0;
for Nu=1:Nring_used
  Nr=SHIMMING_OPT.ring_to_be_used(Nu);   %Number of used Ring
  for Ns=1:SHIMMING_GEO.Nsector_Shim
    for Nm=1:SHIMMING_OPT.PM_add_per_sector
      ii=ii+1;
      if mag_yes_no(ii,1)==1
        ipos=pos_in_sector(ii,1);
        angolo=SHIMMING_GEO.angoloslot(Ns,ipos);
        if SHIMMING_GEO.MatrixShimming(Nr,Ns,ipos) == 0
          MAGNETI_GEO.Ndipoli_tot=MAGNETI_GEO.Ndipoli_tot+1;
          MAGNETI_GEO.xdip(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_GEO.Raggio_shim*cos(angolo/180*pi);
          MAGNETI_GEO.ydip(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_GEO.Raggio_shim*sin(angolo/180*pi);
          MAGNETI_GEO.zdip(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_GEO.ring_Zposition(Nr);
          MAGNETI_GEO.angle(MAGNETI_GEO.Ndipoli_tot,1)=2*angolo+dir_mag_shim(ii,1);
          MAGNETI_GEO.cubo_dimU2(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.cubo_dimU2_shim_ref;
          MAGNETI_GEO.cubo_dimV2(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.cubo_dimV2_shim_ref;
          MAGNETI_GEO.cubo_dimW2(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.cubo_dimW2_shim_ref;
          MAGNETI_GEO.mat_number(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.PM_material;
          SHIMMING_GEO.MatrixShimming(Nr,Ns,ipos)=SHIMMING_OPT.Run;
        end
        if SHIMMING_OPT.Zaxis_symmetry
          Nr2=Nr+SHIMMING_GEO.Nring_tot/2;
          if SHIMMING_GEO.MatrixShimming(Nr2,Ns,ipos) == 0
            MAGNETI_GEO.Ndipoli_tot=MAGNETI_GEO.Ndipoli_tot+1;
            MAGNETI_GEO.xdip(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_GEO.Raggio_shim*cos(angolo/180*pi);
            MAGNETI_GEO.ydip(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_GEO.Raggio_shim*sin(angolo/180*pi);
            MAGNETI_GEO.zdip(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_GEO.ring_Zposition(Nr2);
            MAGNETI_GEO.angle(MAGNETI_GEO.Ndipoli_tot,1)=2*angolo+dir_mag_shim(ii,1);
            MAGNETI_GEO.cubo_dimU2(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.cubo_dimU2_shim_ref;
            MAGNETI_GEO.cubo_dimV2(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.cubo_dimV2_shim_ref;
            MAGNETI_GEO.cubo_dimW2(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.cubo_dimW2_shim_ref;
            MAGNETI_GEO.mat_number(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.PM_material;
            SHIMMING_GEO.MatrixShimming(Nr2,Ns,ipos)=SHIMMING_OPT.Run;
          end
        end
      end
    end
  end
end
return
end



function deviation_tot_ppm = deviation_in_sphere3(x_input)
global MAGNETI_GEO
global TEMPERATURE
global MAGNETI_STATO
global MATERIALI
global SIMUL_DATA
global SHIMMING_OPT
global POINTS

nvars=size(x_input,2);
nnn=nvars/3;

mag_yes_no=transpose(x_input(1,1:nnn));
pos_in_sector=transpose(x_input(1,nnn+1:2*nnn));
dir_mag_shim=transpose(x_input(1,2*nnn+1:nvars));

%------------------------------------------------------
%  Input geometrico configurazione shimming
%------------------------------------------------------
build_dipole_shim(mag_yes_no,pos_in_sector,dir_mag_shim);

if MAGNETI_GEO.Ndipoli_tot > 0
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
      recompute,MAGNETI_GEO,MAGNETI_STATO,MATERIALI,POINTS,L,U,P,Tnoto2,SHIMMING_OPT.Mu0,output_JHmag);
  B5 = util.estrai_output(SIMUL_DATA,BFIELD);
else
  B5=zeros(length(SHIMMING_OPT.B5),1);
end    
B5_tot=SHIMMING_OPT.B5+B5;
[~,~,~,DEV1_ppm,DEV2_ppm] = util.variability(B5_tot);

if strcmpi(SHIMMING_OPT.DevOpt,'DEV1')
  deviation_tot_ppm=DEV1_ppm;
elseif strcmpi(SHIMMING_OPT.DevOpt,'DEV2')
  deviation_tot_ppm=DEV2_ppm;
end

fprintf('Added magnets: %d - deviation_tot_ppm: %.1f\n',MAGNETI_GEO.Ndipoli_tot,deviation_tot_ppm);


return
end