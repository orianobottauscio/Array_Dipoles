
classdef read_labels

methods (Static)


function [MONTECARLO] = readMONTECARLO(toml_data,MONTECARLO)

MONTECARLO.nround=toml_data.montecarlo.estraction;

if isfield(toml_data.montecarlo,'JHmaterial')
  MONTECARLO.JHmaterialFlag=toml_data.montecarlo.JHmaterial.active;
  MONTECARLO.JHmaterialVar=toml_data.montecarlo.JHmaterial.value;
  MONTECARLO.JHmaterialDistribution=toml_data.montecarlo.JHmaterial.distribution;
else
  MONTECARLO.Jresidua_relativoFlag=false;
end

if isfield(toml_data.montecarlo,'Jmagnet')
  MONTECARLO.JmagnetFlag=toml_data.montecarlo.Jmagnet.active;
  MONTECARLO.JmagnetVar=toml_data.montecarlo.Jmagnet.value;
  MONTECARLO.JmagnetDistribution=toml_data.montecarlo.Jmagnet.distribution;
else
  MONTECARLO.JresiduaFlag=false;
end

if isfield(toml_data.montecarlo,'Angle')
  MONTECARLO.AngleFlag=toml_data.montecarlo.Angle.active;
  MONTECARLO.AngleVar=toml_data.montecarlo.Angle.value;
  MONTECARLO.AngleDistribution=toml_data.montecarlo.Angle.distribution;
else
  MONTECARLO.AngleFlag=false;
end

if isfield(toml_data.montecarlo,'Zring')
  MONTECARLO.ZringFlag=toml_data.montecarlo.Zring.active;
  MONTECARLO.ZRingVar=toml_data.montecarlo.Zring.value;
  MONTECARLO.ZRingDistribution=toml_data.montecarlo.Zring.distribution;
else
  MONTECARLO.ZringFlag=false;
end

if isfield(toml_data.montecarlo,'XMagnets')
  MONTECARLO.XMagnetsFlag=toml_data.montecarlo.XMagnets.active;
  MONTECARLO.XMagnetsVar=toml_data.montecarlo.XMagnets.value;
  MONTECARLO.XMagnetsDistribution=toml_data.montecarlo.XMagnets.distribution;
else
  MONTECARLO.XMagnetsFlag=false;
end
    
if isfield(toml_data.montecarlo,'YMagnets')
  MONTECARLO.YMagnetsFlag=toml_data.montecarlo.YMagnets.active;
  MONTECARLO.YMagnetsVar=toml_data.montecarlo.YMagnets.value;
  MONTECARLO.YMagnetsDistribution=toml_data.montecarlo.YMagnets.distribution;
else
  MONTECARLO.YMagnetsFlag=false;
end

if isfield(toml_data.montecarlo,'ZMagnets')
  MONTECARLO.ZMagnetsFlag=toml_data.montecarlo.ZMagnets.active;
  MONTECARLO.ZMagnetsVar=toml_data.montecarlo.ZMagnets.value;
  MONTECARLO.ZMagnetsDistribution=toml_data.montecarlo.ZMagnets.distribution;
else
  MONTECARLO.ZMagnetsFlag=false;
end
return
end


function [MAGNETI_GEO] = readMODGEO(toml_data,scale,MAGNETI_GEO)
MAGNETI_GEO.active=toml_data.pm.modgeo.active;
if ~MAGNETI_GEO.active
  return
end
MAGNETI_GEO.modifica = toml_data.pm.modgeo.code;
zmod=false;
if MAGNETI_GEO.modifica == 1 || MAGNETI_GEO.modifica == 3
  zmod=true;
end
angmod=false;
if MAGNETI_GEO.modifica == 2 || MAGNETI_GEO.modifica == 3
  angmod=true;
end
if zmod
  MAGNETI_GEO.zlista=toml_data.pm.modgeo.zpos*scale;
  for nd=1:MAGNETI_GEO.Ndipoli_tot
    Nr=MAGNETI_GEO.group(nd,1);
    MAGNETI_GEO.zdip(nd,1)=MAGNETI_GEO.zlista(Nr);
  end
end
if angmod
  MAGNETI_GEO.file_momento_torcente=toml_data.pm.modgeo.file;
  aus=load(MAGNETI_GEO.file_momento_torcente);
  MAGNETI_GEO.momento_torcente=aus.momento_total;
  clear aus;
  MAGNETI_GEO.d_angle=toml_data.pm.modgeo.angle;
  MAGNETI_GEO.angle=MAGNETI_GEO.angle+MAGNETI_GEO.d_angle.*sign(MAGNETI_GEO.momento_torcente);
end
return
end



function [SIMUL_DATA] = readSIMUL(toml_data,SIMUL_DATA)
SIMUL_DATA.Type=upper(toml_data.simul.type);
SIMUL_DATA.Bcomponent=toml_data.simul.component;
if isfield(toml_data.simul,'outputfile')
    SIMUL_DATA.OutputFile=toml_data.simul.outputfile;
else
    SIMUL_DATA.OutputFile="standard";
end
% Read type of magnet model
SIMUL_DATA.reaction=toml_data.simul.reaction;
SIMUL_DATA.NL=toml_data.simul.NL;
return
end


function [TORQUE] = readTORQUE(toml_data,TORQUE)
if isfield(toml_data.simul,'torque')
    TORQUE.flag=toml_data.simul.torque.compute;
    TORQUE.OutputFile=toml_data.simul.torque.outputfile;
else
    TORQUE.flag=false;
end    
return
end


function [TEMPERATURE] = readTEMPERATURE(toml_data,TEMPERATURE)
TEMPERATURE.TempVar=toml_data.temperature.type;
if TEMPERATURE.TempVar == 1
  TEMPERATURE.Tactual=toml_data.temperature.value;
elseif TEMPERATURE.TempVar == 2
  TEMPERATURE.file=toml_data.temperature.file;
end
return
end

function [MAGNETI_GEO] = readPM(toml_data,MAGNETI_GEO,scale)
MAGNETI_GEO.Nlist=length(toml_data.pm.mat_code);
MAGNETI_GEO.Ndipoli_tot=0;
MAGNETI_GEO.material=toml_data.pm.mat_code;
MAGNETI_GEO.fileGEO=toml_data.pm.files;

for n=1:MAGNETI_GEO.Nlist
% Read file with magnets geometry
  [MAGNETI_GEO] = read_labels.readMAG_GEO(n,scale,MAGNETI_GEO);
end
return
end


function [MATERIALI] = readMATERIAL(toml_data,MATERIALI)
MATERIALI.Nlist=length(toml_data.material.tref);
MATERIALI.Lista=transpose(1:MATERIALI.Nlist);
MATERIALI.JHcurve=toml_data.material.JHcurve;
if strcmpi(MATERIALI.JHcurve,'YES')
  MATERIALI.Tref=transpose(toml_data.material.tref);
  MATERIALI.JvT=transpose(toml_data.material.jvt);
  MATERIALI.HvT(:,1)=transpose(toml_data.material.hvt1);
  MATERIALI.HvT(:,2)=transpose(toml_data.material.hvt2);
  MATERIALI.fileHJ=toml_data.material.files;
  for n=1:MATERIALI.Nlist
  % Read file with JH characteristic
    [MATERIALI] = read_labels.readHJ(n,MATERIALI);
  end
else
  MATERIALI.MuR=transpose(toml_data.material.MuR);
  MATERIALI.JrNL=transpose(toml_data.material.Jresidua);
  MATERIALI.HJcurve_all=0;    %created only for compatibility
  MATERIALI.HJcurve_all_points=0;
end
return
end

function [MATERIALI] = readHJ(mat_number,MATERIALI)
fileNL=string(MATERIALI.fileHJ(mat_number));
IDfileNL = fopen(fileNL,'r');
sizeA = [2 Inf];
HJcurve = fscanf(IDfileNL,'%f %f',sizeA);
HJcurve=transpose(HJcurve);
fclose(IDfileNL);
npNL=size(HJcurve,1);
for i=1:npNL-1
  if HJcurve(i,1)*HJcurve(i+1,1) <= 0
      ii=i;
      break
  end
end
DJ=HJcurve(ii+1,2)-HJcurve(ii,2);
DH=HJcurve(ii+1,1)-HJcurve(ii,1);
MATERIALI.JrNL(mat_number,1)=HJcurve(ii,2)-HJcurve(ii,1)*DJ/DH;   %Valore di Jr riferito alla caratteristica lineare
for i=1:npNL-1
  if HJcurve(i,2) <= 0
      ii=i;
      break
  end
end
MATERIALI.HcNL(mat_number,1)=HJcurve(ii,1);   %Valore di HC riferito alla caratteristica lineare
MATERIALI.HSat(mat_number,1)=HJcurve(end,1);
MATERIALI.JSat(mat_number,1)=HJcurve(end,2);
MATERIALI.MuR(mat_number,1)=(MATERIALI.JSat(mat_number,1)-MATERIALI.JrNL(mat_number,1))/...
    MATERIALI.HSat(mat_number,1)/(4e-7*pi);
MATERIALI.HJcurve_all_points(mat_number,1)=size(HJcurve,1);
MATERIALI.HJcurve_all(1:size(HJcurve,1),1:2,mat_number)=HJcurve;
return
end


function [POINTS] = readPOINTS(toml_data,POINTS,scale)
POINTS.file=toml_data.points.file;
IDfilePoint = fopen(POINTS.file,'r');
sizeA = [3 Inf];
CPoints = fscanf(IDfilePoint,'%f,%f,%f',sizeA);
CPoints=transpose(CPoints);
fclose(IDfilePoint);
POINTS.Npoint=size(CPoints,1);
POINTS.xp=CPoints(:,1)*scale;
POINTS.yp=CPoints(:,2)*scale;
POINTS.zp=CPoints(:,3)*scale;
return
end

function [MAGNETI_GEO] = readMAG_GEO(n,scale,MAGNETI_GEO)
%------------------------------------------------------------------
% Readling geometrical data of magnets
%------------------------------------------------------------------
fileGEO=string(MAGNETI_GEO.fileGEO(n));
fileID = fopen(fileGEO,'r');
riga=fgetl(fileID);
valori=textscan(riga,'%f %f %f','Delimiter',{' ',','});
dimU=valori{1}*scale;
dimV=valori{2}*scale;
dimW=valori{3}*scale;
%
riga=fgetl(fileID);
valori=textscan(riga,'%d');
Ndipoli=valori{1};
for nd=1:Ndipoli
   riga=fgetl(fileID);
   valori=textscan(riga,'%f %f %f %f','Delimiter',{' ',','});
   MAGNETI_GEO.Ndipoli_tot=MAGNETI_GEO.Ndipoli_tot+1;
   MAGNETI_GEO.xdip(MAGNETI_GEO.Ndipoli_tot,1)=valori{1}*scale;
   MAGNETI_GEO.ydip(MAGNETI_GEO.Ndipoli_tot,1)=valori{2}*scale;
   MAGNETI_GEO.zdip(MAGNETI_GEO.Ndipoli_tot,1)=valori{3}*scale;
   MAGNETI_GEO.angle(MAGNETI_GEO.Ndipoli_tot,1)=valori{4};
   MAGNETI_GEO.cubo_dimU2(MAGNETI_GEO.Ndipoli_tot,1)=dimU;
   MAGNETI_GEO.cubo_dimV2(MAGNETI_GEO.Ndipoli_tot,1)=dimV;
   MAGNETI_GEO.cubo_dimW2(MAGNETI_GEO.Ndipoli_tot,1)=dimW;
   MAGNETI_GEO.mat_number(MAGNETI_GEO.Ndipoli_tot,1)=MAGNETI_GEO.material(n);
   MAGNETI_GEO.group(MAGNETI_GEO.Ndipoli_tot,1)=n;
end
fclose(fileID);
return
end



end   %Method

end   %Classdef