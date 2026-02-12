
classdef read_labels

methods (Static)


function [scale] = readUnit(toml_data)
scale=1.0;
if isfield(toml_data,'unit')
  if isfield(toml_data.unit,'scale')
    scale = toml_data.unit.scale;
  end
end
return
end

function [ierr,MAGNETI_GEO] = readPM(toml_data,MAGNETI_GEO,scale)
ierr=0;
if ~isfield(toml_data,'pm')
    fprintf('Field [pm] not present\n');
    ierr=1;
    return
end
if ~isfield(toml_data.pm,'mat_code')
    fprintf('Field [pm.mat_code] not present\n');
    ierr=1;
    return
end
if ~isfield(toml_data.pm,'files')
    fprintf('Field [pm.files] not present\n');
    ierr=1;
    return
end
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

function [ierr,MAGNETI_GEO,MAGNETI_GEOshim] = readPMshim(toml_data,MAGNETI_GEO,MAGNETI_GEOshim,scale)
ierr=0;
if ~isfield(toml_data.pm,'shim')
    return
end
if ~isfield(toml_data.pm.shim,'mat_code')
    fprintf('Field [pm.shim.mat_code] not present\n');
    ierr=1;
    return
end
if ~isfield(toml_data.pm.shim,'files')
    fprintf('Field [pm.shim.files] not present\n');
    ierr=1;
    return
end
MAGNETI_GEOshim.Nlist=length(toml_data.pm.shim.mat_code);
MAGNETI_GEO.Nlist=MAGNETI_GEO.Nlist+MAGNETI_GEOshim.Nlist;
MAGNETI_GEOshim.material=toml_data.pm.shim.mat_code;
MAGNETI_GEO.material=[MAGNETI_GEO.material MAGNETI_GEOshim.material];
MAGNETI_GEOshim.fileGEO=toml_data.pm.shim.files;
MAGNETI_GEO.fileGEO=[MAGNETI_GEO.fileGEO MAGNETI_GEOshim.fileGEO];

for n=MAGNETI_GEO.Nlist-MAGNETI_GEOshim.Nlist+1:MAGNETI_GEO.Nlist
% Read file with shimming magnets geometry
  [MAGNETI_GEO] = read_labels.readMAG_GEO(n,scale,MAGNETI_GEO);
end
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

function [ierr,MAGNETI_GEO] = readMODGEO(toml_data,scale,MAGNETI_GEO)
ierr=0;
if ~isfield(toml_data.pm,'modgeo')
  MAGNETI_GEO.active=false;
  return
end
if ~isfield(toml_data.pm.modgeo,'active')
  MAGNETI_GEO.active=false;
  return
end
MAGNETI_GEO.active=toml_data.pm.modgeo.active;
if ~MAGNETI_GEO.active
  return
end
if ~isfield(toml_data.pm.modgeo,'code')
  fprintf('Field [pm.modgeo.code] not present\n');
  ierr=1;
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
  if ~isfield(toml_data.pm.modgeo,'zpos')
    fprintf('Field [pm.modgeo.zpos] not present\n');
    ierr=1;
    return
  end
  MAGNETI_GEO.zlista=toml_data.pm.modgeo.zpos*scale;
  for nd=1:MAGNETI_GEO.Ndipoli_tot
    Nr=MAGNETI_GEO.group(nd,1);
    if ~isnan(MAGNETI_GEO.zlista(Nr))
      MAGNETI_GEO.zdip(nd,1)=MAGNETI_GEO.zlista(Nr);
    end
  end
end
if angmod
  if ~isfield(toml_data.pm.modgeo,'file')
    fprintf('Field [pm.modgeo.file] not present\n');
    ierr=1;
    return
  end
  MAGNETI_GEO.file_momento_torcente=toml_data.pm.modgeo.file;
  aus=load(MAGNETI_GEO.file_momento_torcente);
  MAGNETI_GEO.momento_torcente=aus.momento_total(:,1);
  clear aus;
  if ~isfield(toml_data.pm.modgeo,'angle')
    fprintf('Field [pm.modgeo.angle] not present\n');
    ierr=1;
    return
  end
  MAGNETI_GEO.d_angle=toml_data.pm.modgeo.angle;
  MAGNETI_GEO.angle=MAGNETI_GEO.angle+MAGNETI_GEO.d_angle.*sign(MAGNETI_GEO.momento_torcente);
end
return
end

function [ierr,MATERIALI] = readMATERIAL(toml_data,MATERIALI)
ierr=0;
if ~isfield(toml_data,'material')
  fprintf('Field [material] not present\n');
  ierr=1;
  return
end
if ~isfield(toml_data.material,'tref')
  fprintf('Field [material.tref] not present\n');
  ierr=1;
  return
end
MATERIALI.Nlist=length(toml_data.material.tref);
MATERIALI.Lista=transpose(1:MATERIALI.Nlist);

if isfield(toml_data.material,'JHcurve')
  MATERIALI.JHcurve=toml_data.material.JHcurve;
else
  MATERIALI.JHcurve=false;
end
if MATERIALI.JHcurve
  if ~isfield(toml_data.material,'jvt')
    fprintf('Field [material.jvt] not present\n');
    ierr=1;
    return
  end
  if ~isfield(toml_data.material,'hvt1')
    fprintf('Field [material.hvt1] not present\n');
    ierr=1;
    return
  end
  if ~isfield(toml_data.material,'hvt2')
    fprintf('Field [material.hvt2] not present\n');
    ierr=1;
    return
  end
  if ~isfield(toml_data.material,'files')
    fprintf('Field [material.files] not present\n');
    ierr=1;
    return
  end
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
  if ~isfield(toml_data.material,'MuR')
    fprintf('Field [material.MuR] not present\n');
    ierr=1;
    return
  end
  if ~isfield(toml_data.material,'Jresidua')
    fprintf('Field [material.Jresidua] not present\n');
    ierr=1;
    return
  end
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
HJcurve = fscanf(IDfileNL,'%f,%f',sizeA);
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

function [ierr,TEMPERATURE] = readTEMPERATURE(toml_data,TEMPERATURE)
ierr=0;
if ~isfield(toml_data,'temperature')
   TEMPERATURE.TempVar=0;
   return
end
if ~isfield(toml_data.temperature,'type')
    fprintf('Field [temperature.type] not present\n');
    ierr=1;
    return
end
TEMPERATURE.TempVar=toml_data.temperature.type;
if TEMPERATURE.TempVar == 1
  if ~isfield(toml_data.temperature,'value')
    fprintf('Field [temperature.value] not present\n');
    ierr=1;
    return
  end
  TEMPERATURE.Tactual=toml_data.temperature.value;
elseif TEMPERATURE.TempVar == 2
  if ~isfield(toml_data.temperature,'file')
    fprintf('Field [temperature.file] not present\n');
    ierr=1;
    return
  end
  TEMPERATURE.file=toml_data.temperature.file;
end
return
end

function [ierr,POINTS] = readPOINTS(toml_data,POINTS,scale)
ierr=0;
if ~isfield(toml_data,'points')
    fprintf('Field [points] not present\n');
    ierr=1;
    return
end
if ~isfield(toml_data.points,'file')
    fprintf('Field [points.file] not present\n');
    ierr=1;
    return
end
POINTS.file=toml_data.points.file;

if isfield(toml_data.points,'typepoint')
  POINTS.typepoint=toml_data.points.typepoint;
else
  POINTS.typepoint='V';
end

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

function [ierr,SIMUL_DATA] = readSIMUL(toml_data,SIMUL_DATA)
ierr=0;
if ~isfield(toml_data,'simul')
    fprintf('Field [simul] not present\n');
    ierr=1;
    return
end
if ~isfield(toml_data.simul,'type')
    fprintf('Field [simul.type] not present\n');
    ierr=1;
    return
end
if ~isfield(toml_data.simul,'component')
    fprintf('Field [simul.component] not present\n');
    ierr=1;
    return
end
SIMUL_DATA.Type=upper(toml_data.simul.type);
SIMUL_DATA.Bcomponent=toml_data.simul.component;
% Read type of magnet model
if isfield(toml_data.simul,'reaction')
   SIMUL_DATA.reaction=toml_data.simul.reaction;
else
   SIMUL_DATA.reaction=true;
end
if isfield(toml_data.simul,'NL')
    SIMUL_DATA.NL=toml_data.simul.NL;
else
    SIMUL_DATA.NL=true;
end
return
end

function [ierr,TORQUE] = readTORQUE(toml_data,TORQUE)
ierr=0;
if isfield(toml_data,'torque')
  if ~isfield(toml_data.torque,'compute')
    fprintf('Field [torque.compute] not present\n');
    ierr=1;
    return
  end
  TORQUE.flag=toml_data.torque.compute;
else
  TORQUE.flag=false;
end    
if TORQUE.flag
  if isfield(toml_data.torque,'all')
    TORQUE.all=toml_data.torque.all;
  else
    TORQUE.all=false;
  end
  TORQUE.OutputFile=toml_data.torque.outputfile;
  if ~isfield(toml_data.torque,'outputfile')
    fprintf('Field [torque.outputfile] not present\n');
    ierr=1;
    return
  end
  TORQUE.OutputFile=toml_data.torque.outputfile;
end
return
end

function [ierr,MONTECARLO] = readMONTECARLO(toml_data,MONTECARLO,scale)
ierr=0;
if ~isfield(toml_data,'montecarlo')
    fprintf('Field [montecarlo] not present\n');
    ierr=1;
    return
end
if ~isfield(toml_data.montecarlo,'extractions')
    fprintf('Field [montecarlo.extractions] not present\n');
    ierr=1;
    return
end
MONTECARLO.nround=toml_data.montecarlo.extractions;
if isfield(toml_data.montecarlo,'JHmaterial')
  if ~isfield(toml_data.montecarlo.JHmaterial,'active')
      fprintf('Field [montecarlo.JHmaterial.active] not present\n');
      ierr=1;
      return
  end
  if ~isfield(toml_data.montecarlo.JHmaterial,'value')
      fprintf('Field [montecarlo.JHmaterial.value] not present\n');
      ierr=1;
      return
  end
  if ~isfield(toml_data.montecarlo.JHmaterial,'distribution')
      fprintf('Field [montecarlo.JHmaterial.distribution] not present\n');
      ierr=1;
      return
  end
  MONTECARLO.JHmaterialFlag=toml_data.montecarlo.JHmaterial.active;
  MONTECARLO.JHmaterialVar=toml_data.montecarlo.JHmaterial.value;
  MONTECARLO.JHmaterialDistribution=toml_data.montecarlo.JHmaterial.distribution;
else
  MONTECARLO.JHmaterialFlag=false;
end

if isfield(toml_data.montecarlo,'Jmagnet')
  if ~isfield(toml_data.montecarlo.Jmagnet,'active')
      fprintf('Field [montecarlo.Jmagnet.active] not present\n');
      ierr=1;
      return
  end
  if ~isfield(toml_data.montecarlo.Jmagnet,'value')
      fprintf('Field [montecarlo.Jmagnet.value] not present\n');
      ierr=1;
      return
  end
  if ~isfield(toml_data.montecarlo.Jmagnet,'distribution')
      fprintf('Field [montecarlo.Jmagnet.distribution] not present\n');
      ierr=1;
      return
  end
  MONTECARLO.JmagnetFlag=toml_data.montecarlo.Jmagnet.active;
  MONTECARLO.JmagnetVar=toml_data.montecarlo.Jmagnet.value;
  MONTECARLO.JmagnetDistribution=toml_data.montecarlo.Jmagnet.distribution;
else
  MONTECARLO.JmagnetFlag=false;
end

if isfield(toml_data.montecarlo,'Angle')
  if ~isfield(toml_data.montecarlo.Angle,'active')
      fprintf('Field [montecarlo.Angle.active] not present\n');
      ierr=1;
      return
  end
  if ~isfield(toml_data.montecarlo.Angle,'value')
      fprintf('Field [montecarlo.Angle.value] not present\n');
      ierr=1;
      return
  end
  if ~isfield(toml_data.montecarlo.Angle,'distribution')
      fprintf('Field [montecarlo.Angle.distribution] not present\n');
      ierr=1;
      return
  end
  MONTECARLO.AngleFlag=toml_data.montecarlo.Angle.active;
  MONTECARLO.AngleVar=toml_data.montecarlo.Angle.value;
  MONTECARLO.AngleDistribution=toml_data.montecarlo.Angle.distribution;
else
  MONTECARLO.AngleFlag=false;
end

if isfield(toml_data.montecarlo,'Zring')
  if ~isfield(toml_data.montecarlo.Zring,'active')
      fprintf('Field [montecarlo.Zring.active] not present\n');
      ierr=1;
      return
  end
  if ~isfield(toml_data.montecarlo.Zring,'value')
      fprintf('Field [montecarlo.Zring.value] not present\n');
      ierr=1;
      return
  end
  if ~isfield(toml_data.montecarlo.Zring,'distribution')
      fprintf('Field [montecarlo.Zring.distribution] not present\n');
      ierr=1;
      return
  end
  MONTECARLO.ZringFlag=toml_data.montecarlo.Zring.active;
  MONTECARLO.ZRingVar=toml_data.montecarlo.Zring.value*scale;
  MONTECARLO.ZRingDistribution=toml_data.montecarlo.Zring.distribution;
else
  MONTECARLO.ZringFlag=false;
end

if isfield(toml_data.montecarlo,'XMagnets')
  if ~isfield(toml_data.montecarlo.XMagnets,'active')
      fprintf('Field [montecarlo.XMagnets.active] not present\n');
      ierr=1;
      return
  end
  if ~isfield(toml_data.montecarlo.XMagnets,'value')
      fprintf('Field [montecarlo.XMagnets.value] not present\n');
      ierr=1;
      return
  end
  if ~isfield(toml_data.montecarlo.XMagnets,'distribution')
      fprintf('Field [montecarlo.XMagnets.distribution] not present\n');
      ierr=1;
      return
  end
  MONTECARLO.XMagnetsFlag=toml_data.montecarlo.XMagnets.active;
  MONTECARLO.XMagnetsVar=toml_data.montecarlo.XMagnets.value*scale;
  MONTECARLO.XMagnetsDistribution=toml_data.montecarlo.XMagnets.distribution;
else
  MONTECARLO.XMagnetsFlag=false;
end
    
if isfield(toml_data.montecarlo,'YMagnets')
  if ~isfield(toml_data.montecarlo.YMagnets,'active')
      fprintf('Field [montecarlo.YMagnets.active] not present\n');
      ierr=1;
      return
  end
  if ~isfield(toml_data.montecarlo.YMagnets,'value')
      fprintf('Field [montecarlo.YMagnets.value] not present\n');
      ierr=1;
      return
  end
  if ~isfield(toml_data.montecarlo.YMagnets,'distribution')
      fprintf('Field [montecarlo.YMagnets.distribution] not present\n');
      ierr=1;
      return
  end
  MONTECARLO.YMagnetsFlag=toml_data.montecarlo.YMagnets.active;
  MONTECARLO.YMagnetsVar=toml_data.montecarlo.YMagnets.value*scale;
  MONTECARLO.YMagnetsDistribution=toml_data.montecarlo.YMagnets.distribution;
else
  MONTECARLO.YMagnetsFlag=false;
end

if isfield(toml_data.montecarlo,'ZMagnets')
  if ~isfield(toml_data.montecarlo.ZMagnets,'active')
      fprintf('Field [montecarlo.ZMagnets.active] not present\n');
      ierr=1;
      return
  end
  if ~isfield(toml_data.montecarlo.ZMagnets,'value')
      fprintf('Field [montecarlo.ZMagnets.value] not present\n');
      ierr=1;
      return
  end
  if ~isfield(toml_data.montecarlo.ZMagnets,'distribution')
      fprintf('Field [montecarlo.ZMagnets.distribution] not present\n');
      ierr=1;
      return
  end
  MONTECARLO.ZMagnetsFlag=toml_data.montecarlo.ZMagnets.active;
  MONTECARLO.ZMagnetsVar=toml_data.montecarlo.ZMagnets.value*scale;
  MONTECARLO.ZMagnetsDistribution=toml_data.montecarlo.ZMagnets.distribution;
else
  MONTECARLO.ZMagnetsFlag=false;
end
return
end

end   %Method

end   %Classdef