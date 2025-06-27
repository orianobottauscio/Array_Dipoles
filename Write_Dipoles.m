function Write_Dipoles(varargin)
%-------------------------------------------------------------------------
% Write magnet data into .txt file
%
% Author: O. Bottauscio (first version: 2025)
%-------------------------------------------------------------------------
Code_Version='1.0';

fprintf('Code: Write_Dipoles, Version: %s\n',Code_Version);
DirRun=pwd;   %Directory di run

% Input from prompt line
n_input=size(varargin,2);
if n_input >= 1
  input_data = varargin{1};
else
  fprintf('Input file missing\n')
end
if ~contains(input_data,'.mat')
  input_data=append(input_data,'.mat');
end

output_data=erase(input_data,'.mat');
output_data=append(output_data,'.txt');

aus=load(input_data);
MAGNETI_GEO=aus.MAGNETI_GEO;

if ~all(MAGNETI_GEO.cubo_dimU2(:,1)==MAGNETI_GEO.cubo_dimU2(end,1)) ||...
   ~all(MAGNETI_GEO.cubo_dimV2(:,1)==MAGNETI_GEO.cubo_dimV2(end,1)) ||...
   ~all(MAGNETI_GEO.cubo_dimW2(:,1)==MAGNETI_GEO.cubo_dimW2(end,1))
  fprintf('Magnets having different size\n')
  return
end

scale=1000;   %From m to mm
%-------Ricostruzione e memorizzazione posizioni magneti shimming--------------------
fileSH=append(DirRun,'\',output_data);
fid_SH=fopen(fileSH,'w');
fprintf(fid_SH,'%d,%d,%d\n',MAGNETI_GEO.cubo_dimU2(1)*scale,...
      MAGNETI_GEO.cubo_dimV2(1)*scale,MAGNETI_GEO.cubo_dimW2(1)*scale);
fprintf(fid_SH,'%d\n',MAGNETI_GEO.Ndipoli_tot);
%
for ii=1:MAGNETI_GEO.Ndipoli_tot
  fprintf(fid_SH,'%.2f,%.2f,%.2f,%.2f\n',MAGNETI_GEO.xdip(ii,1)*scale,MAGNETI_GEO.ydip(ii,1)*scale,...
                                 MAGNETI_GEO.zdip(ii,1)*scale,MAGNETI_GEO.angle(ii,1));
end
fclose(fid_SH);
fprintf('Geometrical data of shimming magnets saved in file: %s\n',fileSH);
end


