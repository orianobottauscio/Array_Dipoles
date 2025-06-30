function Write_Dipoles(varargin)
%-------------------------------------------------------------------------
% Write magnet data into .txt file
%
% Author: O. Bottauscio (first version: 2025)
%-------------------------------------------------------------------------
Code_Version='1.0';
fprintf('Code: Write_Dipoles, Version: %s\n',Code_Version);

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
aus=load(input_data);
MAGNETI_GEO=aus.MAGNETI_GEO;
scale=1000;   %From m to mm

run=0;
nd_ini=1;
ciclo=true;
while ciclo
  run=run+1;
  for nd=nd_ini+1:MAGNETI_GEO.Ndipoli_tot
    if MAGNETI_GEO.cubo_dimU2(nd,1) ~= MAGNETI_GEO.cubo_dimU2(nd-1,1) ||...
       MAGNETI_GEO.cubo_dimV2(nd,1) ~= MAGNETI_GEO.cubo_dimV2(nd-1,1) ||...
       MAGNETI_GEO.cubo_dimW2(nd,1) ~= MAGNETI_GEO.cubo_dimW2(nd-1,1)

%-------Ricostruzione e memorizzazione posizioni magneti shimming--------------------
      ntot=nd-nd_ini;
      fileSH=append(output_data,'_',num2str(run),'.txt');
      fid_SH=fopen(fileSH,'w');
      fprintf(fid_SH,'%d,%d,%d\n',MAGNETI_GEO.cubo_dimU2(nd_ini)*scale,...
            MAGNETI_GEO.cubo_dimV2(nd_ini)*scale,MAGNETI_GEO.cubo_dimW2(nd_ini)*scale);
      fprintf(fid_SH,'%d\n',ntot);
      %
      for ii=nd_ini:nd-1
        fprintf(fid_SH,'%.2f,%.2f,%.2f,%.2f\n',MAGNETI_GEO.xdip(ii,1)*scale,MAGNETI_GEO.ydip(ii,1)*scale,...
                                       MAGNETI_GEO.zdip(ii,1)*scale,MAGNETI_GEO.angle(ii,1));
      end
      fclose(fid_SH);
      fprintf('Geometrical data of shimming magnets saved in file: %s\n',fileSH);
      nd_ini=nd;
      break
    end
    if nd == MAGNETI_GEO.Ndipoli_tot
        ciclo=false;
    end
  end
end
if run == 1
    fileSH=append(output_data,'.txt');
else
    fileSH=append(output_data,'_',num2str(run),'.txt');
end    
fid_SH=fopen(fileSH,'w');
fprintf(fid_SH,'%d,%d,%d\n',MAGNETI_GEO.cubo_dimU2(nd_ini)*scale,...
      MAGNETI_GEO.cubo_dimV2(nd_ini)*scale,MAGNETI_GEO.cubo_dimW2(nd_ini)*scale);
ntot=MAGNETI_GEO.Ndipoli_tot-nd_ini+1;
fprintf(fid_SH,'%d\n',ntot);
%
for ii=nd_ini:MAGNETI_GEO.Ndipoli_tot
  fprintf(fid_SH,'%.2f,%.2f,%.2f,%.2f\n',MAGNETI_GEO.xdip(ii,1)*scale,MAGNETI_GEO.ydip(ii,1)*scale,...
                                 MAGNETI_GEO.zdip(ii,1)*scale,MAGNETI_GEO.angle(ii,1));
end
fclose(fid_SH);
fprintf('Geometrical data of shimming magnets saved in file: %s\n',fileSH);

end


