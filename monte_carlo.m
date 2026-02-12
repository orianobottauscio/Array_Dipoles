classdef monte_carlo

methods (Static)


function [x_output,evaluation_optMC] = MC_opt(f,x_input,N_estrazioni,Var_input,nvars)
a=-Var_input;
b=Var_input;
evaluation_optMC=1e36;
for n=1:N_estrazioni
  r = a +(b-a)*rand(nvars,1);
  xx=x_input+transpose(r);
  yy = feval(f,xx);
  if yy < evaluation_optMC
    evaluation_optMC=yy;
    x_output=xx;
  end
end
return
end



function [MONTECARLO] = Run_Monte_Carlo(write_every,MONTECARLO,MAGNETI_GEO,MATERIALI,MAGNETI_STATO,SIMUL_DATA,POINTS,Mu0)

if strcmpi(POINTS.typepoint,'V')   %Volume data in DSV
  MONTECARLO.raccolta_bmean=zeros(MONTECARLO.nround,1);
  MONTECARLO.raccolta_bmin=zeros(MONTECARLO.nround,1);
  MONTECARLO.raccolta_bmax=zeros(MONTECARLO.nround,1);
  MONTECARLO.raccolta_DEV1=zeros(MONTECARLO.nround,1);
  MONTECARLO.raccolta_DEV2=zeros(MONTECARLO.nround,1);
elseif strcmpi(POINTS.typepoint,'S')   %Surface data in DSV
  MONTECARLO.raccolta_b0=zeros(MONTECARLO.nround,1);
  MONTECARLO.raccolta_bmin=zeros(MONTECARLO.nround,1);
  MONTECARLO.raccolta_bmax=zeros(MONTECARLO.nround,1);
  MONTECARLO.raccolta_DEV3=zeros(MONTECARLO.nround,1);
end

% Save the original values of magnet geometric data
xdip_SAVE=MAGNETI_GEO.xdip(1:MAGNETI_GEO.Ndipoli_tot,1);
ydip_SAVE=MAGNETI_GEO.ydip(1:MAGNETI_GEO.Ndipoli_tot,1);
zdip_SAVE=MAGNETI_GEO.zdip(1:MAGNETI_GEO.Ndipoli_tot,1);
angle_SAVE=MAGNETI_GEO.angle(1:MAGNETI_GEO.Ndipoli_tot,1);

first=true;
% Initialization
L=0;
U=0;
P=0;
Tnoto2=0;
output_JHmag=false;
%
for nnr = 1:MONTECARLO.nround
  if mod(nnr,write_every) == 0
    fprintf('N. round: %d over %d\n',nnr,MONTECARLO.nround);
  end

  recompute=false;
%Recover the original values of magnet geometric data
  MAGNETI_GEO.xdip(1:MAGNETI_GEO.Ndipoli_tot,1)=xdip_SAVE;
  MAGNETI_GEO.ydip(1:MAGNETI_GEO.Ndipoli_tot,1)=ydip_SAVE;
  MAGNETI_GEO.zdip(1:MAGNETI_GEO.Ndipoli_tot,1)=zdip_SAVE;
  MAGNETI_GEO.angle(1:MAGNETI_GEO.Ndipoli_tot,1)=angle_SAVE;
%  
%...................JH Characteristics.................................
%
  if MONTECARLO.JHmaterialFlag
% Estraction to correct the JH magnetic curves
    if strcmp(MONTECARLO.JHmaterialDistribution,'Uniform')
      vmin=1.0-MONTECARLO.JHmaterialVar;
      vmax=1.0+MONTECARLO.JHmaterialVar;
      for n=1:MATERIALI.Nlist
        mat_number=MATERIALI.Lista(n,1);
        Jresidua_relativo(mat_number,1)=random(MONTECARLO.JHmaterialDistribution,vmin,vmax);
      end
    elseif strcmp(MONTECARLO.JHmaterialDistribution,'Normal')
      for n=1:MATERIALI.Nlist
        mat_number=MATERIALI.Lista(n,1);
        Jresidua_relativo(mat_number,1)=random(MONTECARLO.JHmaterialDistribution,1.0,MONTECARLO.JHmaterialVar);
      end
    end
  else
    for n=1:MATERIALI.Nlist
      mat_number=MATERIALI.Lista(n,1);
      Jresidua_relativo(mat_number,1)=1.0;   %Not active therefore JH curves are not rescaled
    end
  end
% Set Jr and Hc on the basis of p.u. input value and magnet temperature
  [MAGNETI_STATO] = dipoles_sub.assegna_valori_dipoles(MAGNETI_GEO,MATERIALI,MAGNETI_STATO,Jresidua_relativo);
%
%...................Jr single magnet.................................
%
  if MONTECARLO.JmagnetFlag
% Estraction Jr of each single magnet
    for nd=1:MAGNETI_GEO.Ndipoli_tot
      Jr_random=random(MONTECARLO.JmagnetDistribution,MAGNETI_STATO.Jresidua_magneti(nd,1),MONTECARLO.JmagnetVar);
      MAGNETI_STATO.Jresidua_magneti(nd,1)=Jr_random;
      MAGNETI_STATO.Hc_magneti(nd,1)=MAGNETI_STATO.Hc_magneti(nd,1)*Jr_random/MAGNETI_STATO.Jresidua_magneti(nd,1);  %Rescale Hc on the basis of Jr random value
    end
  end
%  
%...................Magnet angle.................................
%
  if MONTECARLO.AngleFlag
% Estraction of the angle of each single magnet
    dangle=zeros(MAGNETI_GEO.Ndipoli_tot,1);
    for nd=1:MAGNETI_GEO.Ndipoli_tot
      dangle(nd,1)=random(MONTECARLO.AngleDistribution,-MONTECARLO.AngleVar,MONTECARLO.AngleVar);
    end
    MAGNETI_GEO.angle=MAGNETI_GEO.angle+dangle(1:MAGNETI_GEO.Ndipoli_tot,1);
    recompute=true;
  end
%  
%...................Position of magnet groups along their axis...........................
%
  if MONTECARLO.ZringFlag
% Estraction of the z position of each magnet group
    dZVar=zeros(MAGNETI_GEO.Nlist,1);
    for Nr=1:MAGNETI_GEO.Nlist
      dZVar(Nr,1)=random(MONTECARLO.ZRingDistribution,-MONTECARLO.ZRingVar,MONTECARLO.ZRingVar);
    end
    for nd=1:MAGNETI_GEO.Ndipoli_tot
      Nr=MAGNETI_GEO.group(nd,1);
      MAGNETI_GEO.zdip(nd,1)=MAGNETI_GEO.zdip(nd,1)+dZVar(Nr,1);
    end   
    recompute=true;
  end
%  
%.............Position X of each single magnet.........................
%
  if MONTECARLO.XMagnetsFlag
% Estraction of the x position of each single magnet
    dVar=zeros(MAGNETI_GEO.Ndipoli_tot,1);
    for nd=1:MAGNETI_GEO.Ndipoli_tot
      dVar(nd,1)=random(MONTECARLO.XMagnetsDistribution,-MONTECARLO.XMagnetsVar,MONTECARLO.XMagnetsVar);
    end
    MAGNETI_GEO.xdip=MAGNETI_GEO.xdip+dVar(1:MAGNETI_GEO.Ndipoli_tot,1);
    recompute=true;
  end
%  
%.............Position Y of each single magnet.........................
%
  if MONTECARLO.YMagnetsFlag
% Estraction of the y position of each single magnet
    dVar=zeros(MAGNETI_GEO.Ndipoli_tot,1);
    for nd=1:MAGNETI_GEO.Ndipoli_tot
      dVar(nd,1)=random(MONTECARLO.YMagnetsDistribution,-MONTECARLO.YMagnetsVar,MONTECARLO.YMagnetsVar);
    end
    MAGNETI_GEO.ydip=MAGNETI_GEO.ydip+dVar(1:MAGNETI_GEO.Ndipoli_tot,1);
    recompute=true;
  end
%  
%.............Position Z of each single magnet.........................
%
  if MONTECARLO.ZMagnetsFlag
% Estraction of the z position of each single magnet
    dVar=zeros(MAGNETI_GEO.Ndipoli_tot,1);
    for nd=1:MAGNETI_GEO.Ndipoli_tot
      dVar(nd,1)=random(MONTECARLO.ZMagnetsDistribution,-MONTECARLO.ZMagnetsVar,MONTECARLO.ZMagnetsVar);
    end
    MAGNETI_GEO.zdip=MAGNETI_GEO.zdip+dVar(1:MAGNETI_GEO.Ndipoli_tot,1);
    recompute=true;
  end
%
  [L,U,P,Tnoto2,BFIELD,MAGNETI_STATO] = dipoles_sub.output_dipoles(first,SIMUL_DATA.reaction,...
      SIMUL_DATA.NL,recompute,MAGNETI_GEO,MAGNETI_STATO,MATERIALI,POINTS,L,U,P,Tnoto2,Mu0,output_JHmag);
  first=false;
  B5 = util.estrai_output(SIMUL_DATA,BFIELD);
  [MONTECARLO] = monte_carlo.get_results(nnr,MONTECARLO,B5,POINTS);
end

return
end


function   [MONTECARLO] = get_results(nnr,MONTECARLO,B5,POINTS)

if strcmpi(POINTS.typepoint,'V')   %Volume data in DSV
  [Bmin,Bmax,Bmean,DEV1_ppm,DEV2_ppm] = util.variability(B5);
  MONTECARLO.raccolta_bmean(nnr,1)=Bmean*1000;
  MONTECARLO.raccolta_bmin(nnr,1)=Bmin*1000;
  MONTECARLO.raccolta_bmax(nnr,1)=Bmax*1000;
  MONTECARLO.raccolta_DEV1(nnr,1)=DEV1_ppm;
  MONTECARLO.raccolta_DEV2(nnr,1)=DEV2_ppm;
elseif strcmpi(POINTS.typepoint,'S')   %Surface data in DSV
  [Bmin,Bmax,B0,DEV3_ppm] = util.variability3(B5);
  MONTECARLO.raccolta_b0(nnr,1)=B0*1000;
  MONTECARLO.raccolta_bmin(nnr,1)=Bmin*1000;
  MONTECARLO.raccolta_bmax(nnr,1)=Bmax*1000;
  MONTECARLO.raccolta_DEV3(nnr,1)=DEV3_ppm;
end
    
return
end



end   %Method

end   %Classdef