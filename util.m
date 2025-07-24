classdef util

methods (Static)

%--------------------------------------------------------------------
% Utility
%--------------------------------------------------------------------
function summary_MC(MONTECARLO)
  fprintf('MC estraction: %d\n',MONTECARLO.nround);
  if MONTECARLO.JHmaterialFlag
    fprintf('JHmaterialVar (pu): %d\n',MONTECARLO.JHmaterialVar);
  end
  if MONTECARLO.JmagnetFlag
    fprintf('JmagnetVar (T): %d\n',MONTECARLO.JmagnetVar);
  end
  if MONTECARLO.AngleFlag
    fprintf('AngleVar (degree): %d\n',MONTECARLO.AngleVar);
  end
  if MONTECARLO.ZringFlag
    fprintf('ZRingVar (mm): %d\n',MONTECARLO.ZRingVar*1000);
  end
  if MONTECARLO.XMagnetsFlag
    fprintf('XMagnetsVar (mm): %d\n',MONTECARLO.XMagnetsVar*1000);
  end
  if MONTECARLO.YMagnetsFlag
    fprintf('YMagnetsVar (mm): %d\n',MONTECARLO.YMagnetsVar*1000);
  end
  if MONTECARLO.ZMagnetsFlag
    fprintf('ZMagnetsVar (mm): %d\n',MONTECARLO.ZMagnetsVar*1000);
  end
return
end

function [yMean,yStd,CI] = statistical_data(y)
cov=95;   %Coverage 95%
t1=(1-cov/100)/2;
t2=1-t1;
[yMean,yStd,CI,~] = util.statistical_value(y,t1,t2);
return
end

function [yMean,yStd,CI,ySEM] = statistical_value(y,t1,t2)

N = length(y);
yMean = mean(y);     % Mean Of All Experiments At Each Value Of ‘x’
yStd = std(y);
ySEM = yStd/sqrt(N);                              % Compute SStandard Error Of The Mean’ Of All Experiments At Each Value Of ‘x’
CI = yMean+yStd*tinv([t1 t2], N-1);                    % Calculate 95% Probability Intervals Of t-Distribution

return
end

function [Bmin,Bmax,Bmean,DEV1_ppm,DEV2_ppm] = variability(B5)
  Bmin=min(abs(B5(:,1)));
  Bmax=max(abs(B5(:,1)));
  Bmean=mean(abs(B5(:,1)));
  Bstd=std(abs(B5(:,1)));
  DEV2_ppm=Bstd/Bmean*1e6;
  DEV1_ppm=(Bmax-Bmin)/Bmean*1e6;
return
end


function B5 = estrai_output(SIMUL_DATA,BFIELD)
if strcmpi(SIMUL_DATA.Bcomponent,'x')
 B5=BFIELD.Bx5;
elseif strcmpi(SIMUL_DATA.Bcomponent,'y')
 B5=BFIELD.By5;
elseif strcmpi(SIMUL_DATA.Bcomponent,'z')
 B5=BFIELD.Bz5;
elseif strcmpi(SIMUL_DATA.Bcomponent,'mag')
 B5=sqrt(BFIELD.Bx5.^2+BFIELD.By5.^2+BFIELD.Bz5.^2);
end
return
end



function Hc = Hc_versus_T(HvT,T,Tref,Href)
Hc = Href.*(1+HvT(1).*(T-Tref)+HvT(2).*(T-Tref).^2);

return
end

function Jr = Jr_versus_T(JvT,T,Tref,Jref)
Jr = Jref.*(1+JvT(1).*(T-Tref));
return
end

function [xpl,ypl,zpl] = global_to_local_par(xp,yp,zp,drif)
%===================================================================
%     Il Punto di coordinate (XP,YP,ZP) nel sistema master
%     viene trasformato nel punto di coordinate (XPL,YPL,ZPL)
%     del sistema locale definito da DRIF
%===================================================================
[phi,the,psi] = util.eultre(drif);
[r] = util.matrot(phi,the,psi);
r=transpose(r);   %matrice trasposta per trasformazione inversa
vaus1(1,:)=xp(:,1)-drif(1,1);
vaus1(2,:)=yp(:,1)-drif(1,2);
vaus1(3,:)=zp(:,1)-drif(1,3);
vaus2=r*vaus1;
xpl(:,1)=vaus2(1,:);
ypl(:,1)=vaus2(2,:);
zpl(:,1)=vaus2(3,:);
return
end

function [xpl,ypl,zpl] = global_to_local(xp,yp,zp,drif)
%===================================================================
%     Il Punto di coordinate (XP,YP,ZP) nel sistema master
%     viene trasformato nel punto di coordinate (XPL,YPL,ZPL)
%     del sistema locale definito da DRIF
%===================================================================
[phi,the,psi] = util.eultre(drif);
[r] = util.matrot(phi,the,psi);
r=transpose(r);   %matrice trasposta per trasformazione inversa
vaus1(1,1)=xp-drif(1,1);
vaus1(2,1)=yp-drif(1,2);
vaus1(3,1)=zp-drif(1,3);
vaus2=r*vaus1;
xpl=vaus2(1);
ypl=vaus2(2);
zpl=vaus2(3);
return
end

function [xp,yp,zp] = local_to_global(xpl,ypl,zpl,drif)
%===================================================================
%     Il Punto di coordinate (XPL,YPL,ZPL) nel sistema locale
%     definito da DRIF viene trasformato nel punto di
%     coordinate (XP,YP,ZP) del sistema master
%===================================================================
[phi,the,psi] = util.eultre(drif);
[r] = util.matrot(phi,the,psi);
vaus1(1,1)=xpl;
vaus1(2,1)=ypl;
vaus1(3,1)=zpl;
vaus2=r*vaus1;
xp=vaus2(1)+drif(1,1);
yp=vaus2(2)+drif(1,2);
zp=vaus2(3)+drif(1,3);
return
end

function [vx,vy,vz] = vector_local_to_global_par(vxl,vyl,vzl,drif)
%===================================================================
%       Dato un vettore di componenti (VXL,VYL,VZL) definite
%       nel sistema locale individuato da DRIF viene trasformato nel
%       vettore di componenti (VX,VY,VZ) del sistema master
%===================================================================
[phi,the,psi] = util.eultre(drif);
[r] = util.matrot(phi,the,psi);
vaus1(1,:)=vxl(:,1);
vaus1(2,:)=vyl(:,1);
vaus1(3,:)=vzl(:,1);
vaus2=r*vaus1;
vx(:,1)=vaus2(1,:);
vy(:,1)=vaus2(2,:);
vz(:,1)=vaus2(3,:);
return
end

function [vx,vy,vz] = vector_local_to_global(vxl,vyl,vzl,drif)
%===================================================================
%       Dato un vettore di componenti (VXL,VYL,VZL) definite
%       nel sistema locale individuato da DRIF viene trasformato nel
%       vettore di componenti (VX,VY,VZ) del sistema master
%===================================================================
[phi,the,psi] = util.eultre(drif);
[r] = util.matrot(phi,the,psi);
vaus1(1,1)=vxl;
vaus1(2,1)=vyl;
vaus1(3,1)=vzl;
vaus2=r*vaus1;
vx=vaus2(1);
vy=vaus2(2);
vz=vaus2(3);
return
end

function [phi,the,psi] = eultre(drif)
%===================================================================================
%     ROUTINE PER DETERMINARE GLI ANGOLI DI EULERO CORRISPONDENTI
%     AD UN SISTEMA DI RIFERIMENTO CARTESIANO DI
%     ORIGINE                   XO,YO,ZO
%     PUNTO SU ASSE ZETA        XZ,YZ,ZZ
%     PUNTO SU ASSE X           XX,YX,ZX
%
%     GLI ANGOLI SONO:
%
%     PHI ROTAZIONE ATTORNO ALL'ASSE Z DEL RIFERIMENTO PRINCIPALE
%     THE ROTAZIONE ATTORNO ALL'ASSE X DEL SISTEMA RUOTATO DI PHI
%     PSI ROTAZIONE ATTORNO ALL'ASSE Z DEL SISTEMA RUOTATO DI PHI E THE
%
%     NOTAZIONE GOLDSTEIN
%===================================================================================
%     TRASLAZIONE NELL'ORIGINE
uz(1:3,1)=drif(2,1:3)-drif(1,1:3);
ux(1:3,1)=drif(3,1:3)-drif(1,1:3);
%     NORMALIZZAZIONE DEI VETTORI PZ E PX
ux=ux/norm(ux);
uz=uz/norm(uz);
%     CALCOLO ANGOLO PHI DALLA DIREZIONE DEL PRODOTTO VETTORE
an(1,1)=-uz(2,1);
an(2,1)=uz(1,1);
an(3,1)=0;
phi = atan2(an(2,1),an(1,1));
%     RICAVA THE DAL PRODOTTO SCALARE DI UZ E VERSORE ASSE Z
eta1=-uz(1,1)*sin(phi)+uz(2,1)*cos(phi);
if eta1 <= 0
    the=acos(uz(3,1));
else
    the=-acos(uz(3,1));
end
%     RICAVA PSI DAL PRODOTTO SCALARE DI UX CON IL VERSORE DELL'ASSE X RUOTATO
an(1,1)=cos(phi);
an(2,1)=sin(phi);
an(3,1)=0;
eta1=ux(1,1)*(-cos(the)*sin(phi))+ux(2,1)*cos(the)*cos(phi)+ux(3,1)*sin(the);
ango=dot(ux,an);
if eta1 >= 0
    psi=acos(ango);
else
    psi=-acos(ango);
end

end

function [r] = matrot(phi,the,psi)
%=========================================================================
%     ROUTINE CHE CALCOLA LA MATRICE DI ROTAZIONE DATI GLI ANGOLI
%     DI EULERO NOTAZIONE GOLDSTEIN
%
%     {X}=[R]{X'}
%=========================================================================
      r(1,1)=cos(psi)*cos(phi)-cos(the)*sin(phi)*sin(psi);
      r(1,2)=-sin(psi)*cos(phi)-cos(the)*sin(phi)*cos(psi);
      r(1,3)=sin(the)*sin(phi);
%r
      r(2,1)=cos(psi)*sin(phi)+cos(the)*cos(phi)*sin(psi);
      r(2,2)=-sin(psi)*sin(phi)+cos(the)*cos(phi)*cos(psi);
      r(2,3)=-sin(the)*cos(phi);
%r
      r(3,1)=sin(the)*sin(psi);
      r(3,2)=sin(the)*cos(psi);
      r(3,3)=cos(the);
end

end   %Method

end   %Classdef