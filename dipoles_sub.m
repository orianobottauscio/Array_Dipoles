classdef dipoles_sub

methods (Static)

%--------------------------------------------------------------------
% Functions for dipole treatment
%--------------------------------------------------------------------
function [MAGNETI_STATO] = assegna_valori_dipoles(MAGNETI_GEO,MATERIALI,MAGNETI_STATO,Jresidua_relativo)
MAGNETI_STATO.Jresidua_magneti=zeros(MAGNETI_GEO.Ndipoli_tot,1);
MAGNETI_STATO.Hc_magneti=zeros(MAGNETI_GEO.Ndipoli_tot,1);
MAGNETI_STATO.MuR=zeros(MAGNETI_GEO.Ndipoli_tot,1);
MAGNETI_STATO.JrNL=zeros(MAGNETI_GEO.Ndipoli_tot,1);
MAGNETI_STATO.HcNL=zeros(MAGNETI_GEO.Ndipoli_tot,1);

if MATERIALI.JHcurve
  for nd=1:MAGNETI_GEO.Ndipoli_tot
    nmat=MAGNETI_GEO.mat_number(nd,1);
    MAGNETI_STATO.MuR(nd,1)=MATERIALI.MuR(nmat,1);
    MAGNETI_STATO.JrNL(nd,1)=MATERIALI.JrNL(nmat,1);
    MAGNETI_STATO.HcNL(nd,1)=MATERIALI.HcNL(nmat,1);
    Jresidua_ref=MATERIALI.JrNL(nmat,1)*Jresidua_relativo(nmat,1);
    Hc_ref=MATERIALI.HcNL(nmat,1)*Jresidua_relativo(nmat,1);
    Tref=MATERIALI.Tref(nmat,1);
    JvT=MATERIALI.JvT(nmat,1);
    MAGNETI_STATO.Jresidua_magneti(nd,1) = util.Jr_versus_T(JvT,MAGNETI_STATO.Temperature_magneti(nd,1),Tref,Jresidua_ref);  
    HvT=MATERIALI.HvT(nmat,1:2);
    MAGNETI_STATO.Hc_magneti(nd,1) = util.Hc_versus_T(HvT,MAGNETI_STATO.Temperature_magneti(nd,1),Tref,Hc_ref); 
  end
else
  for nd=1:MAGNETI_GEO.Ndipoli_tot
    nmat=MAGNETI_GEO.mat_number(nd,1);
    MAGNETI_STATO.MuR(nd,1)=MATERIALI.MuR(nmat,1);
    MAGNETI_STATO.Jresidua_magneti(nd,1)=MATERIALI.JrNL(nmat,1)*Jresidua_relativo(nmat,1);
  end
end
return
end



function [L,U,P,Tnoto2,BFIELD,MAGNETI_STATO] = output_dipoles(first,reaction,NL,recompute,MAGNETI_GEO,MAGNETI_STATO,...
          MATERIALI,POINTS,L,U,P,Tnoto2,Mu0,output_JHmag)
if (first || recompute) && reaction
    [L,U,P,Tnoto2] = dipoles_sub.sistema(MAGNETI_GEO.Ndipoli_tot,MAGNETI_GEO.cubo_dimU2,MAGNETI_GEO.cubo_dimV2,MAGNETI_GEO.cubo_dimW2,...
                     Mu0,MAGNETI_STATO.MuR,MAGNETI_GEO.angle,MAGNETI_GEO.xdip,MAGNETI_GEO.ydip,MAGNETI_GEO.zdip);
end
if reaction
    [J,H] = dipoles_sub.compute_J_dipoles(reaction,NL,MAGNETI_GEO.mat_number,MATERIALI.HJcurve_all,MATERIALI.HJcurve_all_points,MAGNETI_GEO.Ndipoli_tot,...
        MAGNETI_STATO.Jresidua_magneti,MAGNETI_STATO.Hc_magneti,MAGNETI_STATO.JrNL,MAGNETI_STATO.HcNL,Mu0,MAGNETI_STATO.MuR,L,U,P,Tnoto2);
else
    [J,H] = dipoles_sub.compute_J_dipoles(reaction,NL,MAGNETI_GEO.mat_number,MATERIALI.HJcurve_all,MATERIALI.HJcurve_all_points,MAGNETI_GEO.Ndipoli_tot,...
        MAGNETI_STATO.Jresidua_magneti,MAGNETI_STATO.Hc_magneti,MAGNETI_STATO.JrNL,MAGNETI_STATO.HcNL,Mu0,MAGNETI_STATO.MuR,[],[],[],[]);
end
[BFIELD.Bx5,BFIELD.By5,BFIELD.Bz5] = dipoles_sub.compute_field(POINTS.Npoint,POINTS.xp,POINTS.yp,POINTS.zp,...
        MAGNETI_GEO.Ndipoli_tot,MAGNETI_GEO.xdip,MAGNETI_GEO.ydip,MAGNETI_GEO.zdip,MAGNETI_GEO.angle,...
        MAGNETI_GEO.cubo_dimU2,MAGNETI_GEO.cubo_dimV2,MAGNETI_GEO.cubo_dimW2,J,Mu0);

if output_JHmag
  MAGNETI_STATO.Jmag=J;
  MAGNETI_STATO.Hmag=H;
end

return
end

function [J,H] = compute_J_dipoles(reaction,NL,mat_number,HJcurve_all,HJcurve_all_points,Ndipoli_tot,...
                 Jresidua_magneti,Hc_magneti,JrNL,HcNL,Mu0,MuR,L,U,P,Tnoto2)
errorFL_lim=1e-7;
iterFP_max=500;
writeFP=false;
if reaction
  if NL
    Jaus=zeros(Ndipoli_tot,1);
    residual=zeros(Ndipoli_tot,1);
    for iterFP=1:iterFP_max
      aa=Tnoto2(:,1:Ndipoli_tot).*(Jresidua_magneti+residual);
      Tnoto3=sum(aa,2);
      y=L\(P*Tnoto3);
      H=U\y;       %Valori di H per ciascun dipolo
      clear y
%.........Adattamento caratteristica NL
      for np=1:Ndipoli_tot
        clear HJJ
        HJJ(:,1)=HJcurve_all(1:HJcurve_all_points(mat_number(np)),1,mat_number(np))/HcNL(np,1)*Hc_magneti(np,1);    %Aggiusta caratt. NL sulla base del valore di Hc effettivo
        HJJ(:,2)=HJcurve_all(1:HJcurve_all_points(mat_number(np)),2,mat_number(np))/JrNL(np,1)*Jresidua_magneti(np,1);    %Aggiusta caratt. NL sulla base del valore di Jr effettivo
        Jaus(np,1)=interp1(HJJ(:,1),HJJ(:,2),H(np,1),'linear','extrap');
      end
      residual=Jaus-Jresidua_magneti-Mu0.*MuR.*H;
      if iterFP > 1
          CL2=sqrt(sum(abs(H.^2-Hold.^2)));
          CMM=sqrt(sum(H.^2));
          CNO=CL2/CMM;
          if writeFP
             fprintf('IterFP: %d Norma L2: %e Valor Medio: %e Errore: %e\n',iterFP,CL2/Ndipoli_tot,CMM/Ndipoli_tot,CNO);
          end
          if CNO < errorFL_lim
              break
          end
      end
      Hold=H;
    end
    J=Jresidua_magneti+Mu0.*MuR.*H+residual;   %Valori di J per ciascun dipolo
  else
    aa=Tnoto2(:,1:Ndipoli_tot).*Jresidua_magneti;
    Tnoto3=sum(aa,2);
    y=L\(P*Tnoto3);
    H=U\y;       %Valori di H per ciascun dipolo
    clear y
    J=Jresidua_magneti+Mu0.*MuR.*H;   %Valori di J per ciascun dipolo
  end
else
  J=Jresidua_magneti;
  H=(J-Jresidua_magneti)./(Mu0.*MuR);
end

return
end


function [L,U,P,Tnoto2] = sistema(Ndipoli_tot,cubo_dimU2,cubo_dimV2,cubo_dimW2,...
                                  Mu0,MuR,angle,xdip,ydip,zdip)
%----------------------------------------------------------------------------------
% Costruzione sistema algebrico per caso reaction
%----------------------------------------------------------------------------------
  parallelo = 1;
  Matrix=zeros(Ndipoli_tot,Ndipoli_tot);
  Tnoto2=zeros(Ndipoli_tot,Ndipoli_tot);
  if parallelo == 1
    Nsn=(cubo_dimU2(1:Ndipoli_tot,1).*cubo_dimW2(1:Ndipoli_tot,1))./...
        (cubo_dimU2(1:Ndipoli_tot,1).*cubo_dimV2(1:Ndipoli_tot,1)+...
         cubo_dimV2(1:Ndipoli_tot,1).*cubo_dimW2(1:Ndipoli_tot,1)+...
         cubo_dimW2(1:Ndipoli_tot,1).*cubo_dimU2(1:Ndipoli_tot,1));
    for nd=1:Ndipoli_tot
      Matrix(nd,nd)=-Nsn(nd)*Mu0*MuR(nd,1)-Mu0;
    end
%
    Un(:,1)=cos(pi./180.*angle(1:Ndipoli_tot,1));
    Un(:,2)=sin(pi./180.*angle(1:Ndipoli_tot,1));
    Un(:,3)=0;
    for md=1:Ndipoli_tot
      xp=xdip(1:Ndipoli_tot,1);   %Baricentro dipolo (ricevente contributo)
      yp=ydip(1:Ndipoli_tot,1);
      zp=zdip(1:Ndipoli_tot,1);
%  
      xp(md,1)=1000.;   %Valore fittizio per evitare singolarità di md su nd
      yp(md,1)=1000.;
      zp(md,1)=1000.;
%  
      U2m(1)=sin(pi/180*angle(md,1));
      U2m(2)=-cos(pi/180*angle(md,1));
      U2m(3)=0;
      [gx,gy,gz]=dipoles_sub.geometrico_plus2_par(xdip(md,1),ydip(md,1),zdip(md,1),U2m,...
         cubo_dimU2(md,1),cubo_dimV2(md,1),cubo_dimW2(md,1),xp,yp,zp);
      gx(md,1)=0.0;
      gy(md,1)=0.0;
      gz(md,1)=0.0;
%      
      g=gx.*Un(:,1)+gy.*Un(:,2);
      Tnoto2(1:Ndipoli_tot,md)=-Mu0.*g(1:Ndipoli_tot);
      Matrix(1:Ndipoli_tot,md)=Matrix(1:Ndipoli_tot,md)+Mu0.*g(1:Ndipoli_tot).*Mu0.*MuR(1:Ndipoli_tot);
    end
    for nd=1:Ndipoli_tot
      Tnoto2(nd,nd)=Nsn(nd);
    end
  else
    for nd=1:Ndipoli_tot
      xp=xdip(nd,1);   %Baricentro dipolo (ricevente contributo)
      yp=ydip(nd,1);
      zp=zdip(nd,1);
      Un(1)=cos(pi/180*angle(nd,1));
      Un(2)=sin(pi/180*angle(nd,1));
      Un(3)=0;
      Nsn=(cubo_dimU2(nd,1)*cubo_dimW2(nd,1))/...
          (cubo_dimU2(nd,1)*cubo_dimV2(nd,1)+...
           cubo_dimV2(nd,1)*cubo_dimW2(nd,1)+...
           cubo_dimW2(nd,1)*cubo_dimU2(nd,1));
      Tnoto2(nd,nd)=Nsn;
      Matrix(nd,nd)=-Nsn*Mu0*MuR(nd,1)-Mu0;
      for md=1:Ndipoli_tot
        if md ~= nd
          U2m(1)=sin(pi/180*angle(md,1));
          U2m(2)=-cos(pi/180*angle(md,1));
          U2m(3)=0;
          [gx,gy,~]=dipoles_sub.geometrico_plus2(xdip(md,1),ydip(md,1),zdip(md,1),U2m,...
             cubo_dimU2(md,1),cubo_dimV2(md,1),cubo_dimW2(md,1),xp,yp,zp);
          g=gx*Un(1)+gy*Un(2);
          Tnoto2(nd,md)=-Mu0*g;
          Matrix(nd,md)=Matrix(nd,md)+Mu0*g*Mu0*MuR(md,1);
        end
      end
    end
  end
  [L,U,P] = lu(Matrix);
return
end

%--------------------------------------------------------------------
% Dipole Plus sequenziale (Hegel-Herbert)
%--------------------------------------------------------------------
function [gx,gy,gz]=geometrico_plus2_par(xdip,ydip,zdip,U2,dimU,dimV,dimW,xp,yp,zp)
% Versione parallela con set di punti di calcolo in input
%  Sistema magnete rispetto a master
  drif(1,1)=xdip;
  drif(1,2)=ydip;
  drif(1,3)=zdip;
  drif(2,1)=xdip;
  drif(2,2)=ydip;
  drif(2,3)=zdip+dimW/2;
  drif(3,1)=xdip+U2(1);
  drif(3,2)=ydip+U2(2);
  drif(3,3)=zdip;

%coord. punto di calcolo nel sistema magnete
  [uP,vP,wP] = util.global_to_local_par(xp,yp,zp,drif);

%Componenti geometriche in sistema locale (dipolo sorgente)
  dd=length(xp);
  gu=zeros(dd,1);
  gv=zeros(dd,1);
  gw=zeros(dd,1);
  for k=1:2
    for l=1:2
      for m=1:2
        a=sqrt((uP+(-1).^k*dimU/2).^2+(vP+(-1).^l*dimV/2).^2+(wP+(-1).^m.*dimW/2).^2);
        gu=gu+(-1).^(k+l+m).*log(wP+(-1).^m.*dimW/2+a);
        gw=gw+(-1).^(k+l+m).*log(uP+(-1).^k.*dimU/2+a);
        a1=vP+(-1).^l.*dimV/2;
        a2=uP+(-1).^k.*dimU/2;
        b1=abs(uP+(-1).^k.*dimU/2).*(wP+(-1).^m.*dimW/2);
        b2=abs(vP+(-1).^l.*dimV/2).*a;
        a3=atan(b1./b2);
        gv=gv-(-1).^(k+l+m).*(a1.*a2)./(abs(a1).*abs(a2)).*a3;
      end
    end
  end
  coef=1.0/(4*pi);
  Mu0=4e-7*pi;
  gu=gu.*coef/Mu0;
  gv=gv.*coef/Mu0;
  gw=gw.*coef/Mu0;
  [gx(:,1),gy(:,1),gz(:,1)] = util.vector_local_to_global_par(gu,gv,gw,drif);
end

function [gx,gy,gz]=geometrico_plus2(xdip,ydip,zdip,U2,dimU,dimV,dimW,xp,yp,zp)
%  Sistema magnete rispetto a master
  drif(1,1)=xdip;
  drif(1,2)=ydip;
  drif(1,3)=zdip;
  drif(2,1)=xdip;
  drif(2,2)=ydip;
  drif(2,3)=zdip+dimW/2;
  drif(3,1)=xdip+U2(1);
  drif(3,2)=ydip+U2(2);
  drif(3,3)=zdip;

%coord. punto di calcolo nel sistema magnete
  [uP,vP,wP] = util.global_to_local(xp,yp,zp,drif);

%Componenti B in sistema locale (dipolo sorgente)
  gu=0.0;
  gv=0.0;
  gw=0.0;
  for k=1:2
    for l=1:2
      for m=1:2
        a=sqrt((uP+(-1)^k*dimU/2)^2+(vP+(-1)^l*dimV/2)^2+(wP+(-1)^m*dimW/2)^2);
        gu=gu+(-1)^(k+l+m)*log(wP+(-1)^m*dimW/2+a);
        gw=gw+(-1)^(k+l+m)*log(uP+(-1)^k*dimU/2+a);
        a1=vP+(-1)^l*dimV/2;
        a2=uP+(-1)^k*dimU/2;
        b1=abs(uP+(-1)^k*dimU/2)*(wP+(-1)^m*dimW/2);
        b2=abs(vP+(-1)^l*dimV/2)*a;
        a3=atan(b1/b2);
        gv=gv-(-1)^(k+l+m)*(a1*a2)/(abs(a1)*abs(a2))*a3;
      end
    end
  end
  coef=1.0/(4*pi);
  Mu0=4e-7*pi;
  gu=gu*coef/Mu0;
  gv=gv*coef/Mu0;
  gw=gw*coef/Mu0;

  [gx,gy,gz] = util.vector_local_to_global(gu,gv,gw,drif);

end


function [Bx5,By5,Bz5] = compute_field(Npoint,xp,yp,zp,Ndipoli_tot,xdip,ydip,zdip,angle,cubo_dimU2,cubo_dimV2,cubo_dimW2,J,Mu0)
  Hx5=zeros(Npoint,1);
  Hy5=zeros(Npoint,1);
  Hz5=zeros(Npoint,1);
  for nd=1:Ndipoli_tot
    U2(1)=sin(pi/180*angle(nd,1));
    U2(2)=-cos(pi/180*angle(nd,1));
    U2(3)=0;
    [Hx5,Hy5,Hz5]=dipoles_sub.field_plus2_par(Hx5,Hy5,Hz5,xdip(nd,1),ydip(nd,1),zdip(nd,1),U2,...
         cubo_dimU2(nd,1),cubo_dimV2(nd,1),cubo_dimW2(nd,1),J(nd,1),xp,yp,zp);
  end
  Bx5=Mu0*Hx5;
  By5=Mu0*Hy5;
  Bz5=Mu0*Hz5;
return
end

%--------------------------------------------------------------------
% Dipole Plus parallelo (Hegel-Herbert)
%--------------------------------------------------------------------
function [Hx5,Hy5,Hz5]=field_plus2_par(Hx5,Hy5,Hz5,xdip,ydip,zdip,U2,dimU,dimV,dimW,valore_Jr,xp,yp,zp)
%  Sistema magnete rispetto a master
  drif(1,1)=xdip;
  drif(1,2)=ydip;
  drif(1,3)=zdip;
  drif(2,1)=xdip;
  drif(2,2)=ydip;
  drif(2,3)=zdip+dimW/2;
  drif(3,1)=xdip+U2(1);
  drif(3,2)=ydip+U2(2);
  drif(3,3)=zdip;

%coord. punto di calcolo nel sistema magnete
  [uP,vP,wP] = util.global_to_local_par(xp,yp,zp,drif);

%Componenti B in sistema locale
  dd=length(xp);
  Hu=zeros(dd,1);
  Hv=zeros(dd,1);
  Hw=zeros(dd,1);
  for k=1:2
    for l=1:2
      for m=1:2
        a=sqrt((uP+(-1).^k*dimU/2).^2+(vP+(-1).^l*dimV/2).^2+(wP+(-1).^m.*dimW/2).^2);
        Hu=Hu+(-1).^(k+l+m).*log(wP+(-1).^m.*dimW/2+a);
        Hw=Hw+(-1).^(k+l+m).*log(uP+(-1).^k.*dimU/2+a);
        a1=vP+(-1).^l.*dimV/2;
        a2=uP+(-1).^k.*dimU/2;
        b1=abs(uP+(-1).^k.*dimU/2).*(wP+(-1).^m.*dimW/2);
        b2=abs(vP+(-1).^l.*dimV/2).*a;
        a3=atan(b1./b2);
        Hv=Hv-(-1).^(k+l+m).*(a1.*a2)./(abs(a1).*abs(a2)).*a3;
      end
    end
  end
  coef=valore_Jr./(4*pi);
  Mu0=4e-7*pi;
  Hu=Hu.*coef/Mu0;
  Hv=Hv.*coef/Mu0;
  Hw=Hw.*coef/Mu0;
  [Hx(:,1),Hy(:,1),Hz(:,1)] = util.vector_local_to_global_par(Hu,Hv,Hw,drif);
  Hx5(:,1)=Hx5(:,1)+Hx(:,1);
  Hy5(:,1)=Hy5(:,1)+Hy(:,1);
  Hz5(:,1)=Hz5(:,1)+Hz(:,1);
end


function [Hx,Hy,Hz]=field_plus2_serial(xdip,ydip,zdip,U2,dimU,dimV,dimW,valore_Jr,xp,yp,zp)
% Calcolo contributo di un solo magnete in un solo punto
%  Sistema magnete rispetto a master
  drif(1,1)=xdip;
  drif(1,2)=ydip;
  drif(1,3)=zdip;
  drif(2,1)=xdip;
  drif(2,2)=ydip;
  drif(2,3)=zdip+dimW/2;
  drif(3,1)=xdip+U2(1);
  drif(3,2)=ydip+U2(2);
  drif(3,3)=zdip;

%coord. punto di calcolo nel sistema magnete
  [uP,vP,wP] = util.global_to_local_par(xp,yp,zp,drif);

%Componenti B in sistema locale
  dd=length(xp);
  Hu=zeros(dd,1);
  Hv=zeros(dd,1);
  Hw=zeros(dd,1);
  for k=1:2
    for l=1:2
      for m=1:2
        a=sqrt((uP+(-1).^k*dimU/2).^2+(vP+(-1).^l*dimV/2).^2+(wP+(-1).^m.*dimW/2).^2);
        Hu=Hu+(-1).^(k+l+m).*log(wP+(-1).^m.*dimW/2+a);
        Hw=Hw+(-1).^(k+l+m).*log(uP+(-1).^k.*dimU/2+a);
        a1=vP+(-1).^l.*dimV/2;
        a2=uP+(-1).^k.*dimU/2;
        b1=abs(uP+(-1).^k.*dimU/2).*(wP+(-1).^m.*dimW/2);
        b2=abs(vP+(-1).^l.*dimV/2).*a;
        a3=atan(b1./b2);
        Hv=Hv-(-1).^(k+l+m).*(a1.*a2)./(abs(a1).*abs(a2)).*a3;
      end
    end
  end
  coef=valore_Jr./(4*pi);
  Mu0=4e-7*pi;
  Hu=Hu.*coef/Mu0;
  Hv=Hv.*coef/Mu0;
  Hw=Hw.*coef/Mu0;
  [Hx,Hy,Hz] = util.vector_local_to_global_par(Hu,Hv,Hw,drif);
end


function  [TORQUE] = compute_torque(flag_all,TORQUE,MAGNETI_GEO,MAGNETI_STATO,Mu0)
TORQUE.momento_total=size(MAGNETI_GEO.Ndipoli_tot,1);
Dipole_group_ini=zeros(MAGNETI_GEO.Nlist,1);
Dipole_group_fin=zeros(MAGNETI_GEO.Nlist,1);
fin=0;
for Nr=1:MAGNETI_GEO.Nlist
  Dipole_group_ini(Nr,1)=fin+1;
  for nd=fin+1:MAGNETI_GEO.Ndipoli_tot
    if MAGNETI_GEO.group(nd,1)>Nr
      fin=nd-1;
      Dipole_group_fin(Nr,1)=fin;
      break
    end
  end
end
Dipole_group_fin(MAGNETI_GEO.Nlist,1)=MAGNETI_GEO.Ndipoli_tot;
for Nr=1:MAGNETI_GEO.Nlist
  i1=Dipole_group_ini(Nr,1);
  i2=Dipole_group_fin(Nr,1);
  nnn=i2-i1+1;
  Hx5=zeros(nnn,1);
  Hy5=zeros(nnn,1);
  Hz5=zeros(nnn,1);
  pos_angolare=zeros(nnn,1);  
  ii=0;
  for np=i1:i2
    ii=ii+1;
    xp=MAGNETI_GEO.xdip(np,1);
    yp=MAGNETI_GEO.ydip(np,1);
    zp=MAGNETI_GEO.zdip(np,1);
    aa=atan2(yp,xp);
    if aa<0
        aa=aa+2*pi;
    end
    pos_angolare(ii,1)=aa/pi*180;
    if flag_all
       nd_ini=1;
       nd_fin=MAGNETI_GEO.Ndipoli_tot;
    else
       nd_ini=i1;
       nd_fin=i2;
    end
    for nd=nd_ini:nd_fin
      if nd ~= np  
        U2(1)=sin(pi/180*MAGNETI_GEO.angle(nd,1));
        U2(2)=-cos(pi/180*MAGNETI_GEO.angle(nd,1));
        U2(3)=0;
        [Hx,Hy,Hz]=dipoles_sub.field_plus2_serial(MAGNETI_GEO.xdip(nd,1),MAGNETI_GEO.ydip(nd,1),MAGNETI_GEO.zdip(nd,1),U2,...
             MAGNETI_GEO.cubo_dimU2(nd,1),MAGNETI_GEO.cubo_dimV2(nd,1),MAGNETI_GEO.cubo_dimW2(nd,1),MAGNETI_STATO.Jmag(nd,1),xp,yp,zp);
        Hx5(ii,1)=Hx5(ii,1)+Hx;
        Hy5(ii,1)=Hy5(ii,1)+Hy;
        Hz5(ii,1)=Hz5(ii,1)+Hz;
      end
    end
  end
  clear BB
  BB(:,1)=Mu0*Hx5;
  BB(:,2)=Mu0*Hy5;
  BB(:,3)=Mu0*Hz5;
  clear U1
  U1(:,1)=cos(pi/180*MAGNETI_GEO.angle(i1:i2,1));
  U1(:,2)=sin(pi/180*MAGNETI_GEO.angle(i1:i2,1));
  U1(:,3)=zeros(nnn,1);
  clear JJ
  JJ(:,1)=MAGNETI_STATO.Jmag(i1:i2,1).*U1(:,1);
  JJ(:,2)=MAGNETI_STATO.Jmag(i1:i2,1).*U1(:,2);
  JJ(:,3)=MAGNETI_STATO.Jmag(i1:i2,1).*U1(:,3);
  clear tt
  tt=cross(JJ,BB);
  TORQUE.momento_total(i1:i2,1)=tt(:,3).*MAGNETI_GEO.cubo_dimU2(i1:i2,1).*MAGNETI_GEO.cubo_dimV2(i1:i2,1).*MAGNETI_GEO.cubo_dimW2(i1:i2,1)/Mu0;
end
return
end


end   %Method

end   %Classdef