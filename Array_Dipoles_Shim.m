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
global Neval_function
global Neval_function_write
global Method_Actual
global PM_position
global PM_presence_save


Code_Version='2.1.1';

fprintf('Code: Array_Dipoles_Shim, Version: %s\n',Code_Version);
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
[ierr,POINTS,typemis] = input_DSVMeas(toml_data,POINTS);
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
[ierr,SHIMMING_OPT,SHIMMING_GEO,MAGNETI_GEO_PREVIOUS] = input_shimming_opt(toml_data,SHIMMING_OPT,SHIMMING_GEO,POINTS,scale,input_data);
if ierr==1
    return
end

SHIMMING_OPT.Mu0=Mu0;


if strcmpi(typemis,'V')   %Volume data in DSV
  [Bmin,Bmax,Bmean,DEV1_ppm,DEV2_ppm] = util.variability(SHIMMING_OPT.B5);
  fprintf('Volume measurements in DSV - \n');
  fprintf('Bmin,Bmax,Bmean [mT]: %f %f %f\n',Bmin*1000,Bmax*1000,Bmean*1000);
  fprintf('DEV1 [ppm]: %f\n',DEV1_ppm);
  fprintf('DEV2 [ppm]: %f\n',DEV2_ppm);
  if strcmpi(SHIMMING_OPT.DevOpt,'DEV3')
      fprintf('Opt. DEV3 not compatible with volume measurements\n');
      return
  end
elseif strcmpi(typemis,'S')   %Surface data in DSV
  [Bmin,Bmax,B0,DEV3_ppm] = util.variability3(SHIMMING_OPT.B5);
  fprintf('Surface measurements in DSV - \n');
  fprintf('Bmin,Bmax,B0 [mT]: %f %f %f\n',Bmin*1000,Bmax*1000,B0*1000);
  fprintf('DEV3 [ppm]: %f\n',DEV3_ppm);
  if strcmpi(SHIMMING_OPT.DevOpt,'DEV1') || strcmpi(SHIMMING_OPT.DevOpt,'DEV2')
      fprintf('Opt. DEV1/DEV2 not compatible with surface measurements\n');
      return
  end
end

%Parametro usato per ottimizzazione
if strcmpi(SHIMMING_OPT.DevOpt,'DEV1')
  SHIMMING_OPT.original_deviation=DEV1_ppm;
elseif strcmpi(SHIMMING_OPT.DevOpt,'DEV2')
  SHIMMING_OPT.original_deviation=DEV2_ppm;
elseif strcmpi(SHIMMING_OPT.DevOpt,'DEV3')
  SHIMMING_OPT.original_deviation=DEV3_ppm;
elseif strcmpi(SHIMMING_OPT.DevOpt,'dB')
  SHIMMING_OPT.original_deviation=(Bmax-Bmin)*1e6;
end  


% Define the positions on the PMs for the optimization run
[SHIMMING_OPT] = define_PM_positions(SHIMMING_OPT,SHIMMING_GEO);

% One block of PMs with PM rotation (0:360). If the angle is > 360 the PM
% is not present
nvars=SHIMMING_OPT.Ndipoli_tot_shim_max;
x_input_Min(1,1:nvars)=SHIMMING_OPT.PM_AngleVar(1);
x_input_Max(1,1:nvars)=SHIMMING_OPT.PM_AngleVar(2);

% Tipo di funzione obiettivo
if SIMUL_DATA.reaction || SIMUL_DATA.NL
   fprintf('Optimization with non-ideal magnets\n');
   ideal_dipole=0;
   Jresidua_magneti=0;   %Non usato
else
   fprintf('Optimization with ideal magnets\n');
   ideal_dipole=1;
   Temperature_magneti=TEMPERATURE.Tactual;
   nmat=SHIMMING_OPT.PM_material;
   Jresidua_ref=MATERIALI.JrNL(nmat,1);
   Tref=MATERIALI.Tref(nmat,1);
   JvT=MATERIALI.JvT(nmat,1);
   Jresidua_magneti = util.Jr_versus_T(JvT,Temperature_magneti,Tref,Jresidua_ref);  
end

f = @(x_input) deviation_in_sphere3(x_input,MAGNETI_GEO,TEMPERATURE,MAGNETI_STATO,MATERIALI,SIMUL_DATA,...
                                            SHIMMING_OPT,SHIMMING_GEO,POINTS,ideal_dipole,Jresidua_magneti);

%---------------------------------------------------------------
%  Section of Global Algorithm optimization
%---------------------------------------------------------------
SHIMMING_OPT.GA=false;
SHIMMING_OPT.SW=false;
SHIMMING_OPT.SA=false;
SHIMMING_OPT.WO=false;
PM_position=true;   %N. PM variabile
Tolerance=SHIMMING_OPT.Tolerance;
Neval_function=0;
Neval_function_write=true;
if strcmpi(SHIMMING_OPT.Global_Algorithm,'GA')
  MaxGenerations1=SHIMMING_OPT.MaxGenerations1;
  MaxStallGenerations=SHIMMING_OPT.MaxStallGenerations;
  PopulationSize=SHIMMING_OPT.PopulationSize;
  Rate=SHIMMING_OPT.MutationRate;
  Method_Actual='(GA)';
  start_optimization = tic;
  fprintf('----OPTIMIZATION with GA----\n');
  fprintf('N. of parameters: %d\n',nvars);
  MaxGenerations=MaxGenerations1*nvars;
  if Rate == 0
    options = optimoptions('ga','MaxGenerations',MaxGenerations,...
                                'MaxStallGenerations',MaxStallGenerations,...
                                'PopulationSize',PopulationSize,...
                                'FunctionTolerance',Tolerance);
  else
    options = optimoptions('ga','MaxGenerations',MaxGenerations,...
                                'MaxStallGenerations',MaxStallGenerations,...
                                'PopulationSize',PopulationSize,...
                                'FunctionTolerance',Tolerance,...
                                'MutationFcn', {@mutationuniform, Rate});
  end
  [x_input,evaluation_optGA,exitflagGA,OutputGA] = ga(f,nvars,[],[],[],[],x_input_Min,x_input_Max,[],[],options);
  x_input_good=x_input;    %Save the output of GA
  computational_time_GA=toc(start_optimization);
  evaluation_opt=evaluation_optGA;
  fprintf('Number of generations: %d\n', OutputGA.generations);
  fprintf('Computational time: %.1f\n',computational_time_GA);
  fprintf('Optimal value of Deviation after GA: %.1f\n',evaluation_optGA);
  SHIMMING_OPT.GA=true;
elseif strcmpi(SHIMMING_OPT.Global_Algorithm,'SW')
  Method_Actual='(SW)';
  start_optimization = tic;
  fprintf('----OPTIMIZATION with ParticleSwarm----\n');
  fprintf('N. of parameters: %d\n',nvars);
  options = optimoptions('particleswarm','FunctionTolerance',Tolerance);
%  options = optimoptions('particleswarm','FunctionTolerance',Tolerance,'HybridFcn',@fmincon);
  [x_input,~,exitflagSW,OutputSW] = particleswarm(f,nvars,x_input_Min,x_input_Max,options);
  evaluation_optSW = deviation_in_sphere3(x_input,MAGNETI_GEO,TEMPERATURE,MAGNETI_STATO,MATERIALI,SIMUL_DATA,...
                                                  SHIMMING_OPT,SHIMMING_GEO,POINTS,ideal_dipole,Jresidua_magneti);   %Effettua una rivalutazione della funzione obiettivo per casi nonlinear
  x_input_good=x_input;    %Save the output of SW
  computational_time_GA=toc(start_optimization);
  evaluation_opt=evaluation_optSW;
  fprintf('Computational time: %.1f\n',computational_time_GA);
  fprintf('Optimal value of Deviation after SW: %.1f\n',evaluation_optSW);
  SHIMMING_OPT.SW=true;
elseif strcmpi(SHIMMING_OPT.Global_Algorithm,'SA')
  Method_Actual='(SA)';
  start_optimization = tic;
  fprintf('----OPTIMIZATION with Simulated Annealing----\n');
  fprintf('N. of parameters: %d\n',nvars);
%  options = optimoptions('simulannealbnd','FunctionTolerance',Tolerance,'HybridFcn',@fmincon);
  options = optimoptions('simulannealbnd','FunctionTolerance',Tolerance);
  [x_input,~,exitflagSA,OutputSA] = simulannealbnd(f,x_input_start,x_input_Min,x_input_Max,options);
  evaluation_optSA = deviation_in_sphere3(x_input,MAGNETI_GEO,TEMPERATURE,MAGNETI_STATO,MATERIALI,SIMUL_DATA,...
                                                  SHIMMING_OPT,SHIMMING_GEO,POINTS,ideal_dipole,Jresidua_magneti);   %Effettua una rivalutazione della funzione obiettivo per casi nonlinear
  x_input_good=x_input;    %Save the output of SA
  computational_time_GA=toc(start_optimization);
  evaluation_opt=evaluation_optSA;
  fprintf('Computational time: %.1f\n',computational_time_GA);
  fprintf('Optimal value of Deviation after SA: %.1f\n',evaluation_optSA);
  SHIMMING_OPT.SA=true;
elseif strcmpi(SHIMMING_OPT.Global_Algorithm,'WO')
  Method_Actual='(WO)';
  start_optimization = tic;
  fprintf('----OPTIMIZATION with Whale Optimization ----\n');
  fprintf('N. of parameters: %d\n',nvars);
  SearchAgents_no=30; % Number of search agents
  Max_iteration=500; % Maximum numbef of iterations
  [evaluation_optWO,x_input_good,WOA_cg_curve]=woa.WOA(SearchAgents_no,Max_iteration,x_input_Min,x_input_Max,nvars,f);
  computational_time_GA=toc(start_optimization);
  evaluation_opt=evaluation_optWO;
  fprintf('Computational time: %.1f\n',computational_time_GA);
  fprintf('Optimal value of Deviation after WO: %.1f\n',evaluation_optWO);
  SHIMMING_OPT.WO=true;
elseif strcmpi(SHIMMING_OPT.Global_Algorithm,'Null')
  computational_time_GA=0;
  x_input_good=(x_input_Min+x_input_Max)/2;   %Serve per input metodi successivi (se non viene usato GA)
  evaluation_opt = f(x_input_good);
  fprintf('Starting Deviation value: %.1f\n',evaluation_opt);
end    

%---------------------------------------------------------------
% Congela le posizioni dei PMs che portano a valore ottimale
% della Global Optimization (PM_presence_save)
% Riduce il numero di variabili in input (No. PMs fissato)
%---------------------------------------------------------------
if ~strcmpi(SHIMMING_OPT.Global_Algorithm,'Null')
  [~,PM_presence_save] = from_input_to_PM(x_input_good);
  if sum(PM_presence_save)>0
    PM_position=false;   %N. PM fissato
    x_input_good_rid=x_input_good(PM_presence_save>0);
    nvars_rid=sum(PM_presence_save);
    fprintf('No. Added magnets after Global Optimization (kept fixed): %d\n',nvars_rid);
  else
    x_input_good_rid=x_input_good;
    nvars_rid=nvars;
    fprintf('No added magnets after Global Optimization\n');
  end
else
  x_input_good_rid=x_input_good;
  nvars_rid=nvars;
end
%---------------------------------------------------------------
%  Section of Local Algorithm optimization
%  In questa fase, non modifica più il numero di PMs usati
%  che sono quelli stabiliti dalla variabile 'PM_presence_save'
%---------------------------------------------------------------
SHIMMING_OPT.FM=false;
SHIMMING_OPT.PS=false;
Var_input=SHIMMING_OPT.Local_Var_input;  %Variazioni angolo per Local Optimization

if strcmpi(SHIMMING_OPT.Local_Algorithm,'PS')
  Method_Actual='(PS)';
  Neval_function=0;
  Neval_function_write=true;
  start_optimization = tic;
  fprintf('----OPTIMIZATION with Patternsearch----\n');
  fprintf('N. of parameters: %d\n',nvars_rid);
  options = optimoptions('patternsearch','FunctionTolerance',Tolerance,'StepTolerance',Tolerance);
%  options = optimoptions('patternsearch','FunctionTolerance',Tolerance,'StepTolerance',Tolerance,Algorithm="nups");
  x_input_start = x_input_good_rid;
  [x_input_rid,~,exitflagPS,OutputPS] = patternsearch(f,x_input_start,[],[],[],[],...
                                    x_input_start-Var_input,x_input_start+Var_input,[],options);
  evaluation_optPS = deviation_in_sphere3(x_input_rid,MAGNETI_GEO,TEMPERATURE,MAGNETI_STATO,MATERIALI,SIMUL_DATA,...
                                                  SHIMMING_OPT,SHIMMING_GEO,POINTS,ideal_dipole,Jresidua_magneti);   %Effettua una rivalutazione della funzione obiettivo per casi nonlinear
  computational_time_LA=toc(start_optimization);
  fprintf('Computational time: %.1f\n',computational_time_LA);
  fprintf('Optimal value of Deviation after PS: %.1f\n',evaluation_optPS);
  if evaluation_optPS > evaluation_opt   %Se al termine di Patternsearch risultato peggiore, ripristina risultato GA
    fprintf('WARNING: Patternsearch provides worst result: recovered GA solution\n');
  else
    evaluation_opt=evaluation_optPS;
    x_input_good_rid = x_input_rid;
  end
  SHIMMING_OPT.PS=true;
elseif strcmpi(SHIMMING_OPT.Local_Algorithm,'FM')
  Method_Actual='(FM)';
  Neval_function=0;
  Neval_function_write=true;
  start_optimization = tic;
  fprintf('----OPTIMIZATION with Fmincon----\n');
  fprintf('N. of parameters: %d\n',nvars_rid);
  options = optimoptions('fmincon','ObjectiveLimit',Tolerance,'Algorithm','active-set','MaxFunctionEvaluations',5000);
%  options = optimoptions('fmincon','ObjectiveLimit',Tolerance,'MaxFunctionEvaluations',5000);
  x_input_start = x_input_good_rid;
  [x_input_rid,~,exitflagFM,OutputFM] = fmincon(f,x_input_start,[],[],[],[],...
      x_input_start-Var_input,x_input_start+Var_input,[],options);
  evaluation_optFM = deviation_in_sphere3(x_input_rid,MAGNETI_GEO,TEMPERATURE,MAGNETI_STATO,MATERIALI,SIMUL_DATA,...
                                                  SHIMMING_OPT,SHIMMING_GEO,POINTS,ideal_dipole,Jresidua_magneti);   %Effettua una rivalutazione della funzione obiettivo per casi nonlinear
  computational_time_LA=toc(start_optimization);
  fprintf('Computational time: %.1f\n',computational_time_LA);
  fprintf('Optimal value of Deviation after FM: %.1f\n',evaluation_optFM);
  if evaluation_optFM > evaluation_opt   %Se al termine di Fmincon risultato peggiore, ripristina il precedente
    fprintf('WARNING: Fmincon provides worst result: recovered best solution\n');
  else
    evaluation_opt=evaluation_optFM;
    x_input_good_rid = x_input_rid;
  end
  SHIMMING_OPT.FM=true;
elseif contains(SHIMMING_OPT.Local_Algorithm,'+MC')
  Max_iter=200;
  Max_iter_noChange=5;
  Max_fval=700;
  N_estrazioni=SHIMMING_OPT.MC_extraction;
  Neval_function=0;
  start_optimization = tic;
  x_input_start = x_input_good_rid;
  if strcmpi(SHIMMING_OPT.Local_Algorithm,'LSQ+MC')
    Method_Actual='(LSQ+MC)';
    Neval_function_write=false;
    fprintf('----OPTIMIZATION with LSQ+MC----\n');
    fprintf('N. of parameters: %d\n',nvars_rid);
    options = optimoptions('lsqnonlin','Algorithm','trust-region-reflective','FunctionTolerance',...
        Tolerance,'MaxFunctionEvaluations',Max_fval,'display','off','FinDiffRelStep',1e-4);
    iter_noChange=0;
    iter=0;
    while iter_noChange<=Max_iter_noChange && iter <= Max_iter
      iter=iter+1;
      [x_input_rid,~,exitflagLS,~] = lsqnonlin(f,x_input_start,x_input_start-Var_input,x_input_start+Var_input,...
                                 [],[],[],[],[],options);
      [x_output,evaluation_optMC] = monte_carlo.MC_opt(f,x_input_rid,N_estrazioni,Var_input,nvars_rid);
      if evaluation_optMC < evaluation_opt
        evaluation_opt = evaluation_optMC;
        x_input_good_rid=x_output;
        iter_noChange=0;
      else
        iter_noChange=iter_noChange+1;
      end
      x_input_start = x_output;
      fprintf('Iter No. %d - Deviation value: %.1f\n',iter,evaluation_opt);
    end
    computational_time_LA=toc(start_optimization);
    fprintf('Computational time: %.1f\n',computational_time_LA);
    fprintf('Optimal value of Deviation after LSQ+MC: %.1f\n',evaluation_opt);
  elseif strcmpi(SHIMMING_OPT.Local_Algorithm,'PS+MC')
    Method_Actual='(PS+MC)';
    Neval_function_write=false;
    fprintf('----OPTIMIZATION with PS+MC----\n');
    fprintf('N. of parameters: %d\n',nvars_rid);
    options = optimoptions('patternsearch','FunctionTolerance',Tolerance,'StepTolerance',Tolerance);
    iter_noChange=0;
    iter=0;
    while iter_noChange<=Max_iter_noChange && iter <= Max_iter
      iter=iter+1;
      [x_input_rid,~,exitflagPS,~] = patternsearch(f,x_input_start,[],[],[],[],...
                                 x_input_start-Var_input,x_input_start+Var_input,[],options);
      [x_output,evaluation_optMC] = monte_carlo.MC_opt(f,x_input_rid,N_estrazioni,Var_input,nvars_rid);
      if evaluation_optMC < evaluation_opt
        evaluation_opt = evaluation_optMC;
        x_input_good_rid=x_output;
        iter_noChange=0;
      else
        iter_noChange=iter_noChange+1;
      end
      x_input_start = x_output;
      fprintf('Iter No. %d - Deviation value: %.1f\n',iter,evaluation_opt);
    end
    computational_time_LA=toc(start_optimization);
    fprintf('Computational time: %.1f\n',computational_time_LA);
    fprintf('Optimal value of Deviation after PS+MC: %.1f\n',evaluation_opt);
  elseif strcmpi(SHIMMING_OPT.Local_Algorithm,'FM+MC')
    Method_Actual='(FM+MC)';
    Neval_function_write=false;
    fprintf('----OPTIMIZATION with FM+MC----\n');
    fprintf('N. of parameters: %d\n',nvars_rid);
    options = optimoptions('fmincon','ObjectiveLimit',Tolerance,'MaxFunctionEvaluations',5000);
    iter_noChange=0;
    iter=0;
    while iter_noChange<=Max_iter_noChange && iter <= Max_iter
      iter=iter+1;
      [x_input_rid,~,exitflagFM,~] = fmincon(f,x_input_start,[],[],[],[],...
                                 x_input_start-Var_input,x_input_start+Var_input,[],options);
      [x_output,evaluation_optMC] = monte_carlo.MC_opt(f,x_input_rid,N_estrazioni,Var_input,nvars_rid);
      if evaluation_optMC < evaluation_opt
        evaluation_opt = evaluation_optMC;
        x_input_good_rid=x_output;
        iter_noChange=0;
      else
        iter_noChange=iter_noChange+1;
      end
      x_input_start = x_output;
      fprintf('Iter No. %d - Deviation value: %.1f\n',iter,evaluation_opt);
    end
    computational_time_LA=toc(start_optimization);
    fprintf('Computational time: %.1f\n',computational_time_LA);
    fprintf('Optimal value of Deviation after FM+MC: %.1f\n',evaluation_opt);
  end
elseif strcmpi(SHIMMING_OPT.Local_Algorithm,'MC')
  N_estrazioni=SHIMMING_OPT.MC_extraction;
  Method_Actual='(MC)';
  Neval_function=0;
  Neval_function_write=false;
  start_optimization = tic;
  fprintf('----OPTIMIZATION with MC----\n');
  fprintf('N. of parameters: %d\n',nvars_rid);
  x_input_start = x_input_good_rid;

  for iter=1:4
    [x_output,evaluation_optMC] = monte_carlo.MC_opt(f,x_input_start,N_estrazioni,Var_input,nvars_rid);
    if evaluation_optMC < evaluation_opt
        evaluation_opt = evaluation_optMC;
        x_input_good_rid=x_output;
        x_input_start = x_output;
    end
    fprintf('Iter No. %d - Deviation value: %.1f\n',iter,evaluation_opt);
    Var_input=Var_input/2;
  end
  computational_time_LA=toc(start_optimization);
  fprintf('Computational time: %.1f\n',computational_time_LA);
  fprintf('Optimal value of Deviation after MC: %.1f\n',evaluation_opt);
end

%=======================================================================
computational_time=computational_time_GA+computational_time_LA;
fprintf('Total computational time: %.1f\n',computational_time);
%
if strcmpi(SHIMMING_OPT.Global_Algorithm,'GA')
  fprintf('--------------------------------------------------------\n');
  fprintf('Genetic Algorithm\n');
  fprintf('Exit flag: %d\n',exitflagGA);
  fprintf('Optimal value of Deviation after GA: %.2f\n',evaluation_optGA);
  fprintf('--------------------------------------------------------\n');
elseif strcmpi(SHIMMING_OPT.Global_Algorithm,'PS')
  fprintf('--------------------------------------------------------\n');
  fprintf('ParticleSwarm Algorithm\n');
  fprintf('Exit flag: %d\n',exitflagSW);
  fprintf('Optimal value of Deviation after SW: %.2f\n',evaluation_optSW);
  fprintf('--------------------------------------------------------\n');
elseif strcmpi(SHIMMING_OPT.Global_Algorithm,'SA')
  fprintf('--------------------------------------------------------\n');
  fprintf('Simulated Annealing Algorithm\n');
  fprintf('Exit flag: %d\n',exitflagSA);
  fprintf('Optimal value of Deviation after SA: %.2f\n',evaluation_optSA);
  fprintf('--------------------------------------------------------\n');
elseif strcmpi(SHIMMING_OPT.Global_Algorithm,'WO')
  fprintf('--------------------------------------------------------\n');
  fprintf('Whale Optimization Algorithm\n');
  fprintf('Optimal value of Deviation after SA: %.2f\n',evaluation_optWO);
  fprintf('--------------------------------------------------------\n');
end

%=======================================================================
if PM_position
  [dir_mag_shim,PM_presence] = from_input_to_PM(x_input_good_rid);
else
  [dir_mag_shim,PM_presence] = from_input_to_PM2(x_input_good_rid,SHIMMING_OPT,PM_presence_save);
end

% Reconstruct position of PM added in shimming and adjust final field in DSV
[SHIMMING_OPT,MAGNETI_GEO,SHIMMING_GEO] = configuration_at_end(PM_presence,dir_mag_shim,...
                                          MAGNETI_GEO,TEMPERATURE,MAGNETI_STATO,...
                                          MATERIALI,SIMUL_DATA,SHIMMING_OPT,SHIMMING_GEO,POINTS);

fprintf('#######################################################\n');
fprintf('SUMMARY\n');

if MAGNETI_GEO.Ndipoli_tot == 0
  fprintf('NO MAGNETS ADDED AT THIS RUN\n');
  fprintf('Original value deviation is kept\n');
  return
end

if strcmpi(SHIMMING_OPT.DevOpt,'DEV1')
  fprintf('Original value of DEV1 before shimming (ppm): %.2f\n',SHIMMING_OPT.original_deviation);
  [Bmin,Bmax,Bmean,DEV1_ppm,DEV2_ppm] = util.variability(SHIMMING_OPT.B5);
  fprintf('Bmin,Bmax,Bmean [mT]: %.2f %.2f %.2f\n',Bmin*1000,Bmax*1000,Bmean*1000);
  fprintf('dB [uT]: %.3f\n',(Bmax-Bmin)*1e6);
  fprintf('DEV1 [ppm]: %.1f\n',DEV1_ppm);
  fprintf('DEV2 [ppm]: %.1f\n',DEV2_ppm);
elseif strcmpi(SHIMMING_OPT.DevOpt,'DEV2')
  fprintf('Original value of DEV2 before shimming (ppm): %.2f\n',SHIMMING_OPT.original_deviation);
  [Bmin,Bmax,Bmean,DEV1_ppm,DEV2_ppm] = util.variability(SHIMMING_OPT.B5);
  fprintf('Bmin,Bmax,Bmean [mT]: %.2f %.2f %.2f\n',Bmin*1000,Bmax*1000,Bmean*1000);
  fprintf('dB [uT]: %.3f\n',(Bmax-Bmin)*1e6);
  fprintf('DEV1 [ppm]: %.1f\n',DEV1_ppm);
  fprintf('DEV2 [ppm]: %.1f\n',DEV2_ppm);
elseif strcmpi(SHIMMING_OPT.DevOpt,'DEV3')
  fprintf('Original value of DEV3 before shimming (ppm): %.2f\n',SHIMMING_OPT.original_deviation);
  [Bmin,Bmax,B0,DEV3_ppm] = util.variability3(SHIMMING_OPT.B5);
  fprintf('Bmin,Bmax,B0 [mT]: %.2f %.2f %.2f\n',Bmin*1000,Bmax*1000,B0*1000);
  fprintf('dB [uT]: %.3f\n',(Bmax-Bmin)*1e6);
  fprintf('DEV3 [ppm]: %.1f\n',DEV3_ppm);
elseif strcmpi(SHIMMING_OPT.DevOpt,'dB')
  fprintf('Original value of dB before shimming (uT): %.2f\n',SHIMMING_OPT.original_deviation);
  if strcmpi(typemis,'V')   %Volume data in DSV
    [Bmin,Bmax,Bmean,DEV1_ppm,DEV2_ppm] = util.variability(SHIMMING_OPT.B5);
    fprintf('Bmin,Bmax,Bmean [mT]: %.2f %.2f %.2f\n',Bmin*1000,Bmax*1000,Bmean*1000);
    fprintf('dB [uT]: %.3f\n',(Bmax-Bmin)*1e6);
    fprintf('DEV1 [ppm]: %.1f\n',DEV1_ppm);
    fprintf('DEV2 [ppm]: %.1f\n',DEV2_ppm);
  elseif strcmpi(typemis,'S')   %Surface data in DSV
    [Bmin,Bmax,B0,DEV3_ppm] = util.variability3(SHIMMING_OPT.B5);
    fprintf('Bmin,Bmax,B0 [mT]: %.2f %.2f %.2f\n',Bmin*1000,Bmax*1000,B0*1000);
    fprintf('dB [uT]: %.3f\n',(Bmax-Bmin)*1e6);
    fprintf('DEV3 [ppm]: %.1f\n',DEV3_ppm);
  end
end

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
save(filemat,'computational_time','-append');
if SHIMMING_OPT.GA
  save(filemat,'OutputGA','evaluation_optGA','-append');
end
if SHIMMING_OPT.SW
  save(filemat,'OutputSW','evaluation_optSW','-append');
end
if SHIMMING_OPT.SA
  save(filemat,'OutputSA','evaluation_optSA','-append');
end
if SHIMMING_OPT.FM
  save(filemat,'OutputFM','-append');
end
if SHIMMING_OPT.PS
  save(filemat,'OutputPS','-append');
end
if SHIMMING_OPT.WO
  save(filemat,'WOA_cg_curve','evaluation_optWO','-append');
end


% Extract general info for optimization to be saved
if SHIMMING_GEO.sector_type == 1
  INFO_SHIM_GEO = rmfield (SHIMMING_GEO,{'MatrixShimming','Xslot','Yslot'});
elseif SHIMMING_GEO.sector_type == 2
  INFO_SHIM_GEO = rmfield (SHIMMING_GEO,{'MatrixShimming','angoloslot'});
end

INFO_SHIM_OPT = rmfield (SHIMMING_OPT,{'B5'});
save(filemat,'INFO_SHIM_GEO','INFO_SHIM_OPT','-append');
save(filemat,'SIMUL_DATA','TEMPERATURE','-append');
%
if strcmpi(SHIMMING_OPT.DevOpt,'DEV1') || strcmpi(SHIMMING_OPT.DevOpt,'DEV2')
  Bmean_mT=Bmean*1000;
  save(filemat,'DEV1_ppm','DEV2_ppm','Bmean_mT','-append');
elseif strcmpi(SHIMMING_OPT.DevOpt,'DEV3')
  B0_mT=B0*1000;
  save(filemat,'DEV3_ppm','B0_mT','-append');
elseif strcmpi(SHIMMING_OPT.DevOpt,'dB')
  if strcmpi(typemis,'V')
    Bmean_mT=Bmean*1000;
    save(filemat,'DEV1_ppm','DEV2_ppm','Bmean_mT','-append');
  elseif strcmpi(typemis,'S')
    B0_mT=B0*1000;
    save(filemat,'DEV3_ppm','B0_mT','-append');
  end
end
save(filemat,'MATERIALI','POINTS','-append');   %Servono per poter usare il programma Check_Shim_configuration

fprintf('Results of shimming saved in file: %s\n',filemat);
fprintf('#######################################################\n');

if evaluation_opt > SHIMMING_OPT.original_deviation
  fprintf('WARNING: the starting value of Deviation before this optimization step was better\n');
  fprintf('DELETE this run\n');
  fprintf('Original value of Deviation before shimming: %.2f\n',SHIMMING_OPT.original_deviation);
  fprintf('Deviation after this step: %.1f\n',evaluation_opt);
end

end


function [dir_mag_shim,PM_presence] = from_input_to_PM(x_input)
PM_presence=ones(length(x_input),1);
PM_presence(x_input>360)=0;
dir_mag_shim=transpose(x_input);
return
end

function [dir_mag_shim,PM_presence] = from_input_to_PM2(x_input_rid,SHIMMING_OPT,PM_presence_save)
PM_presence=PM_presence_save;
dir_mag_shim=zeros(SHIMMING_OPT.Ndipoli_tot_shim_max,1);
dir_mag_shim(PM_presence>0)=x_input_rid;
return
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


function [SHIMMING_OPT,MAGNETI_GEO,SHIMMING_GEO] = configuration_at_end(PM_presence,dir_mag_shim,...
                                                   MAGNETI_GEO,TEMPERATURE,MAGNETI_STATO,...
                                                   MATERIALI,SIMUL_DATA,SHIMMING_OPT,SHIMMING_GEO,POINTS)
[MAGNETI_GEO,SHIMMING_GEO] = build_dipole_shim_end(PM_presence,dir_mag_shim,MAGNETI_GEO,SHIMMING_OPT,SHIMMING_GEO);
if MAGNETI_GEO.Ndipoli_tot == 0
  return
end
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


function [ierr,POINTS,typemis] = input_DSVMeas(toml_data,POINTS)
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
if ~isfield(toml_data.shimming.measurement,'typemis')
    fprintf('Field [shimming.pm.measurement.typemis] not present\n');
    ierr=1;
    return
end
typemis=toml_data.shimming.measurement.typemis;
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


function [ierr,SHIMMING_OPT,SHIMMING_GEO,MAGNETI_GEO_PREVIOUS] = input_shimming_opt(toml_data,SHIMMING_OPT,SHIMMING_GEO,POINTS,scale,input_data)
ierr=0;
if ~isfield(toml_data.shimming,'opt_rules')
    fprintf('Field [shimming.opt_rules] not present\n');
    ierr=1;
    return
end
if ~isfield(toml_data.shimming.opt_rules,'PM_AngleVar')
    fprintf('Field [shimming.opt_rules.PM_AngleVar] not present\n');
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

if isfield(toml_data.shimming.opt_rules,'ring_to_be_used')
  SHIMMING_OPT.ring_to_be_used=toml_data.shimming.opt_rules.ring_to_be_used;
  aa=SHIMMING_OPT.ring_to_be_used > SHIMMING_GEO.Nring_tot;
  if any(aa)
    fprintf('Error in data of Field [shimming.opt_rules.ring_to_be_used]\n');
    ierr=1;
    return
  end
else
  SHIMMING_OPT.ring_to_be_used=1:1:SHIMMING_GEO.Nring_tot;
end
SHIMMING_OPT.Nring_used=length(SHIMMING_OPT.ring_to_be_used);

if isfield(toml_data.shimming.opt_rules,'sector_to_be_used')
  SHIMMING_OPT.sector_to_be_used=toml_data.shimming.opt_rules.sector_to_be_used;
  aa=SHIMMING_OPT.sector_to_be_used > SHIMMING_GEO.Nsector_available;
  if any(aa)
    fprintf('Error in data of Field [shimming.opt_rules.sector_to_be_used]\n');
    ierr=1;
    return
  end
else
  SHIMMING_OPT.sector_to_be_used=1:1:SHIMMING_GEO.Nsector_available;
end
SHIMMING_OPT.Nsector_used=length(SHIMMING_OPT.sector_to_be_used);

if isfield(toml_data.shimming.opt_rules,'slot_to_be_used')
  SHIMMING_OPT.slot_to_be_used=toml_data.shimming.opt_rules.slot_to_be_used;
  aa=SHIMMING_OPT.slot_to_be_used > SHIMMING_GEO.Nslot_per_sector;
  if any(aa)
    fprintf('Error in data of Field [shimming.opt_rules.Nslot_per_sector]\n');
    ierr=1;
    return
  end
else
  SHIMMING_OPT.slot_to_be_used=1:1:SHIMMING_GEO.Nslot_per_sector;
end
SHIMMING_OPT.Nslot_used=length(SHIMMING_OPT.slot_to_be_used);


SHIMMING_OPT.PM_AngleVar=toml_data.shimming.opt_rules.PM_AngleVar;

SHIMMING_OPT.cubo_dimU2_shim_ref=toml_data.shimming.opt_rules.PM_size(1)*scale;
SHIMMING_OPT.cubo_dimV2_shim_ref=toml_data.shimming.opt_rules.PM_size(2)*scale;
SHIMMING_OPT.cubo_dimW2_shim_ref=toml_data.shimming.opt_rules.PM_size(3)*scale;
SHIMMING_OPT.PM_material=toml_data.shimming.opt_rules.PM_material;


if isfield(toml_data.shimming.opt_rules,'dev_to_be_optimized')
  SHIMMING_OPT.DevOpt=toml_data.shimming.opt_rules.dev_to_be_optimized;
else
  SHIMMING_OPT.DevOpt='dB';
end

if isfield(toml_data.shimming.opt_rules,'shimming_Dir')
  SHIMMING_OPT.ShimmingDir=toml_data.shimming.opt_rules.shimming_Dir;
else
  SHIMMING_OPT.ShimmingDir=erase(input_data,'.toml');
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
if isfield(toml_data.shimming.opt_rules,'PopulationSize')
  SHIMMING_OPT.PopulationSize=toml_data.shimming.opt_rules.PopulationSize;
else
  SHIMMING_OPT.PopulationSize=200;
end
if isfield(toml_data.shimming.opt_rules,'MutationRate')
  SHIMMING_OPT.MutationRate=toml_data.shimming.opt_rules.MutationRate;   %Mutation function uniform for rate dato
else
  SHIMMING_OPT.MutationRate=0';  %Mutation function di default (gaussian)
end

if isfield(toml_data.shimming.opt_rules,'Tolerance')
  SHIMMING_OPT.Tolerance=toml_data.shimming.opt_rules.Tolerance;
else
  SHIMMING_OPT.Tolerance=1e-4;
end

% Select among the available Global Algorithms: GA, SW, SA
if isfield(toml_data.shimming.opt_rules,'Global_Algorithm')
  SHIMMING_OPT.Global_Algorithm=toml_data.shimming.opt_rules.Global_Algorithm;
else
  SHIMMING_OPT.Global_Algorithm='SW';   %Default: SW
end

% Select among the available Local Algorithms: Patternsearch,Fmincon,LSQ+MC 
if isfield(toml_data.shimming.opt_rules,'Local_Algorithm')
  SHIMMING_OPT.Local_Algorithm=toml_data.shimming.opt_rules.Local_Algorithm;
else
  SHIMMING_OPT.Local_Algorithm='LSQ+MC';  %Default: LSQ+MC
end

if isfield(toml_data.shimming.opt_rules,'Local_Var_input')
  SHIMMING_OPT.Local_Var_input=toml_data.shimming.opt_rules.Local_Var_input;
else
  SHIMMING_OPT.Local_Var_input=1;  %Default: 1 deg
end

if isfield(toml_data.shimming.opt_rules,'MC_extraction')
  SHIMMING_OPT.MC_extraction=toml_data.shimming.opt_rules.MC_extraction;
else
  SHIMMING_OPT.MC_extraction=5000;  %Default: 5000
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
if ~isfield(toml_data.shimming.pm,'radius')
    fprintf('Field [shimming.pm.radius] not present\n');
    ierr=1;
    return
end
if ~isfield(toml_data.shimming.pm,'Nsector_available')
    fprintf('Field [shimming.pm.Nsector_available] not present\n');
    ierr=1;
    return
end

if ~isfield(toml_data.shimming.pm,'sector_type')
    fprintf('Field [shimming.pm.sector_type] not present\n');
    ierr=1;
    return
end

if ~isfield(toml_data.shimming.pm,'Nslot_per_sector')
    fprintf('Field [shimming.pm.Nslot_per_sector] not present\n');
    ierr=1;
    return
end

if ~isfield(toml_data.shimming.pm,'ring_Zposition')
    fprintf('Field [shimming.pm.ring_Zposition] not present\n');
    ierr=1;
    return
end
SHIMMING_GEO.Raggio_shim=toml_data.shimming.pm.radius*scale;
SHIMMING_GEO.Nsector_available=toml_data.shimming.pm.Nsector_available;
SHIMMING_GEO.Nslot_per_sector=toml_data.shimming.pm.Nslot_per_sector;
SHIMMING_GEO.sector_type=toml_data.shimming.pm.sector_type;

if SHIMMING_GEO.sector_type == 1
% settore rettangolare
  if ~isfield(toml_data.shimming.pm,'sector_width')
      fprintf('Field [shimming.pm.sector_width] not present\n');
      ierr=1;
      return
  end
  SHIMMING_GEO.Ampiezza_settore=toml_data.shimming.pm.sector_width*scale;
elseif SHIMMING_GEO.sector_type == 2
% settore ad arco
  if ~isfield(toml_data.shimming.pm,'sector_angular_width')
      fprintf('Field [shimming.pm.sector_angular_width] not present\n');
      ierr=1;
      return
  end
  SHIMMING_GEO.Ampiezza_angolare_settore=toml_data.shimming.pm.sector_angular_width;
end

SHIMMING_GEO.ring_Zposition=toml_data.shimming.pm.ring_Zposition*scale;
SHIMMING_GEO.Nring_tot=length(SHIMMING_GEO.ring_Zposition);

Dang_shim=360./SHIMMING_GEO.Nsector_available;
aus=0:1:SHIMMING_GEO.Nsector_available-1;
Pos_Angular_shim=Dang_shim*aus;
% Structure which keep trace of position filled with magnets (the value corresponds to the Run number)
SHIMMING_GEO.MatrixShimming=zeros(SHIMMING_GEO.Nring_tot,SHIMMING_GEO.Nsector_available,SHIMMING_GEO.Nslot_per_sector);  %Matrice che conterrà posizioni occupate
%
if SHIMMING_GEO.sector_type == 1
  DV_shim=SHIMMING_GEO.Ampiezza_settore/SHIMMING_GEO.Nslot_per_sector;   %Distanza di separazione tra magneti del settore (v locale)
  SHIMMING_GEO.Xslot=zeros(SHIMMING_GEO.Nsector_available,SHIMMING_GEO.Nslot_per_sector);  %position X of each slot of the sectors
  SHIMMING_GEO.Yslot=zeros(SHIMMING_GEO.Nsector_available,SHIMMING_GEO.Nslot_per_sector);  %position Y of each slot of the sectors
  for Ns=1:SHIMMING_GEO.Nsector_available
    ang=Pos_Angular_shim(Ns)/180*pi;
    xb_sector=SHIMMING_GEO.Raggio_shim*cos(ang);
    yb_sector=SHIMMING_GEO.Raggio_shim*sin(ang);
    drif(1,1)=xb_sector;
    drif(1,2)=yb_sector;
    drif(1,3)=0.0;
    drif(2,1)=xb_sector;
    drif(2,2)=yb_sector;
    drif(2,3)=1.0;
    x2_sector=1.2*SHIMMING_GEO.Raggio_shim*cos(ang);
    y2_sector=1.2*SHIMMING_GEO.Raggio_shim*sin(ang);
    drif(3,1)=x2_sector;
    drif(3,2)=y2_sector;
    drif(3,3)=0.0;
%
    for Nm=1:SHIMMING_GEO.Nslot_per_sector
      mm=mod(SHIMMING_GEO.Nslot_per_sector,2);
      if mm==0
% pari
        Nsub=SHIMMING_GEO.Nslot_per_sector/2+0.5;
      else
% dispari
        Nsub=fix(SHIMMING_GEO.Nslot_per_sector/2)+1;
      end
      uloc=0.0;
      vloc=(Nm-Nsub)*DV_shim;
      wloc=0.0;
      [xxs,yys,~] = util.local_to_global(uloc,vloc,wloc,drif);
      SHIMMING_GEO.Xslot(Ns,Nm)=xxs;
      SHIMMING_GEO.Yslot(Ns,Nm)=yys;
    end
  end
elseif SHIMMING_GEO.sector_type == 2
  Dalfa_shim=SHIMMING_GEO.Ampiezza_angolare_settore/SHIMMING_GEO.Nslot_per_sector;   %Angolo in gradi di separazione tra magneti del settore
  SHIMMING_GEO.angoloslot=zeros(SHIMMING_GEO.Nsector_available,SHIMMING_GEO.Nslot_per_sector);  %Angular position of each slot of the sectors
  for Ns=1:SHIMMING_GEO.Nsector_available
    for Nm=1:SHIMMING_GEO.Nslot_per_sector
      mm=mod(SHIMMING_GEO.Nslot_per_sector,2);
      if mm==0
        Nsub=SHIMMING_GEO.Nslot_per_sector/2+0.5;
        SHIMMING_GEO.angoloslot(Ns,Nm)=Pos_Angular_shim(Ns)+(Nm-Nsub)*Dalfa_shim;
      else
        Nsub=fix(SHIMMING_GEO.Nslot_per_sector/2)+1;
        SHIMMING_GEO.angoloslot(Ns,Nm)=Pos_Angular_shim(Ns)+(Nm-Nsub)*Dalfa_shim;
      end
    end
  end
end
return
end



function [SHIMMING_OPT] = define_PM_positions(SHIMMING_OPT,SHIMMING_GEO)
%----------------------------------------------------------------------
% Define PM positions in the considered optimization run
%----------------------------------------------------------------------
ii=0;
for N1=1:SHIMMING_OPT.Nring_used
  Nr=SHIMMING_OPT.ring_to_be_used(N1);   %Number of used Ring
  for N2=1:SHIMMING_OPT.Nsector_used
    Ns=SHIMMING_OPT.sector_to_be_used(N2);
    for N3=1:SHIMMING_OPT.Nslot_used
      Nm=SHIMMING_OPT.slot_to_be_used(N3);
      if SHIMMING_GEO.MatrixShimming(Nr,Ns,Nm) == 0
        ii=ii+1;
        if SHIMMING_GEO.sector_type == 1
% Settori rettangolari
          SHIMMING_OPT.PM_posX(ii,1)=SHIMMING_GEO.Xslot(Ns,Nm);
          SHIMMING_OPT.PM_posY(ii,1)=SHIMMING_GEO.Yslot(Ns,Nm);
          SHIMMING_OPT.PM_posZ(ii,1)=SHIMMING_GEO.ring_Zposition(Nr);
          SHIMMING_OPT.RingIndex(ii,1)=Nr;
          SHIMMING_OPT.SectorIndex(ii,1)=Ns;
          SHIMMING_OPT.PositionIndex(ii,1)=Nm;
        elseif SHIMMING_GEO.sector_type == 2
% Settori ad arco
          angolo=SHIMMING_GEO.angoloslot(Ns,Nm);
          SHIMMING_OPT.PM_posX(ii,1)=SHIMMING_GEO.Raggio_shim*cos(angolo/180*pi);
          SHIMMING_OPT.PM_posY(ii,1)=SHIMMING_GEO.Raggio_shim*sin(angolo/180*pi);
          SHIMMING_OPT.PM_posZ(ii,1)=SHIMMING_GEO.ring_Zposition(Nr);
          SHIMMING_OPT.RingIndex(ii,1)=Nr;
          SHIMMING_OPT.SectorIndex(ii,1)=Ns;
          SHIMMING_OPT.PositionIndex(ii,1)=Nm;
        end
      end
    end
  end
end
SHIMMING_OPT.Ndipoli_tot_shim_max=ii;   %Numero massimo di PM che possono essere aggiunti con questo Run di ottimizzazione
return
end




function [MAGNETI_GEO] = build_dipole_shim(PM_presence,dir_mag_shim,MAGNETI_GEO,SHIMMING_OPT,SHIMMING_GEO)
%----------------------------------------------------------------------
% Input dati magneti shimming
%----------------------------------------------------------------------
Nmax=SHIMMING_OPT.Ndipoli_tot_shim_max;
MAGNETI_GEO.xdip=zeros(Nmax,1);
MAGNETI_GEO.ydip=zeros(Nmax,1);
MAGNETI_GEO.zdip=zeros(Nmax,1);
MAGNETI_GEO.angle=zeros(Nmax,1);
MAGNETI_GEO.cubo_dimU2=zeros(Nmax,1);
MAGNETI_GEO.cubo_dimV2=zeros(Nmax,1);
MAGNETI_GEO.cubo_dimW2=zeros(Nmax,1);
MAGNETI_GEO.mat_number=zeros(Nmax,1);
MAGNETI_GEO.Ndipoli_tot=0;
%
for ii=1:Nmax
  if PM_presence(ii) > 0
    MAGNETI_GEO.Ndipoli_tot=MAGNETI_GEO.Ndipoli_tot+1;
    MAGNETI_GEO.xdip(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.PM_posX(ii,1);
    MAGNETI_GEO.ydip(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.PM_posY(ii,1);
    MAGNETI_GEO.zdip(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.PM_posZ(ii,1);
    MAGNETI_GEO.angle(MAGNETI_GEO.Ndipoli_tot,1)=dir_mag_shim(ii,1);
    MAGNETI_GEO.cubo_dimU2(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.cubo_dimU2_shim_ref;
    MAGNETI_GEO.cubo_dimV2(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.cubo_dimV2_shim_ref;
    MAGNETI_GEO.cubo_dimW2(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.cubo_dimW2_shim_ref;
    MAGNETI_GEO.mat_number(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.PM_material;
  end
end
return
end


function [MAGNETI_GEO,SHIMMING_GEO] = build_dipole_shim_end(PM_presence,dir_mag_shim,MAGNETI_GEO,SHIMMING_OPT,SHIMMING_GEO)
MAGNETI_GEO.Ndipoli_tot=0;
for ii=1:SHIMMING_OPT.Ndipoli_tot_shim_max
  if PM_presence(ii) > 0
    Nr=SHIMMING_OPT.RingIndex(ii,1);
    Ns=SHIMMING_OPT.SectorIndex(ii,1);
    Nm=SHIMMING_OPT.PositionIndex(ii,1);
    if SHIMMING_GEO.MatrixShimming(Nr,Ns,Nm) == 0
      MAGNETI_GEO.Ndipoli_tot=MAGNETI_GEO.Ndipoli_tot+1;
      MAGNETI_GEO.xdip(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.PM_posX(ii,1);
      MAGNETI_GEO.ydip(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.PM_posY(ii,1);
      MAGNETI_GEO.zdip(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.PM_posZ(ii,1);
      MAGNETI_GEO.angle(MAGNETI_GEO.Ndipoli_tot,1)=dir_mag_shim(ii,1);
      MAGNETI_GEO.cubo_dimU2(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.cubo_dimU2_shim_ref;
      MAGNETI_GEO.cubo_dimV2(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.cubo_dimV2_shim_ref;
      MAGNETI_GEO.cubo_dimW2(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.cubo_dimW2_shim_ref;
      MAGNETI_GEO.mat_number(MAGNETI_GEO.Ndipoli_tot,1)=SHIMMING_OPT.PM_material;
      SHIMMING_GEO.MatrixShimming(Nr,Ns,Nm)=SHIMMING_OPT.Run;
    end
  end
end
return
end



function deviation = deviation_in_sphere3(x_input,MAGNETI_GEO,TEMPERATURE,MAGNETI_STATO,MATERIALI,SIMUL_DATA,...
                                                  SHIMMING_OPT,SHIMMING_GEO,POINTS,ideal_dipole,Jresidua_magneti)
global Neval_function
global Neval_function_write
global Method_Actual
global PM_position
global PM_presence_save

if PM_position
  [dir_mag_shim,PM_presence] = from_input_to_PM(x_input);
else
  [dir_mag_shim,PM_presence] = from_input_to_PM2(x_input,SHIMMING_OPT,PM_presence_save);
end

%------------------------------------------------------
%  Input geometrico configurazione shimming
%------------------------------------------------------
[MAGNETI_GEO] = build_dipole_shim(PM_presence,dir_mag_shim,MAGNETI_GEO,SHIMMING_OPT,SHIMMING_GEO);

if MAGNETI_GEO.Ndipoli_tot > 0
%-------------------------------------------------------------------
%  Assegnazione parametri magnetici per PM di configurazione shimming
%-------------------------------------------------------------------
  if ideal_dipole == 1
    [BFIELD.Bx5,BFIELD.By5,BFIELD.Bz5]=dipoles_sub.ideal_dipole(POINTS,MAGNETI_GEO,Jresidua_magneti);
%    [BFIELD.Bx5,BFIELD.By5,BFIELD.Bz5] = dipoles_sub.compute_field(POINTS.Npoint,POINTS.xp,POINTS.yp,POINTS.zp,...
%        MAGNETI_GEO.Ndipoli_tot,MAGNETI_GEO.xdip,MAGNETI_GEO.ydip,MAGNETI_GEO.zdip,MAGNETI_GEO.angle,...
%        MAGNETI_GEO.cubo_dimU2,MAGNETI_GEO.cubo_dimV2,MAGNETI_GEO.cubo_dimW2,J,SHIMMING_OPT.Mu0);
    B5 = util.estrai_output(SIMUL_DATA,BFIELD);
  else    
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
  end
else
  B5=zeros(length(SHIMMING_OPT.B5),1);
end    
B5_tot=SHIMMING_OPT.B5+B5;

if strcmpi(SHIMMING_OPT.DevOpt,'DEV1')
  [~,~,~,deviation,~] = util.variability(B5_tot);
elseif strcmpi(SHIMMING_OPT.DevOpt,'DEV2')
  [~,~,~,~,deviation] = util.variability(B5_tot);
elseif strcmpi(SHIMMING_OPT.DevOpt,'DEV3')
  [~,~,~,deviation] = util.variability3(B5_tot);
elseif strcmpi(SHIMMING_OPT.DevOpt,'dB')
  [~,~,deviation] = util.variability_delta(B5_tot);
end

Neval_function=Neval_function+1;
if mod(Neval_function,100)==0 && Neval_function_write
    fprintf('%s Fval: %d Added magnets: %d - Deviation: %.1f\n',Method_Actual,...
    Neval_function,MAGNETI_GEO.Ndipoli_tot,deviation);
end
return
end