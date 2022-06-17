clear all
close all
clc


solverType = {'adj'};
solverType = {'constDirect','sqp','interior-point'};
solverType = {'constDirect','interior-point'};
solverType = {'constDirect'};
ObjectiveType = {'TotalSignal'}
ObjectiveType = {'TotalSignal','SumQuad'}
gpList = [3, 4,5]
gpList = [5]
uncertainList = [3]
snrList = [2,10]
snrList = [2,10,25]
snrList = [2,10,20,25]
snrList = [2,5,10,20]
snrList = [2,5,10,15,20]
% pareto trade off total signal vs MI

myFAList =  repmat([6:6:35],2,1);
myFAList(:,end)  =  35;
myFAList(2,:) =  28;

myFAList =  repmat([3 35],2,1);
myFAList(2,:) =  28;



%numsolves = numel(ObjectiveType)* numel(solverType) * length(gpList) * length(uncertainList) * length(snrList) + (2+size(myFAList,2))*length(snrList) + 2*length(snrList)
numsolves =  (2+size(myFAList,2)+2)*length(snrList) 
solnList(numsolves) = struct('gp',[],'snr',[],'numberuncertain',[],'FaList',[],'solver',[],'objective',[],'plotlabel',[], 'params', [], 'Mxy', [], 'Mz', [],'signuImage',[],'signu',[],'MIval',[]);
icount  = 0;
%% for iobj = 1:numel(ObjectiveType)
%%  for isolver = 1:numel(solverType)
%%   for igp = 1:length(gpList)
%%    for inu = 1:length(uncertainList)
%%     for isnr = 1:length(snrList)
%%        worktmp = load(sprintf('poptNG%dNu%d%s%sSNR%02dHermite.mat',gpList(igp),uncertainList(inu),solverType{isolver},ObjectiveType{iobj},snrList(isnr)));
%%        icount= icount+1;
%%        solnList (icount) = struct('gp',gpList(igp),'snr',snrList(isnr),'numberuncertain',uncertainList(inu),'FaList',worktmp.popt.FaList,'solver',solverType{isolver},'objective',ObjectiveType{iobj},'plotlabel',sprintf('%s%s',solverType{isolver},ObjectiveType{iobj}),'params',worktmp.params, 'Mxy',worktmp.Mxyref, 'Mz',worktmp.Mzref,'signuImage',worktmp.signuImage,'signu',worktmp.signu,'MIval',worktmp.fval);
%%     end
%%    end
%%   end
%%  end
%% end

% compute variance for each SNR for 20/30 pyr/lac FA
hackuncertainList=3
hackgpList=5
for isnr = 1:length(snrList)
   icount= icount+1;
   worktmp = load(sprintf('poptNG%dNu%d%s%sSNR%02dHermite.mat',hackgpList,hackuncertainList,'constDirect','TotalSignal',snrList(isnr)));
   solnList (icount) = struct('gp',-1,'snr',snrList(isnr),'numberuncertain',-1,'FaList',worktmp.params.FaList,'solver','const','objective','Max','plotlabel','constMax','params',worktmp.params, 'Mxy',worktmp.Mxy, 'Mz',worktmp.Mz,'signuImage',worktmp.signuImage,'signu',worktmp.signu,'MIval',NaN);
end

% array info
Ntime    = size(solnList(1).Mz,2);
Nspecies = size(solnList(1).Mz,1);

% establish control
cntrlparams = worktmp.params;
cntrlparams.FaList = ones(size(cntrlparams.FaList))*pi/180;
M0 = [0,0];
model = HPKinetics.NewMultiPoolTofftsGammaVIF();
[t_axiscntrl,Mxycntrl,Mzcntrl] = model.compile(M0.',cntrlparams);

for isnr = 1:length(snrList)
   icount= icount+1;
   solnList (icount) = struct('gp',-1,'snr',snrList(isnr),'numberuncertain',-1,'FaList',cntrlparams.FaList,'solver','control','objective','MI','plotlabel','control','params',cntrlparams, 'Mxy',Mxycntrl, 'Mz',Mzcntrl,'signuImage',solnList(isnr).signuImage,'signu',solnList(isnr).signu,'MIval',NaN);
end


% pareto total signal vs MI
for jjj = 1:size(myFAList,2)
  paretoparams = worktmp.params;
  paretoparams.FaList = repmat(myFAList(:,jjj),1,Ntime )*pi/180;
  M0 = [0,0];
  model = HPKinetics.NewMultiPoolTofftsGammaVIF();
  [t_axispareto,Mxypareto,Mzpareto] = model.compile(M0.',paretoparams);
  
  
  for isnr = 1:length(snrList)
     icount= icount+1;
     solnList (icount) = struct('gp',-1,'snr',snrList(isnr),'numberuncertain',-1,'FaList',paretoparams.FaList,'solver',sprintf('paretoP%dL%d',myFAList(:,jjj)),'objective','Max','plotlabel',sprintf('paretoP%dL%dMax',myFAList(:,jjj)),'params',paretoparams, 'Mxy',Mxypareto, 'Mz',Mzpareto,'signuImage',solnList(isnr).signuImage,'signu',solnList(isnr).signu,'MIval',NaN);
  end
end

%UB/LB for MI solution
hacksolvertype='interior-point'
for isnr = 1:length(snrList)
   worktmpLB = load(sprintf('poptNG%dNu%d%s%sSNR%02dHermite.mat',hackgpList,hackuncertainList,hacksolvertype,ObjectiveType{1},snrList(1)));
   icount= icount+1;
   solnList (icount) = struct('gp',hackgpList,'snr',snrList(isnr),'numberuncertain',hackuncertainList,'FaList',worktmp.popt.FaList,'solver',hacksolvertype,'objective',ObjectiveType{1},'plotlabel',sprintf('%s%sLB',hacksolvertype,ObjectiveType{1} ),'params',worktmp.params, 'Mxy',worktmp.Mxyref, 'Mz',worktmp.Mzref,'signuImage',solnList(isnr).signuImage,'signu',solnList(isnr).signu,'MIval',worktmp.fval);
end

for isnr = 1:length(snrList)
   worktmpUB = load(sprintf('poptNG%dNu%d%s%sSNR%02dHermite.mat',hackgpList,hackuncertainList,hacksolvertype,ObjectiveType{1},snrList(end)));
   icount= icount+1;
   solnList (icount) = struct('gp',hackgpList,'snr',snrList(isnr),'numberuncertain',hackuncertainList,'FaList',worktmp.popt.FaList,'solver',hacksolvertype,'objective',ObjectiveType{1},'plotlabel',sprintf('%s%sUB',hacksolvertype,ObjectiveType{1} ),'params',worktmp.params, 'Mxy',worktmp.Mxyref, 'Mz',worktmp.Mzref,'signuImage',solnList(isnr).signuImage,'signu',solnList(isnr).signu,'MIval',worktmp.fval);
end

% extract timehistory info
num_trials = 25;
timehistory  = zeros(Ntime,Nspecies,num_trials+1,length(solnList) );

for jjj =1:length(solnList)
  % NOTE - image noise is at the single image for single species - signu is for the sum over time for both species ==> divide by Ntime and Nspecies
  imagenoise = solnList(jjj).signuImage;
  %disp([xroi(jjj),yroi(jjj),zroi(jjj)]);
  timehistory(:,1,num_trials+1,jjj) = solnList(jjj).Mxy(1,:);
  timehistory(:,2,num_trials+1,jjj) = solnList(jjj).Mxy(2,:);
  % add noise for num_trials
  for kkk = 1:num_trials
      pyrnoise = imagenoise *randn(Ntime,1);
      timehistory(:,1,kkk,jjj)= timehistory(:,1,num_trials+1,jjj) + pyrnoise;
      lacnoise = imagenoise *randn(Ntime,1);
      timehistory(:,2,kkk,jjj)= timehistory(:,2,num_trials+1,jjj) + lacnoise;
  end
end

%% TR_list = solnList(1).params.TRList 
%% figure(1)
%% plot( TR_list, timehistory(:,1,1,1), 'b', ...
%%       TR_list, timehistory(:,2,1,1), 'b-.', ...
%%       TR_list, timehistory(:,1,1,2), 'g', ...
%%       TR_list, timehistory(:,2,1,2), 'g-.', ...
%%       TR_list, timehistory(:,1,1,3), 'r', ...
%%       TR_list, timehistory(:,2,1,3), 'r-.', ...
%%       TR_list, timehistory(:,1,1,4), 'c', ...
%%       TR_list, timehistory(:,2,1,4), 'c-.')
%% title('noise added')
%% figure(2)
%% plot( TR_list, timehistory(:,1,num_trials+1,1), 'b', ...
%%       TR_list, timehistory(:,2,num_trials+1,1), 'b-.', ...
%%       TR_list, timehistory(:,1,num_trials+1,2), 'g', ...
%%       TR_list, timehistory(:,2,num_trials+1,2), 'g-.', ...
%%       TR_list, timehistory(:,1,num_trials+1,3), 'r', ...
%%       TR_list, timehistory(:,2,num_trials+1,3), 'r-.', ...
%%       TR_list, timehistory(:,1,num_trials+1,4), 'c', ...
%%       TR_list, timehistory(:,2,num_trials+1,4), 'c-.')
%% title('original data')

%% store soln
storekplopt   = zeros(num_trials+1,length(solnList));
storekveqpopt = zeros(num_trials+1,length(solnList));
storeT1Popt   = zeros(num_trials+1,length(solnList));
storeT1Lopt   = zeros(num_trials+1,length(solnList));
storet0opt    = zeros(num_trials+1,length(solnList));

% num optimization vars
numberParameters = 3

loadDataHistory = true
if loadDataHistory 
  load('recovervalidate.mat') 
  idplot = (num_trials+1)*length(solnList);
else
  for idesign = 1:length(solnList)
  %for idesign =13:16
     % setup optimization variables
     switch (numberParameters)
         case(1) 
           kpl   = optimvar('kpl','LowerBound',0);
           kveqp = solnList(idesign ).params.PerfusionTerms(1) ;
           T1P   = solnList(idesign ).params.T1s(1);
           T1L   = solnList(idesign ).params.T1s(2);
           t0    = solnList(idesign ).params.t0(1);
           % ground truth 
           xstar.kpl        = solnList(idesign).params.ExchangeTerms(1,2);
         case(2) 
           kpl   = optimvar('kpl','LowerBound',0);
           kveqp = optimvar('kveqp','LowerBound',0);
           T1P   = solnList(idesign ).params.T1s(1);        
           T1L   = solnList(idesign ).params.T1s(2);        
           t0    = solnList(idesign ).params.t0(1);
           % ground truth 
           xstar.kpl        = solnList(idesign).params.ExchangeTerms(1,2); 
           xstar.kveqp      = solnList(idesign ).params.PerfusionTerms(1) ;
         case(3) 
           kpl   = optimvar('kpl','LowerBound',0);
           kveqp = optimvar('kveqp','LowerBound',0);
           t0  = optimvar('t0' ,'LowerBound',0);
           T1P   = solnList(idesign ).params.T1s(1);        
           T1L   = solnList(idesign ).params.T1s(2);        
           % ground truth 
           xstar.kpl        = solnList(idesign).params.ExchangeTerms(1,2); 
           xstar.kveqp      = solnList(idesign ).params.PerfusionTerms(1) ; 
           xstar.t0         = solnList(idesign ).params.t0(1);
         case(4) 
           kpl   = optimvar('kpl','LowerBound',0);
           kveqp = optimvar('kveqp','LowerBound',0);
           T1P   = optimvar('T1P','LowerBound',0);
           T1L   = optimvar('T1L','LowerBound',0);
           t0    = solnList(idesign ).params.t0(1);
           % ground truth 
           xstar.kpl        = solnList(idesign).params.ExchangeTerms(1,2);
           xstar.kveqp      = solnList(idesign ).params.PerfusionTerms(1);
           xstar.T1P        = solnList(idesign ).params.T1s(1);
           xstar.T1L        = solnList(idesign ).params.T1s(2);
         case(5) 
           kpl = optimvar('kpl','LowerBound',0);
           kveqp = optimvar('kveqp','LowerBound',0);
           T1P = optimvar('T1P','LowerBound',0);
           T1L = optimvar('T1L','LowerBound',0);
           t0  = optimvar('t0' ,'LowerBound',0);
           % ground truth 
           xstar.kpl        = solnList(idesign).params.ExchangeTerms(1,2); 
           xstar.kveqp      = solnList(idesign ).params.PerfusionTerms(1);
           xstar.T1P        = solnList(idesign ).params.T1s(1);        
           xstar.T1L        = solnList(idesign ).params.T1s(2);        
           xstar.t0         = solnList(idesign ).params.t0(1);
     end
     statevariable  = optimexpr([Nspecies,Ntime]);
     statesignal    = optimexpr([Nspecies,Ntime]);
  
  
     disp('build constraints')
     A0    = solnList(idesign ).params.scaleFactor(1);
     statevariable(:,1 )=0;
     statesignal(:,1 )=0;
     for iii = 1:Ntime-1
       disp([Nspecies*(iii-1)+1 , Nspecies*(iii-1)+2 ,Ntime* Nspecies])
       klpqp =    0 ;     % @cmwalker where do I get this from ? 
       % 
       currentTR = solnList(idesign ).params.TRList(iii+1) - solnList(idesign ).params.TRList(iii);
       nsubstep = 5;
       deltat = currentTR /nsubstep ;
       integratedt = [solnList(idesign ).params.TRList(iii):deltat:solnList(idesign ).params.TRList(iii+1)] +deltat/2  ;
       integrand = fcn2optimexpr(@(bolusstart)A0*gampdf(integratedt(1:nsubstep )'- bolusstart ,solnList(idesign ).params.gammaPdfA(1),solnList(idesign ).params.gammaPdfB(1)) ,t0);
       % >> syms a  kpl d currentTR    T1P kveqp T1L 
       % >> expATR = expm([a,  0; kpl, d ] * currentTR )
       % 
       % expATR =
       % [                                     exp(a*currentTR),                0]
       % [(kpl*exp(a*currentTR) - kpl*exp(currentTR*d))/(a - d), exp(currentTR*d)]
       % 
       % >> a = -1/T1P - kpl - kveqp
       % >> d = -1/T1L
       % >> eval(expATR)
       % ans =
       % [                                                              exp(-currentTR*(kpl + kveqp + 1/T1P)),                   0]
       % [(kpl*exp(-currentTR/T1L) - kpl*exp(-currentTR*(kpl + kveqp + 1/T1P)))/(kpl + kveqp - 1/T1L + 1/T1P), exp(-currentTR/T1L)]
       %    
       %expATR = fcn2optimexpr(@expm,A*currentTR );
       % A = [-1/T1P - kpl - kveqp,  0; kpl, -1/T1L ];
       expATR = [ exp(-currentTR*(kpl + kveqp/solnList(idesign ).params.volumeFractions + 1/T1P)),                   0; (kpl*exp(-currentTR/T1L) - kpl*exp(-currentTR*(kpl + kveqp/solnList(idesign ).params.volumeFractions + 1/T1P)))/(kpl + kveqp/solnList(idesign ).params.volumeFractions - 1/T1L + 1/T1P), exp(-currentTR/T1L)];
       % mid-point rule integration
       aifterm = kveqp/solnList(idesign ).params.volumeFractions * deltat * [ exp((-1/T1P - kpl - kveqp/solnList(idesign ).params.volumeFractions)*(solnList(idesign ).params.TRList(iii+1)-deltat*[.5:1:nsubstep]-solnList(idesign ).params.TRList(iii)) );
     (kpl*exp((-1/T1P - kpl - kveqp/solnList(idesign ).params.volumeFractions)*(solnList(idesign ).params.TRList(iii+1)-deltat*[.5:1:nsubstep]-solnList(idesign ).params.TRList(iii)) ) - kpl*exp(-1/T1L *(solnList(idesign ).params.TRList(iii+1)-deltat*[.5:1:nsubstep]-solnList(idesign ).params.TRList(iii))))/((-1/T1P - kpl - kveqp/solnList(idesign ).params.volumeFractions) + 1/T1L )] * integrand ;
       statevariable(:,iii+1) =  expATR *( statevariable(:,iii ))   + aifterm     ;
       % 
       % TODO TODO - NOTE  that a simple change to minimize the expression within fcn2optimexpr is the difference between the code running or NOT
       % DOES NOT RUN with this EXPRESSION: statesignal(:,iii+1)   =  fcn2optimexpr(@(newbolusstart)sin(solnList(idesign ).FaList(:,iii+1)).* ( solnList(idesign ).params.volumeFractions * statevariable(:,iii+1) +  (1-solnList(idesign ).params.volumeFractions) * [A0;0] * gampdf( solnList(idesign ).params.TRList(iii+1) - newbolusstart  , solnList(idesign ).params.gammaPdfA(1) , solnList(idesign ).params.gammaPdfB(1) ) ),t0 ) ;
       statesignal(:,iii+1)   =  sin(solnList(idesign ).FaList(:,iii+1)).* ( solnList(idesign ).params.volumeFractions * statevariable(:,iii+1) +  (1-solnList(idesign ).params.volumeFractions) * [A0;0] * fcn2optimexpr(@(newbolusstart)gampdf( solnList(idesign ).params.TRList(iii+1) - newbolusstart  , solnList(idesign ).params.gammaPdfA(1) , solnList(idesign ).params.gammaPdfB(1) ),t0 ) ) ;
       statevariable(:,iii+1) =  cos(solnList(idesign ).FaList(:,iii+1)).* statevariable(:,iii+1);
     end
     truthstate  = evaluate(statevariable ,xstar);
     truthsignal = evaluate(statesignal ,xstar);
  
     %% 
     % Create an optimization problem using these converted optimization expressions.
     disp('build objective function')
     %  for idtrial = num_trials+1:num_trials+1
     for idtrial = 1:num_trials+1
         tic;
         disp(sprintf('design = %d, trial = %d',idesign, idtrial )); solnList(idesign)
         mycostfcn = sum( (statesignal(1,:)'- timehistory(:,1,idtrial,idesign ) ).^2) + sum( (statesignal(2,:)'- timehistory(:,2,idtrial,idesign ) ).^2);
  
         disp('create optim prob')
         convprob = optimproblem('Objective',mycostfcn );
         %% 
         % View the new problem.
         
         %show(convprob)
         % problem = prob2struct(convprob,'ObjectiveFunctionName','generatedObjectiveRecover');
         %% 
         % Solve the new problem. The solution is essentially the same as before.
  
         %myoptions = optimoptions(@fmincon,'Display','iter-detailed','SpecifyObjectiveGradient',true, 'SpecifyConstraintGradient', true)
         %myoptions = optimoptions(@fmincon,'Display','iter-detailed','MaxFunctionEvaluations' , 3.000000e+04)
         myoptions = optimoptions(@lsqnonlin,'Display','iter-detailed');
             
         %[popt,fval,exitflag,output] = solve(convprob,x0,'Options',myoptions, 'ObjectiveDerivative', 'auto-reverse' , 'ConstraintDerivative', 'auto-reverse', 'solver', 'fmincon' )
         %[popt,fval,exitflag,output] = solve(convprob,x0,'Options',myoptions, 'ObjectiveDerivative', 'finite-differences','solver', 'fmincon' )
         % random initial guess
         switch (numberParameters)
            case(1) 
              x0.kpl        = unifrnd(.06,.24);
            case(2) 
              x0.kpl        = unifrnd(.01,10);
              x0.kveqp      = unifrnd(.02,10);
            case(3) 
              x0.kpl        = unifrnd(.06,.24);
              x0.kveqp      = unifrnd(.02,.08);
              x0.t0         = unifrnd(0  ,8);
            case(4) 
              x0.kpl        = unifrnd(.01,10);
              x0.kveqp      = unifrnd(.02,10);
              x0.T1P        = unifrnd(20 ,40);
              x0.T1L        = unifrnd(15 ,35);
            case(5) 
              x0.kpl        = unifrnd(.01,10);
              x0.kveqp      = unifrnd(.02,10);
              x0.T1P        = unifrnd(20 ,40);
              x0.T1L        = unifrnd(15 ,35);
              x0.t0         = unifrnd(0  ,8);
         end
         initialstate  = evaluate(statevariable ,x0);
         initialsignal = evaluate(statesignal ,x0);
         [popt,fval,exitflag,output] = solve(convprob,x0,'Options',myoptions, 'solver','lsqnonlin' , 'ObjectiveDerivative', 'finite-differences')
         switch (numberParameters)
            case(1) 
              storekplopt(  idtrial,idesign)  = popt.kpl;
              storekveqpopt(idtrial,idesign)  = solnList(idesign ).params.PerfusionTerms(1);
              storeT1Popt(  idtrial,idesign)  = solnList(idesign ).params.T1s(1);
              storeT1Lopt(  idtrial,idesign)  = solnList(idesign ).params.T1s(2);
              storet0opt (  idtrial,idesign)  = solnList(idesign ).params.t0(1);
            case(2) 
              storekplopt(  idtrial,idesign)  = popt.kpl;
              storekveqpopt(idtrial,idesign)  = popt.kveqp;
              storeT1Popt(  idtrial,idesign)  = solnList(idesign ).params.T1s(1); 
              storeT1Lopt(  idtrial,idesign)  = solnList(idesign ).params.T1s(2); 
              storet0opt (  idtrial,idesign)  = solnList(idesign ).params.t0(1);
            case(3) 
              storekplopt(  idtrial,idesign)  = popt.kpl;
              storekveqpopt(idtrial,idesign)  = popt.kveqp;
              storeT1Popt(  idtrial,idesign)  = solnList(idesign ).params.T1s(1); 
              storeT1Lopt(  idtrial,idesign)  = solnList(idesign ).params.T1s(2); 
              storet0opt (  idtrial,idesign)  = popt.t0 ; 
            case(4) 
              storekplopt(  idtrial,idesign)  = popt.kpl;
              storekveqpopt(idtrial,idesign)  = popt.kveqp;
              storeT1Popt(  idtrial,idesign)  = popt.T1P;
              storeT1Lopt(  idtrial,idesign)  = popt.T1L;
              storet0opt (  idtrial,idesign)  = solnList(idesign ).params.t0(1);
            case(5) 
              storekplopt(  idtrial,idesign)  = popt.kpl;
              storekveqpopt(idtrial,idesign)  = popt.kveqp;
              storeT1Popt(  idtrial,idesign)  = popt.T1P;
              storeT1Lopt(  idtrial,idesign)  = popt.T1L;
              storet0opt (  idtrial,idesign)  = popt.t0 ; 
         end
         slnstate  = evaluate(statevariable ,popt);
         slnsignal = evaluate(statesignal ,popt);
  
         if idtrial >= num_trials 
         % plot
         idplot = (num_trials+1)*(idesign-1) + idtrial ;
         handle = figure(idplot )
         plot(solnList(idesign ).params.TRList , timehistory(:,1,idtrial,idesign), 'b', ...
              solnList(idesign ).params.TRList , timehistory(:,2,idtrial,idesign), 'b-.', ...
              solnList(idesign ).params.TRList , timehistory(:,1,num_trials+1,idesign), 'r', ...
              solnList(idesign ).params.TRList , timehistory(:,2,num_trials+1,idesign), 'r-.', ...
              solnList(idesign ).params.TRList , initialsignal(1,:), 'g', ...
              solnList(idesign ).params.TRList , initialsignal(2,:), 'g-.', ...
              solnList(idesign ).params.TRList , truthsignal(1,:), 'm', ...
              solnList(idesign ).params.TRList , truthsignal(2,:), 'm-.', ...
              solnList(idesign ).params.TRList , slnsignal(1,:), 'k', ...
              solnList(idesign ).params.TRList , slnsignal(2,:), 'k-.')
         ylabel('Mxy')
         xlabel('sec')
         title(sprintf('curvefit %s %d %d',solnList(idesign ).solver,solnList(idesign ).snr,idtrial))
         legend('walker+rice','','walker','','ic','','truth','','df','')
         set(gca,'FontSize',16)
         saveas(handle,sprintf('recoverparametersFig%04d',idplot),'png')
         pause(.1)
         end 
  
         toc;
     end 
  
  end 
  
  %save workspace
  save('recovervalidate.mat' ,'storekplopt', 'storekveqpopt', 'storeT1Popt', 'storeT1Lopt', 'storet0opt')
end
           
%   Various line types, plot symbols and colors may be obtained with
%   PLOT(X,Y,S) where S is a character string made from one element
%   from any or all the following 3 columns:
%          b     blue          .     point              -     solid
%          g     green         o     circle             :     dotted
%          r     red           x     x-mark             -.    dashdot
%          c     cyan          +     plus               --    dashed
%          m     magenta       *     star             (none)  no line
%          y     yellow        s     square
%          k     black         d     diamond
%          w     white         v     triangle (down)
%                              ^     triangle (up)
%                              <     triangle (left)
%                              >     triangle (right)
%                              p     pentagram
%                              h     hexagram

numplots = length(solnList) / length(snrList)
textscale = .01 * ones(numplots ,length(snrList));
textscale(1,1) = .015;
textscale(4,1) = .02;
myupperb = .24 * ones(numplots ,length(snrList));
myupperb(5,2) = .20;
myupperb(6,2) = .20;
myplottitle = { solnList(0*length(snrList)+1).plotlabel solnList(1*length(snrList)+1).plotlabel solnList(2*length(snrList)+1).plotlabel solnList(3*length(snrList)+1).plotlabel solnList(4*length(snrList)+1).plotlabel solnList(5*length(snrList)+1).plotlabel } 
myplotlabel = { "" "" "" "" "" "" } 
myplottitle(1) = cellstr("clinic")
myplottitle(3) = cellstr("constant")
myplottitle(4) = cellstr("constant")
myplottitle(5) = cellstr("varying")
myplottitle(6) = cellstr("varying")
myplotlabel(3) = cellstr("K_{OED_{20}}")
myplotlabel(4) = cellstr("K_{OED_{2}}")
myplotlabel(5) = cellstr("K_{OED_{2}}")
myplotlabel(6) = cellstr("K_{OED_{20}}")
for isolver = 1:numplots 
  idplot = idplot+1
  handle = figure(idplot )
  inversestd  = std(storekplopt(1:num_trials,:),0,1)
  inversemean = mean(storekplopt(1:num_trials,:),1)
  plot( snrList , ones(1,length(snrList))*solnList(end).params.ExchangeTerms(1,2) , 'g','linewidth',2) 
  hold
  errorbar(snrList,inversemean((isolver-1)*length(snrList)+1:isolver*length(snrList)),2*inversestd((isolver-1)*length(snrList)+1:isolver*length(snrList)),'s','LineStyle','none', 'Color', 'b','linewidth', 2);
  ylabel('fit kpl (sec^{-1})')
  xlabel('SNR_{data}')
  xlim([0 30])
  ylim([0.05 .25])
  %textscale = max(inversestd((isolver-1)*length(snrList)+1:isolver*length(snrList)))*.1;
  text(snrList,min([inversemean((isolver-1)*length(snrList)+1:isolver*length(snrList));myupperb(isolver,:)])+textscale(isolver,:), sprintfc('\\mu=%6.4f',inversemean((isolver-1)*length(snrList)+1:isolver*length(snrList))) )
  text(snrList,min([inversemean((isolver-1)*length(snrList)+1:isolver*length(snrList));myupperb(isolver,:)])-textscale(isolver,:), sprintfc('\\sigma=%6.4f',inversestd( (isolver-1)*length(snrList)+1:isolver*length(snrList))) )
  text(25,.1, string(myplotlabel(isolver))  )
  title(myplottitle(isolver) )
  legend('truth','fit')
  set(gca,'FontSize',16)
  saveas(handle,sprintf('solversummaryNP%d%s',numberParameters,solnList((isolver-1)*length(snrList)+1).plotlabel),'png')
end

idplot = idplot+1
handle = figure(idplot )
boxplot(  storekplopt(1:num_trials,:), {solnList(:).plotlabel} )
ylim([0 .4])
saveas(handle,'globalboxplot','png')
 
    
% create legend for index
%solnTable = struct2table(solnList)
%solnTable.myindex = [1:60]'

% analysis of variance
constDirectSumQuadSNR25 = 16
constMaxSNR25 = 28
paretoP20L28MaxSNR25 = 52
constDirectSumQuadSNR20 = 15
constMaxSNR20 = 27
paretoP20L28MaxSNR20 = 51

%% % test tha variance of constDirectSumQuadSNR25  is less than constMaxSNR25 
%% [myh,myp,myci,mystats] = vartest2(storekplopt(1:num_trials,constDirectSumQuadSNR25 ),storekplopt(1:num_trials,constMaxSNR25 ),'Tail','left')
%% [myh2,myp2,myci2,mystats2] = vartest2(storekplopt(1:num_trials,paretoP20L28MaxSNR25 ),storekplopt(1:num_trials,constMaxSNR25 ),'Tail','both')
%% [myh3,myp3,myci3,mystats3] = vartest2(storekplopt(1:num_trials,paretoP20L28MaxSNR20 ),storekplopt(1:num_trials,constMaxSNR20 ),'Tail','both')
%% std(storekplopt(1:num_trials,paretoP20L28MaxSNR25 ))
%% std(storekplopt(1:num_trials,constMaxSNR25 ))
%% std(storekplopt(1:num_trials,constDirectSumQuadSNR25 ))
%% 
% A paired samples t-test was performed to compare miles per gallon between fuel treatment and no fuel treatment.  There was a significant difference in miles per gallon between fuel treatment (M = 22.75, SD = 3.25) and no fuel treatment (M = 21, SD = 2.73); t(11) = -2.244, p = .046.

% One sample t-test was conducted to determine whether there is a difference between the results on the Math test and the true population mean (M=86.00). The results indicate a significant difference between the true mean (M=86.00) and the mean Math test score (M=73.08; SD=16.89), [t(36) = -4.651, p = .000]. We, therefore, reject the null hypothesis that there is not a difference between the true mean and the comparison value and conclude that our mean Math test score is significantly different from the true population mean.
