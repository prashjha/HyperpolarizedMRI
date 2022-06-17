% ShowMxyPub Shows the flip angle schedules and the resulting magnetization
% evolutions for the paper
% Author: Chris Walker
% Date: 8/6/2018

clear all
clc

mynewoptions.Algorithm = 'constDirect'
driverHPMIopt(5,3, 2,mynewoptions,'SumQuad',false)
driverHPMIopt(5,3, 5,mynewoptions,'SumQuad',false)
driverHPMIopt(5,3,10,mynewoptions,'SumQuad',false)
driverHPMIopt(5,3,15,mynewoptions,'SumQuad',false)
driverHPMIopt(5,3,20,mynewoptions,'SumQuad',false)
driverHPMIopt(5,3, 2,mynewoptions,'TotalSignal',false)
driverHPMIopt(5,3, 5,mynewoptions,'TotalSignal',false)
driverHPMIopt(5,3,10,mynewoptions,'TotalSignal',false)
driverHPMIopt(5,3,15,mynewoptions,'TotalSignal',false)
driverHPMIopt(5,3,20,mynewoptions,'TotalSignal',false)
%
myoptions = optimoptions(@fmincon,'Display','iter-detailed','SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,'MaxFunctionEvaluations',1e7,'ConstraintTolerance',2.e-9, 'OptimalityTolerance',2.5e-4,'Algorithm','interior-point','StepTolerance',1.000000e-12,'MaxIterations',1000,'PlotFcn',{'optimplotfvalconstr', 'optimplotconstrviolation', 'optimplotfirstorderopt' },'HonorBounds',true, 'Diagnostic','on','FunValCheck','on' )
%%driverHPMIopt(3,3,20,myoptions,'TotalSignal',false)
driverHPMIopt(5,3, 2,myoptions,'TotalSignal',false)
driverHPMIopt(5,3, 5,myoptions,'TotalSignal',false)
driverHPMIopt(5,3,10,myoptions,'TotalSignal',false)
driverHPMIopt(5,3,15,myoptions,'TotalSignal',false)
driverHPMIopt(5,3,20,myoptions,'TotalSignal',false)
driverHPMIopt(5,3, 2,myoptions,'SumQuad',false)
driverHPMIopt(5,3, 5,myoptions,'SumQuad',false)
driverHPMIopt(5,3,10,myoptions,'SumQuad',false)
driverHPMIopt(5,3,15,myoptions,'SumQuad',false)
driverHPMIopt(5,3,20,myoptions,'SumQuad',false)

%% % monitor memory: while [ -e /proc/3291925 ] ; do  top -b -n 1 -p 3291925 >>process.txt ;sleep 60; done  

function driverHPMIopt(NGauss,NumberUncertain,modelSNR,myoptions,ObjectiveType,GaussLegendre )
  %NGauss = 3,NumberUncertain=3,modelSNR=10, ObjectiveType = 'SumQuad'

  NGauss,NumberUncertain,modelSNR,myoptions.Algorithm,ObjectiveType,GaussLegendre
  close all
  %% Tissue Parameters
  T1pmean = [ 30 ]; % s
  T1pstdd = [ 10 ]; % s
  T1plb   = [ 5  ]; % s
  T1pub   = [ 45 ]; % s
  T1lmean = [ 25 ]; % s
  T1lstdd = [ 10 ]; % s
  T1llb   = [  5 ]; % s
  T1lub   = [ 45 ]; % s
  kplmean = [ .15 ];       % s
  kplstdd = [ .03 ];       % s
  kpllb   = [ .01 ];       % s
  kplub   = [ .35 ];       % s
  kvemean = [ 0.05 ];      % s
  kvestdd = [ .01  ];      % s
  kvelb   = [ 0.01 ];      % s
  kveub   = [ 0.20 ];      % s
  t0mean  = [ 4    ];      % s
  t0sttd  = [ 1.3  ];      % s
  t0lb    = [ 0    ];       % s
  t0ub    = [ 7    ];       % s
  alphamean  =  [2.5];
  alphasttd  =  [.3];
  betamean  =  [4.5];
  betasttd  =  [.3];
  tisinput=[T1pmean; T1pstdd; T1lmean; T1lstdd; kplmean; kplstdd; kvemean; kvestdd;t0mean;t0sttd;alphamean; alphasttd; betamean ; betasttd ];
  tisinputlbub=[T1plb; T1pub; T1llb;T1lub; kpllb; kplub; kvelb; kveub;t0lb;t0ub;alphamean-2*alphasttd;alphamean+2*alphasttd; betamean-2*betasttd ;betamean+2*betasttd ];
  
  %% Variable Setup
  Ntime = 30;
  currentTR = 3;
  TR_list = (0:(Ntime-1))*currentTR ;
  M0 = [0,0];
  ve = 0.95;
  %ve = 1.;
  VIF_scale_fact = [100;0];
  bb_flip_angle = 20;
  opts = optimset('lsqcurvefit');
  opts.TolFun = 1e-09;
  opts.TolX = 1e-09;
  opts.Display = 'off';
  params = struct('t0',[t0mean(1);0],'gammaPdfA',[alphamean(1)  ;1],'gammaPdfB',[betamean(1);1],...
      'scaleFactor',VIF_scale_fact,'T1s',[T1pmean(1),T1lmean(1)],'ExchangeTerms',[0,kplmean(1) ;0,0],...
      'TRList',TR_list,'PerfusionTerms',[kvemean(1),0],'volumeFractions',ve,...
      'fitOptions', opts)
  model = HPKinetics.NewMultiPoolTofftsGammaVIF();
  
  
  %% Get true Mz
  %% Choose Excitation Angle
  FAType = {'Const'};
  %% HACK- @cmwalker code for initial conditions - https://github.com/fuentesdt/TumorHPMRI/blob/master/models/gPC/walker/ShowMxyPub.m
  for i = 1:numel(FAType)
      switch (FAType{i})
          case('Const') % Nagashima for lactate const 10 pyruvate
              tic
              E1(1) = exp(-currentTR *(1/T1pmean+kplmean));
              E1(2) = exp(-currentTR /T1lmean);
              for n = 1:Ntime
                  % 20deg for pyruvate 30deg for lactate - currently used in brain
                  flips(2,n) = 30*pi/180;
                  flips(1,n) = 20*pi/180;
              end
              params.FaList = flips ;
      end
  
      
      
  
      tic
      %% Fitting
      [t_axis,Mxy,Mz] = model.compile(M0.',params);
      toc
  end
  
  %% Plot initial guess
  plotinit = true;
  if plotinit
      % plot initial guess
      figure(1)
      plot(TR_list,Mxy(1,:),'b',TR_list,Mxy(2,:),'k')
      ylabel('Const Mxy')
      xlabel('sec')
      figure(2)
      plot(TR_list,params.FaList(1,:)*180/pi,'b',TR_list,params.FaList(2,:)*180/pi,'k')
      ylabel('Const FA (deg) ')
      xlabel('sec')
      figure(3)
      plot(TR_list,Mz(1,:),'b--',TR_list,Mz(2,:),'k--')
      hold
      plot(TR_list,Mz(1,:)./cos(params.FaList(1,:)),'b',TR_list,Mz(2,:)./cos(params.FaList(2,:)),'k')
      ylabel('Const Mz')
      xlabel('sec')
  
      % plot gamma
      jmA0    = VIF_scale_fact(1);
      jmalpha = alphamean(1);
      jmbeta  = betamean(1);
      jmt0    = t0mean(1);
      jmaif   = jmA0  * gampdf(TR_list - jmt0  , jmalpha , jmbeta);
      figure(4)
      plot(TR_list,jmaif ,'b')
      ylabel('aif')
      xlabel('sec')
  end
  
  
  %% optimize MI for TR and FA
  optf = true;
  if optf
  
      tic;
      % setup optimization variables
      Nspecies = 2
      FaList = optimvar('FaList',Nspecies,Ntime,'LowerBound',0, 'UpperBound',35*pi/180);
      TRList = TR_list;
      diffTR = diff(TRList);
  
      % noise calc for signal sum - assume same noise in pyruvate and lactate image
      %    modelSNR = (maxsignallac + maxsignalpyr)/(stdsignalpyr +stdsignallac  )
      %             = (maxsignallac + maxsignalpyr)/(2 * signuImage )
      %
      % image level ==> signuImage = (maxsignallac + maxsignalpyr)/2/modelSNR 
      % signal sum  ==> signu = Ntime * signuImage 
      % noise calc max signal assuming total signal is sum of gaussian RV
      signuImage = (max(Mxy(1,:))+max(Mxy(2,:)))/2/modelSNR;
      % walker paper is peak pyruvate only
      signuImage = max(Mxy(1,:))/modelSNR;
      % variance for Gauss RV is sum. sqrt for std
      signu = sqrt(2* Ntime) * signuImage;
      [x2,xn2,xm2,w2,wn2]=GaussHermiteNDGauss(NGauss,0,signu);
      lqp2=length(xn2{1}(:));
  
      % switch between uniform and Gaussian RV for prior
      if(GaussLegendre)
        QuadratureRule = 'Legendre';
      else
        QuadratureRule = 'Hermite';
      end
      switch (NumberUncertain)
         case(3)
           if(GaussLegendre)
             [x,xn,xm,w,wn]=GaussLegendreNDGauss(NGauss,[tisinputlbub(5:2:9)],[tisinputlbub(6:2:10)]);
           else
             [x,xn,xm,w,wn]=GaussHermiteNDGauss(NGauss,[tisinput(5:2:9)],[tisinput(6:2:10)]);
           end
           T1Pqp   = T1pmean;
           T1Lqp   = T1lmean;
           kplqp   = xn{1}(:);
           klpqp   =    0 ;     % @cmwalker where do I get this from ? 
           kveqp   = xn{2}(:);
           t0qp    = xn{3}(:); 
         case(4)
           if(GaussLegendre)
             [x,xn,xm,w,wn]=GaussLegendreNDGauss(NGauss,[tisinputlbub(1:2:7)],[tisinputlbub(2:2:8)]);
           else
             [x,xn,xm,w,wn]=GaussHermiteNDGauss(NGauss,[tisinput(1:2:7)],[tisinput(2:2:8)]);
           end
           T1Pqp   = xn{1}(:);
           T1Lqp   = xn{2}(:);
           kplqp   = xn{3}(:);
           klpqp   =    0 ;     % @cmwalker where do I get this from ? 
           kveqp   = xn{4}(:);
           t0qp    = t0mean(1); 
         case(5)
           if(GaussLegendre)
             [x,xn,xm,w,wn]=GaussLegendreNDGauss(NGauss,[tisinputlbub(1:2:9)],[tisinputlbub(2:2:10)]);
           else
             [x,xn,xm,w,wn]=GaussHermiteNDGauss(NGauss,[tisinput(1:2:9)],[tisinput(2:2:10)]);
           end
           T1Pqp   = xn{1}(:);
           T1Lqp   = xn{2}(:);
           kplqp   = xn{3}(:);
           klpqp   =    0 ;     % @cmwalker where do I get this from ? 
           kveqp   = xn{4}(:);
           t0qp    = xn{5}(:); 
      end
      %alphaqp = xn{6}(:); 
      %betaqp  = xn{7}(:); 
  
      lqp=length(xn{1}(:));
      %statevariable    = optimvar('state',Ntime,Nspecies,lqp,'LowerBound',0,'UpperBound',.1);
      statevariable  = optimvar('state',Ntime,Nspecies,lqp,'LowerBound',0);
      auxvariable      = optimexpr(   [Ntime,Nspecies,lqp]);
      stateconstraint  = optimconstr(    [Ntime,Nspecies,lqp]);
  
      disp('build state variable')
      
      % >> syms a  kpl d currentTR    T1P kveqp T1L 
      % >> expATR = expm([a,  0; kpl, d ] * currentTR )
      % 
      % expATR =
      % 
      % [                                     exp(a*currentTR),                0]
      % [(kpl*exp(a*currentTR) - kpl*exp(currentTR*d))/(a - d), exp(currentTR*d)]
      % 
      % >> a = -1/T1P - kpl - kveqp
      % >> d = -1/T1L
      % >> eval(expATR)
      % 
      % ans =
      % 
      % [                                                              exp(-currentTR*(kpl + kveqp + 1/T1P)),                   0]
      % [(kpl*exp(-currentTR/T1L) - kpl*exp(-currentTR*(kpl + kveqp + 1/T1P)))/(kpl + kveqp - 1/T1L + 1/T1P), exp(-currentTR/T1L)]
      %    
      %expATR = fcn2optimexpr(@expm,A*currentTR );
      % A = [-1/T1P - kpl - kveqp,  0; kpl, -1/T1L ];
      expATRoneone = exp(-currentTR*(kplqp + kveqp/ve + T1Pqp.^(-1)));
      expATRtwoone = (kplqp.*exp(-currentTR*T1Lqp.^(-1)) - kplqp.*exp(-currentTR*(kplqp + kveqp/ve + T1Pqp.^(-1)))).* (kplqp + kveqp/ve - T1Lqp.^(-1) + T1Pqp.^(-1)).^(-1);
      expATRtwotwo = exp(-currentTR * T1Lqp.^(-1));
      % [expATRoneone(end),0; expATRtwoone(end), expATRtwotwo(end)]
      % IC
      stateconstraint(1,:,:)  = statevariable(1,:,:) ==0;
      auxvariable(1,:,:) =0;
      for iii = 1:Ntime-1
          nsubstep = 5;
          deltat = currentTR /nsubstep ;
          % setup AIF
          integratedt = [TRList(iii):deltat:TRList(iii+1)] +deltat/2  ;
          %integrand = jmA0 * my_gampdf(integratedt(1:nsubstep )'-t0qp,jmalpha,jmbeta) ;
          integrand = jmA0 * gampdf(repmat(integratedt(1:nsubstep )',1,lqp)'- repmat(t0qp,1,nsubstep),jmalpha,jmbeta) ;
          aiftermpyr = deltat * (kveqp/ve).*   exp((- T1Pqp.^(-1) - kplqp - kveqp/ve)*(TRList(iii+1)-deltat*[.5:1:nsubstep]-TRList(iii))) .* integrand ; 
          aiftermlac = deltat * (kveqp/ve).*  ( (-kplqp.*exp((-T1Pqp.^(-1) - kplqp - kveqp/ve)*(TRList(iii+1)-deltat*[.5:1:nsubstep]-TRList(iii)) ) + kplqp.*exp(-T1Lqp.^(-1) *(TRList(iii+1)-deltat*[.5:1:nsubstep]-TRList(iii)))).* ((T1Pqp.^(-1) + kplqp + kveqp/ve) - T1Lqp.^(-1) ).^(-1)   ).* integrand ; 
  
          % setup state as linear constraint
          auxvariable(iii+1,1,:) =  reshape(cos(FaList(1,iii))*expATRoneone.* squeeze( auxvariable(iii,1,: ) ),1,1,lqp ) +  reshape( sum(aiftermpyr,2 ),1,1,lqp) ;
          auxvariable(iii+1,2,:) =  reshape(cos(FaList(2,iii))*expATRtwotwo.* squeeze( auxvariable(iii,2,: ) ),1,1,lqp ) + reshape( sum(aiftermlac,2 ),1,1,lqp) +reshape( cos(FaList(1,iii))*expATRtwoone.* squeeze( auxvariable(iii,1,: )  ),1,1,lqp) ; 
          stateconstraint(iii+1,1,:)  = statevariable(iii+1,1,:) ==  reshape(cos(FaList(1,iii))*expATRoneone.* squeeze( statevariable(iii,1,: ) ),1,1,lqp ) +  reshape( sum(aiftermpyr,2 ),1,1,lqp) ;
          stateconstraint(iii+1,2,:)  = statevariable(iii+1,2,:) ==  reshape(cos(FaList(2,iii))*expATRtwotwo.* squeeze( statevariable(iii,2,: ) ),1,1,lqp ) + reshape( sum(aiftermlac,2 ),1,1,lqp) +reshape( cos(FaList(1,iii))*expATRtwoone.* squeeze( statevariable(iii,1,: )  ),1,1,lqp) ; 
      end
  

      disp('build objective function')
      expandvar  = ones(1,lqp);

      switch (ObjectiveType)
          case('TotalSignal')
            % TODO - repmat does not work well with AD
            % TODO - replace repmat with matrix
            % sumstatevariable = squeeze(sum(repmat(sin(FaList)',1,1,lqp).*statevariable,1));
            sumstatevariable = optimexpr([Nspecies,lqp]);
            for jjj = 1:lqp
               sumstatevariable(:,jjj) =  sum(sin(FaList)'.*(ve*statevariable(:,:,jjj)  + (1-ve) *jmA0  * [gampdf( TRList - t0qp(jjj)  , jmalpha , jmbeta);zeros(1,Ntime)]'  ),1)';

            end 
            diffsumm =(sumstatevariable(1,:)+sumstatevariable(2,:))' * expandvar   - expandvar' * (sumstatevariable(1,:)+sumstatevariable(2,:));
            negHz = 0;
            for jjj=1:lqp2
              znu=xn2{1}(jjj) ;
              % note that the sqrt(pi)^N+1  and 2^N factors from the integration over priors is included in the quadrature routines.
              % this makes is easier to switch between Gaussian and Uniform RV
              negHz = negHz + wn2(jjj) * (wn(:)' * log(exp(-(znu + diffsumm).^2/2/signu^2 - log(signu) -log(2*pi)/2   ) * wn(:)));
            end
          case('SumQuad')
            sumstatevariable = optimexpr([Nspecies,lqp]);
            for jjj = 1:lqp
               sumstatevariable(:,jjj) =  sum(sin(FaList)'.*(ve*statevariable(:,:,jjj)  + (1-ve) *jmA0  * [gampdf( TRList - t0qp(jjj)  , jmalpha , jmbeta);zeros(1,Ntime)]'  ),1)';
            end 
            % max total signal sum
            negHz=-sum(sum(sumstatevariable)) ;
          case('SumTimepoints')
            negHz = 0;
            for jjj=1:lqp2
              znu=xn2{1}(jjj) ;
              for kkk  =1:Ntime 
                 diffsummone = (sin(FaList(1,kkk))*squeeze(statevariable(kkk,1,:)))  * expandvar   - expandvar' * (sin(FaList(1,kkk))*squeeze(statevariable(kkk,1,:))') ;
                 diffsummtwo = (sin(FaList(2,kkk))*squeeze(statevariable(kkk,2,:)))  * expandvar   - expandvar' * (sin(FaList(2,kkk))*squeeze(statevariable(kkk,2,:))') ;
                 negHz = negHz + wn2(jjj) * (wn(:)' * log(exp(-(znu + diffsummone).^2/2/signuImage^2   - (znu + diffsummtwo).^2/2/signuImage^2  ) * wn(:)));
              end
            end
          % https://keplerlounge.com/applied-math/2020/02/13/analytic-min-max.html
          % approximate max
          case('MaxSignal')
            maxstatevariable = optimexpr([Nspecies,lqp]);
            for jjj = 1:lqp
               maxstatevariable(:,jjj) =  mean((sin(FaList)'.*statevariable(:,:,jjj)).^(100),1).^(1/100)';
            end 
            diffsumm =(maxstatevariable(1,:)+maxstatevariable(2,:))' * expandvar   - expandvar' * (maxstatevariable(1,:)+maxstatevariable(2,:));
            negHz = 0;
            for jjj=1:lqp2
              znu=xn2{1}(jjj) ;
              negHz = negHz + wn2(jjj) * (wn(:)' * log(exp(-(znu + diffsumm).^2/2/signu^2 - log(signu) -log(2*pi)/2   ) * wn(:)));
            end
          case('MaxSignalDiff')
            maxstatevariable = optimexpr([Nspecies,lqp]);
            for jjj = 1:lqp
               maxstatevariable(:,jjj) =  mean((sin(FaList)'.*statevariable(:,:,jjj)).^(100),1).^(1/100)';
            end 
            diffsumm =(maxstatevariable(1,:)-maxstatevariable(2,:))' * expandvar   - expandvar' * (maxstatevariable(1,:)-maxstatevariable(2,:));
            negHz = 0;
            for jjj=1:lqp2
              znu=xn2{1}(jjj) ;
              negHz = negHz + wn2(jjj) * (wn(:)' * log(exp(-(znu + diffsumm).^2/2/signu^2 - log(signu) -log(2*pi)/2   ) * wn(:)));
            end
      end
      % MI = H(z) - H(z|P) 
      %  H(z|P)  constant ==> max MI = max H(z) = min -H(z)
      MIGaussObj = negHz;
  
      %% 
      % Create an optimization problem using these converted optimization expressions.
      disp('create optim prob')
      convprob = optimproblem('Objective',MIGaussObj , "Constraints",stateconstraint);
      myidx = varindex(convprob )
      %% 
      % View the new problem.
      %show(convprob)
      problem = prob2struct(convprob,'ObjectiveFunctionName','reducedObjective','ConstraintFunctionName','reducedConstraint');
      

      % compare walker solution at qp
      switch (NumberUncertain)
         case(3)
           optparams = struct('t0',[t0qp(1);0],'gammaPdfA',[alphamean(1)  ;1],'gammaPdfB',[betamean(1);1],...
               'scaleFactor',VIF_scale_fact,'T1s',[T1pmean(1),T1lmean(1)],'ExchangeTerms',[0,kplqp(1) ;0,0],...
               'TRList',TR_list,'PerfusionTerms',[kveqp(1),0],'volumeFractions',ve,...
               'fitOptions', opts);
         case(4)
           error("WIP")
         case(5)
           error("WIP")
      end
      % need uncertainty verification at the same parameters
      refparams = params;
      % solve
      switch (myoptions.Algorithm)
          case('constDirect') 
             % const FA constraints - use direct search instead
             %% Aeq=eye(length(InitialGuess));
             %% Aeq(:,end  )=repmat([0;-1],Ntime,1);
             %% Aeq(:,end-1)=repmat([-1;0],Ntime,1);
             %% Aeq= Aeq(1:end-Nspecies,:); 
             %% beq=zeros(length(InitialGuess)-Nspecies,1);
             [pyrgrid,lacgrid] = meshgrid(0:1:35,0:1:35);
             brutesearch= zeros(size(pyrgrid));
             backspaces = '';
             for iii = 1:length(pyrgrid(:))
                 %disp(sprintf('%d %d ',pyrgrid(iii),lacgrid(iii)));
                 % Print percentage progress
                 percentage = iii/length(pyrgrid(:));
                 perc_str = sprintf('completed %3.1f', percentage);
                 fprintf([backspaces, perc_str]);
                 backspaces = repmat(sprintf('\b'), 1, length(perc_str));

                 x0.FaList = repmat([pi/180*pyrgrid(iii);pi/180*lacgrid(iii)],1,Ntime);
                 x0.state  = evaluate(auxvariable ,x0);
                 brutesearch(iii) = evaluate(MIGaussObj,x0);
             end
             [maxMI,idmax] = max(brutesearch(:));
             [fval,idmin] = min(brutesearch(:));
             popt.FaList = repmat([pi/180*pyrgrid(idmin);pi/180*lacgrid(idmin)],1,Ntime);
             pneg.FaList = repmat([pi/180*pyrgrid(idmax);pi/180*lacgrid(idmax)],1,Ntime);
             optparams.FaList = repmat([pi/180*pyrgrid(idmin);pi/180*lacgrid(idmin)],1,Ntime);
             refparams.FaList = repmat([pi/180*pyrgrid(idmin);pi/180*lacgrid(idmin)],1,Ntime);
  
             handle = figure(5)
             imagesc(brutesearch)
             colorbar
             xlabel('FA Pyruvate - deg')
             ylabel('FA Lactate - deg')
             title(sprintf('max %f min %f',maxMI,fval) )
             text(pyrgrid(idmin)+1,lacgrid(idmin)+1, sprintf('opt %d %d', pyrgrid(idmin), lacgrid(idmin)));
             text(pyrgrid(idmax)+1,lacgrid(idmax)+1, 'control');
          otherwise
             InitialGuess =  [flips(:)];   
             pmin =  [flips(:)*0];     
             pmax =  [flips(:)*0+35*pi/180];
             tolx=1.e-9;
             tolfun=5.e-4;
             maxiter=400;
  
             % truthconstraint = infeasibility(stateconstraint,x0);
             %[popt,fval,exitflag,output] = solve(convprob,x0,'Options',myoptions, 'ConstraintDerivative', 'auto-reverse', 'ObjectiveDerivative', 'auto-reverse' )
             Fx = @(x) MIGHQuadHPTofts(x, problem, myidx,Nspecies,Ntime,auxvariable);
             [designopt,fval,exitflag,output,lambda,grad,hessian] ...
              =fmincon(Fx, InitialGuess ,[],[],[],[],pmin,pmax,[],myoptions);

             handle = figure(5)
             optparams.FaList = reshape(designopt(:),size(params.FaList ));
             refparams.FaList = reshape(designopt(:),size(params.FaList ));
             popt.FaList      = reshape(designopt(:),size(params.FaList ));
      end
      % save convergence history
      set(gca,'FontSize',16)
      saveas(handle,sprintf('historyNG%dNu%d%s%sSNR%02d%s',NGauss,NumberUncertain,myoptions.Algorithm,ObjectiveType,modelSNR,QuadratureRule ),'png')
      popt.state       = evaluate(auxvariable, popt);
      toc;


      [t_axisopt,Mxyopt,Mzopt] = model.compile(M0.',optparams);
      [t_axisref,Mxyref,Mzref] = model.compile(M0.',refparams);
      save(sprintf('poptNG%dNu%d%s%sSNR%02d%s.mat',NGauss,NumberUncertain,myoptions.Algorithm,ObjectiveType,modelSNR,QuadratureRule) ,'fval','popt','params','Mxy','Mz','Mxyref','Mzref','signu','signuImage')
      handle = figure(10)
      plot(params.TRList,Mxyopt(1,:),'b',params.TRList,Mxyopt(2,:),'k')
      ylabel('MI Mxy')
      xlabel('sec'); legend('Pyr','Lac')
      set(gca,'FontSize',16)
      saveas(handle,sprintf('OptMxyNG%dNu%d%s%sSNR%02d%s',NGauss,NumberUncertain,myoptions.Algorithm,ObjectiveType,modelSNR,QuadratureRule),'png')
      handle = figure(11)
      plot(params.TRList,popt.FaList(1,:)*180/pi,'b',params.TRList,popt.FaList(2,:)*180/pi,'k')
      ylabel('MI FA (deg)')
      xlabel('sec'); legend('Pyr','Lac')
      set(gca,'FontSize',16)
      saveas(handle,sprintf('OptFANG%dNu%d%s%sSNR%02d%s',NGauss,NumberUncertain,myoptions.Algorithm,ObjectiveType,modelSNR,QuadratureRule),'png')
      handle = figure(12)
      plot(params.TRList,Mzopt(1,:),'b--',params.TRList,Mzopt(2,:),'k--')
      hold
      plot(params.TRList,cos(optparams.FaList(1,:))'.*popt.state(:,1, 1),'b',params.TRList,cos(optparams.FaList(2,:))'.*popt.state(:,2, 1),'k')
      %if(lqp > 1)
      %  plot(params.TRList,popt.state(:,1, 5),'b',params.TRList,popt.state(:,2, 5),'k')
      %  plot(params.TRList,popt.state(:,1,10),'b',params.TRList,popt.state(:,2,10),'k')
      %  plot(params.TRList,popt.state(:,1,15),'b',params.TRList,popt.state(:,2,15),'k')
      %end
      ylabel('MI Mz ')
      xlabel('sec'); legend('Pyr','Lac')
      set(gca,'FontSize',16)
      saveas(handle,sprintf('OptMzNG%dNu%d%s%sSNR%02d%s',NGauss,NumberUncertain,myoptions.Algorithm,ObjectiveType,modelSNR,QuadratureRule),'png')
  end 


end
function [MIobjfun, MIobjfun_Der]=MIGHQuadHPTofts(xopt,problem,myidx,Nspecies,Ntime,auxvariable)
    x0.FaList = reshape(xopt,Nspecies,Ntime);
    x0.state  = evaluate(auxvariable ,x0);
    Xfull = [ x0.FaList(:); x0.state(:)];
    [MIobjfun,initVals.g] = problem.objective(Xfull);
    [initConst.ineq,initConst.ceq, initConst.ineqGrad,initConst.ceqGrad] = problem.nonlcon(Xfull);
    objectiveGradFA    = initVals.g(myidx.FaList);
    objectiveGradState = initVals.g(myidx.state);
    jacobianFA    = initConst.ceqGrad(myidx.FaList,:);
    jacobianState = initConst.ceqGrad(myidx.state,:);
    adjointvar =-jacobianState \objectiveGradState ;
    MIobjfun_Der = objectiveGradFA +  jacobianFA *   adjointvar ;
end
