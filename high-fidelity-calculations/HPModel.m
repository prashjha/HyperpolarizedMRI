classdef HPModel
    % HPModel - A chemical exchange model assuming two pooled
    % Tofts model of perfusion
    % Parameters:
    %* ExchangeTerms - A Matrix defining chemical Exchange
    %* T1s - A row vector of T1 decay terms
    %* FaList - A matrix of excitation angles
    %* TRList - A matrix of excitation times
    %* t0 - A row vector for delivery delay of each metabolite
    %* gammaPdfA - A row vector for shape term alpha of each metabolite
    %* gammaPdfB - A row vector for shape term beta of each metabolite
    %* ScaleFactor - A row vector for each metabolite's VIF scale factor
    %* PerfusionTerms - A row vector for each metabolite's extravisation rate
    %* volumeFractions - A row vector for each metabolite's volume fraction
    properties
    end
    methods
        function params = defaultParams(self, Nscans)
            % parameters based on discussion Dec 7, 2021
            % parameters based on updated calculations by David (May 2022)
            N = 30;
            if nargin > 1
                N = Nscans; 
            end
            TR = 3; TRList = (0:(N-1))*TR; TRi = TR + (1:(N-1))*0;
            FaList = zeros(2, N);
            for i=1:N 
                FaList(1,i) = 20*pi/180.; FaList(2,i) = 30*pi/180.;
            end
            VIF_t0 = [4; 0]; gammaPdfA = 2.5; gammaPdfB = 4.5; VIF_scale = [100; 0];
            VIF_fn = @(t) VIF_scale .* gampdf(t - VIF_t0, gammaPdfA, gammaPdfB);
            kve = 0.05; ve = 0.95; T1s = [30; 25]; 
            kpl = 0.15; klp = 0;

            % set values
            params = struct( 'TRi', TRi, 'TRList', TRList, ...
                             'FaList', FaList, ...
                             't0', VIF_t0, 'gammaPdfA', [gammaPdfA; 1], ...
                             'gammaPdfB', [gammaPdfB; 1], ...
                             'scaleFactor', VIF_scale, ...
                             'T1s', T1s, ...
                             'ExchangeTerms', [0, kpl; klp, 0], ...
                             'kve', [kve; 0], ...
                             've', [ve; ve], ...
                             'VIF', VIF_fn);
        end

        function paramsOut = parseParams(self,paramsIn)
            % copy
            paramsOut = paramsIn;

            % fill missing data
            if ~isfield(paramsOut, 'N')
                def_params = self.defaultParams();
            else 
                def_params = self.defaultParams(paramsOut.N)
            end
            tmpNames = fieldnames(def_params);
            for i=1:numel(tmpNames)
                if ~isfield(paramsOut, tmpNames{i})
                    paramsOut.(tmpNames{i}) = def_params.(tmpNames{i})
                end
            end

            % update secondary parameters in case model or design parameters
            % have been changed
            N = size(paramsIn.ExchangeTerms,1);
            K = triu(paramsIn.ExchangeTerms)+tril(paramsIn.ExchangeTerms);
            T1 = paramsIn.T1s;
            kve = paramsOut.kve;
            ve = paramsOut.ve;

            % matrix A
            A = zeros(N);
            for i = 1:N
                for j = 1:N
                    if(i == j)
                        A(i,i) = -sum(K(i,:))-1/T1(i) - kve(i)/ve(i);;
                    else
                        % This transposes A to match conventional matrix
                        % multiplication convintions
                        A(i,j) = K(j,i);
                    end
                end
            end
            paramsOut.A = A;

            % Build VIF
            paramsOut.VIF = @(t)paramsOut.scaleFactor.*...
                gampdf(t-paramsOut.t0,paramsOut.gammaPdfA,paramsOut.gammaPdfB);
        end

        % convert TR into time sequence
        function TR = getTR(self, t)
            N = size(t,2);
            TR = zeros(1,N);
            TR(1) = t(1);
            for i=2:N
                TR(i) = t(i) - t(i-1);
            end
            TR = TR(2:end);
        end

        % convert time into TR sequence
        function t = getTime(self, TR)
        % compute time sequence from TR
            N = size(TR,2) + 1;
            t = zeros(1,N);
            for i=2:N
                t(i) = 0;
                for j=2:i
                    t(i) = t(i) + TR(j-1);
                end
            end
        end

        function [A] = getA(self,params)
            params = self.parseParams(params);
            A = params.A;
        end

        function [VIF] = getVIF(self,params)
            params = self.parseParams(params);
            VIF_fn = params.VIF;
            FaList = params.FaList;
            TRList = params.TRList;
            N = size(TRList,2);
            
            VIF = zeros(size(FaList));
            for i = 1:N
                VIF(:,i) = VIF_fn(TRList(i));
            end
        end
        function [TRList,Mxy,Mz] = compile(self,M0,params)
            params = self.parseParams(params);
            [TRList, Mxy, Mz] = self.evaluate(M0,params);
        end
        
        function [TRList,Mxy,Mz,G,dGdTR,dGdFA] = compile_der(self,M0,params)
            params = self.parseParams(params);
            [TRList, Mxy, Mz, G, dGdTR, dGdFA] = self.evaluate_der(M0,params);
        end

        function [TRList,Mxy,Mz,G,dGdTR,dGdFA] = compile_der_const_design(self,M0,params)
            params = self.parseParams(params);
            [TRList, Mxy, Mz, G, dGdTR_var, dGdFA_var] = self.evaluate_der(M0,params);

            dGdTR = zeros(size(params.TRList));
            dGdFA = zeros(size(params.FaList));

            % dGdTR is scalar value
            dGdTR(1,1) = sum(dGdTR_var);

            % dGdFA is vector of two values
            dGdFA = zeros(2,1);
            dGdFA(1,1) = sum(dGdFA_var(1, :));
            dGdFA(2,1) = sum(dGdFA_var(2, :));
        end
    end
        
    methods (Access = private)
         function [TRList, Mxy, Mz] = evaluate(self,M0,params)
            A = params.A;
            TRList = params.TRList;
            FaList = params.FaList;
            VIF_fn = params.VIF;
            kve = params.kve;
            ve = params.ve;
            N = length(TRList);
            
            Mz = zeros(size(FaList));
            Mxy = zeros(size(FaList));
            
            % first step
            Mz(:,1) = M0.*cos(FaList(:,1));
            Mxy(:,1) = (ve.*M0+(1-ve).*VIF_fn(TRList(1))).*sin(FaList(:,1));
            
            %expATR = expm((TRList(2) - TRList(1)) * A);
            
            sel_method = 0; 
            if sel_method == 0
                fun = @(t,y)A*y+(kve./ve).*VIF_fn(t);
                for i = 2:N
                    [~,Y] = ode45(fun,[TRList(i-1),TRList(i)],Mz(:,i-1));
                    Mz(:,i) = Y(end,:)';
                    Mxy(:,i) = sin(FaList(:,i)).*(ve.*Mz(:,i)+...
                        (1-ve).*VIF_fn(TRList(i)));
                    Mz(:,i) = cos(FaList(:,i)).*Mz(:,i);
                end
            elseif sel_method == 1
                % compute TRs from TRList
                TR = zeros(N,1);
                TR(1) = TRList(1);
                for i = 2:N
                    TR(i) = TRList(i) - TRList(i-1);
                end
                
                fun = @(t,tk) expm((tk - t)*A)*((kve./ve).*VIF_fn(t));
                for i = 2:N
                    ti = TRList(i);
                    tim = TRList(i-1);
                    
                    Mz(:,i) = expm(TR(i)*A)*Mz(:,i-1) + integral(@(x) fun(x,ti), tim, ti,'ArrayValued',true,'RelTol',1e-7);
                    Mxy(:,i) = sin(FaList(:,i)).*(ve.*Mz(:,i)+...
                        (1-ve).*VIF_fn(TRList(i)));
                    Mz(:,i) = cos(FaList(:,i)).*Mz(:,i);
                end
            end
         end
         
         function [TRList, Mxy, Mz, G, dGdTR, dGdFA] = evaluate_der(self,M0,params)
            A = params.A;
            TRList = params.TRList;
            FaList = params.FaList;
            VIF_fn = params.VIF;
            kve = params.kve;
            ve = params.ve;
            N = length(TRList);
            gpdfa = params.gammaPdfA(1); % for pyruvate
            gpdfb = params.gammaPdfB(1); % for pyruvate
            
            % compute TRs from TRList
            TR = zeros(N,1);
            TR(1) = TRList(1);
            for i = 2:N
                TR(i) = TRList(i) - TRList(i-1);
            end
            
            Mz = zeros(size(FaList));
            Mz_nocos = zeros(size(FaList));  % without post cosine multiplication for derivative calc
            Mxy = zeros(size(FaList));
            
            % first step
            Mz(:,1) = M0.*cos(FaList(:,1));
            Mz_nocos(:,1) = M0;
            Mxy(:,1) = (ve.*M0+(1-ve).*VIF_fn(TRList(1))).*sin(FaList(:,1));
            
            sel_method = 0; 
            if sel_method == 0
                fun = @(t,y)A*y+(kve./ve).*VIF_fn(t);
                for i = 2:N
                    [~,Y] = ode45(fun,[TRList(i-1),TRList(i)],Mz(:,i-1));
                    Mz(:,i) = Y(end,:)';
                    Mxy(:,i) = sin(FaList(:,i)).*(ve.*Mz(:,i)+(1-ve).*VIF_fn(TRList(i)));
                    Mz_nocos(:,i) = Mz(:,i);
                    Mz(:,i) = cos(FaList(:,i)).*Mz(:,i);
                end
            elseif sel_method == 1
                fun = @(t,tk) expm((tk - t)*A)*((kve./ve).*VIF_fn(t));
                for i = 2:N
                    ti = TRList(i);
                    tim = TRList(i-1);
                    
                    Mz(:,i) = expm(TR(i)*A)*Mz(:,i-1) + integral(@(x) fun(x,ti), tim, ti,'ArrayValued',true,'RelTol',1e-7);
                    Mxy(:,i) = sin(FaList(:,i)).*(ve.*Mz(:,i)+(1-ve).*VIF_fn(TRList(i)));
                    Mz_nocos(:,i) = Mz(:,i);
                    Mz(:,i) = cos(FaList(:,i)).*Mz(:,i);
                end
            end
            
            % adjoint step (adjoint is vector of size N+1, where last 
            % element of vector is zero
            Q = zeros(size(FaList,1), N+1); % adjoint
            Q(:,N+1) = zeros(2,1); % init adjoint
            for i = N:-1:1
                fi = FaList(:,i);
                TRip = 0;
                if i < N
                    TRip = TR(i+1);
                end
                Q(:,i) = cos(fi).*(expm(TRip*A).'*Q(:,i+1)) + ve.*sin(fi);
            end
            
            % signal function and derivatives
            G =  sum(Mxy(:));
            dGdTR = zeros(size(TRList));
            dGdFA = zeros(size(FaList));
            
            % loop for derivatives
            for i = 1:N
                ti = TRList(i);
                TRi = TR(i);
                TRip = 0;
                if i < N
                    TRip = TR(i+1);
                end
                Pi = FaList(1,i);
                Li = FaList(2,i);
                Atri = expm(TRi*A);
                Atrip = expm(TRip*A);
                
                % contribution from der of sin term
                xx = (1-ve).*VIF_fn(ti) + ve.*Mz_nocos(:,i);
                
                ds = zeros(2,1);
                ds(1) = cos(Pi);
                dGdFA(1,i) = dot(ds,xx);
                
                ds(1) = 0;
                ds(2) = cos(Li);
                dGdFA(2,i) = dot(ds,xx);
                
                % contribution from der of VIF term
                if i > 1
                    for j=i:N
                        tj = TRList(j);
                        ds = sin(FaList(:,j));
                        ds(2) = 0; % so that VIF for lactate is ignored
                        % derivative of gamma fn: dt g(t,a,b) =
                        % g(t,a,b)[(a-1)/t - 1/b]
                        dGdTR(:,i) = dGdTR(:,i) + ((gpdfa - 1)/tj - 1/gpdfb)*dot(ds,(1-ve).*VIF_fn(tj));
                    end
                end
                
                % contribution from der of fwd model to dGdFA
                if i < N
                    ds = zeros(2,1);
                    ds(1) = -sin(Pi);
                    dGdFA(1,i) = dGdFA(1,i) + dot(Q(:,i+1), Atrip*(ds.*Mz_nocos(:,i)));
                    
                    ds(1) = 0;
                    ds(2) = -sin(Li);
                    dGdFA(2,i) = dGdFA(2,i) + dot(Q(:,i+1), Atrip*(ds.*Mz_nocos(:,i)));
                end
                
                % contribution to dGdTR
                if i > 1
                    ds = cos(FaList(:,i-1)).*Mz_nocos(:,i-1);
                    dGdTR(:,i) = dGdTR(:,i) + dot(Q(:,i), A*Atri*ds);
               

                    % lastly, add contribution to dGdTR from der VIF term
                    for j=i:N
                        % compute derivative of \int_{t_{i-1}}^{t_i} \exp[A(t_i
                        % - \tau)] VIF[\tau] d\tau
                        tj = TRList(j);
                        tjm = TRList(j-1);

                        % term d t_j/d TR_i (from derivative of upper 
                        % integral limit in Leibniz rule)
                        db = (kve./ve).*VIF_fn(tj);

                        % term d t_{j-1} / d TR_i
                        if j-1 >= i
                            db = db - expm((tj-tjm)*A)*((kve./ve).*VIF_fn(tjm));
                        end

                        % from integral
                        % compute integral from solution instead of integrating again
                        bj = Mz_nocos(:,j) - expm(TR(j)*A)*(cos(FaList(:,j-1)).*Mz_nocos(:,j-1));
                        % or compute by integration
                        %bj = (kve.'/ve).*integral(@(x) fun(x,tj), tjm, tj,'ArrayValued',true);
                        db = db + A*bj;

                        dGdTR(:,i) = dGdTR(:,i) + dot(Q(:,j), db);
                    end
                end
            end
         end
    end
end
    
