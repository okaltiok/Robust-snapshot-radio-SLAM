classdef slam
    methods (Static)
        function [obj,dt] = localization(tx,y,power,los_candidate,params)
            R = @(theta) [cos(theta) -sin(theta); sin(theta) cos(theta)];
            
            m_k = size(y,2);
            
            gain = 10.^(power/10);
            epsilon = params.epsilon;
            eta2_threshold = params.eta2_threshold;
            
            % determine all possible measurement combinations 
            if isempty(los_candidate)
                combinations = nchoosek(1:m_k,4)';
            else
                combinations = [repmat(los_candidate,1,m_k-1); setdiff(1:m_k,los_candidate)];
            end

            % determine candidate orientations
            if isempty(los_candidate)
                theta = linspace(-pi,pi,params.N)';
            else
                % in LoS conditions orientation is given in closed-form
                v = -R(-y(3,los_candidate))*R(tx(3))*[cos(y(2,los_candidate)); sin(y(2,los_candidate))];
                theta = atan2(v(2),v(1));
            end
            
            % compute snapshot SLAM solution
            time_start = tic;
            if params.MEX
                % MEX implementation of the algorithm
                [LL,X_hat,OUTLIERS] = slam_mex(combinations,y,tx,theta,los_candidate,gain,epsilon,eta2_threshold);
            else
                % loop through all orientations and measurement
                % combinations
                N = size(theta,1);
                M = size(combinations,2);
                LL = nan(M,N);
                X_hat = zeros(3,M,N);
                OUTLIERS = nan(m_k,M,N);
                for n = 1:N
                    [MU,ETA,ETA_bar,HH,AA,bb] = slam.compute_parameters(R,tx,y,theta(n),los_candidate,gain);
                    for m = 1:M
                        outliers = slam.determine_consensus_set(combinations(:,m),y,MU,ETA,ETA_bar,HH,AA,bb,epsilon,eta2_threshold);
                        
                        if outliers(1) ~= -1
                            [LL(m,n),X_hat(:,m,n)] = slam.fit_model(outliers,y,MU,ETA,ETA_bar,HH,AA,bb,gain,epsilon,eta2_threshold);
                            OUTLIERS(:,m,n) = outliers;
                        end
                    end
                end
            end
            dt = toc(time_start);
            
            % estimate UE state
            [i,j] = find(LL == min(LL,[],'all'),1,'first');
            x_hat = [X_hat(1:2,i,j)' theta(j) X_hat(3,i,j)]';
            
            obj.xn = x_hat;
            obj.xl = [];
            obj.L = LL;
            obj.inliers = ~OUTLIERS(:,i,j);
            obj.los_candidate = los_candidate;
        end
        
        
        function [MU,ETA,ETA_bar,HH,AA,bb,UI,VI] = compute_parameters(R,tx,y,theta,los_candidate,gamma)
            % compute parameters used by the LS solution
            I = eye(2);
            m_k = size(y,2);
            
            MU = zeros(2,m_k);
            ETA = zeros(2,m_k);
            ETA_bar = zeros(2,m_k);
            HH = zeros(2,3,m_k);
            AA = zeros(3,3,m_k);
            bb = zeros(3,m_k);
            UI = zeros(2,m_k);
            VI = zeros(2,m_k);
            
            rot_tx = R(tx(3));
            rot_rx = R(theta);
            for j = 1:m_k
                tau = y(1,j);
                aod = y(2,j);
                aoa = y(3,j);
                
                ui = rot_tx*[cos(aod); sin(aod)];
                vi = rot_rx*[cos(aoa); sin(aoa)];
                
                mu = tx(1:2)-tau.*vi;
                H = [I -vi];
                
                if j == los_candidate
                    eta = zeros(2,1);
                    eta_bar = zeros(2,1);
                    
                    A = gamma(j).*(H'*H);
                    b = gamma(j).*(H'*mu);
                else
                    eta = (ui+vi);
                    eta_bar = eta./norm(eta);
                    
                    tmp = H'*(I-eta_bar*eta_bar');
                    A = gamma(j).*(tmp*H);
                    b = gamma(j).*(tmp*mu);
                end
                
                MU(:,j) = mu;
                ETA(:,j) = eta;
                ETA_bar(:,j) = eta_bar;
                HH(:,:,j) = H;
                bb(:,j) = b;
                AA(:,:,j) = A;
                UI(:,j) = ui;
                VI(:,j) = vi;
            end
        end
        
        function outliers = determine_consensus_set(combinations,y,MU,ETA,ETA_bar,HH,AA,bb,epsilon,eta2_threshold)
            % compute UE estimate and determine inliers
            
            % compute model
            x_hat = slam.compute_model(AA,bb,combinations);
            
            if slam.is_feasible(x_hat,y,MU,ETA,HH,combinations,eta2_threshold)
                % compute residuals
                nu2 = slam.compute_cost(x_hat,HH,MU,ETA_bar);

                % determine outliers
                outliers = nu2 > epsilon;

                if sum(~outliers) < size(combinations,1)
                    outliers = -1;
                end
            else
                outliers = -1;
            end
        end
        
        function [L,x_hat] = fit_model(outliers,y,MU,ETA,ETA_bar,HH,AA,bb,gamma,epsilon,eta2_threshold)
            % compute UE estimate and cost of the solution
            idx = find(~outliers);
            
            % compute model
            x_hat = slam.compute_model(AA,bb,idx);
            
            if slam.is_feasible(x_hat,y,MU,ETA,HH,idx,eta2_threshold)
                % compute residuals
                nu2 = slam.compute_cost(x_hat,HH,MU,ETA_bar);
                
                % compute cost
                nu2(outliers == 1) = epsilon;
                L = gamma*nu2;
            else
                L = nan;
            end
        end

        function x_hat = compute_model(AA,bb,idx)
            % estimate UE state, conditional on orientation
            b = zeros(3,1);
            A = zeros(3,3);
            for j = idx'
                b = b + bb(:,j);
                A = A + AA(:,:,j);
            end
            invA = linearAlgebra.matrixInverse(A);
            x_hat = invA*b;
        end
        
        
        function nu2 = compute_cost(x_hat,HH,MU,ETA_bar)
            % compute cost of UE estimate
            m_k = size(HH,3);
            nu2 = zeros(m_k,1);
            for j = 1:m_k
                tmp = HH(:,:,j)*x_hat - MU(:,j);
                nu = tmp - ETA_bar(:,j)'*tmp*ETA_bar(:,j);
                nu2(j) = nu'*nu;
            end
        end
        
        
        function flag = is_feasible(x_hat,y,MU,ETA,HH,idx,eta2_threshold)
            % check if propagation distance is negative
            if sum(y(1,:) - x_hat(3) < 0) > 0
                flag = false;
                return
            end

            [~,los_idx] = min(y(1,:));
            
            flag = true;
            % compute gamma and check if 0 <= gamma <= 1
            for j = idx'
                eta2 = ETA(:,j)'*ETA(:,j);
                if j == los_idx && eta2 < eta2_threshold
                    gamma = 1;
                else
                    gamma = ETA(:,j)'*(HH(:,:,j)*x_hat-MU(:,j)) ...
                        /   ((y(1,j)-x_hat(3))*eta2);
                end

                if gamma < 0 || gamma > 1
                    flag = false;
                    return
                end
            end
        end
        
        
        function [obj,dt] = mapping(obj,tx,y,power,params)
            R = @(theta) [cos(theta) -sin(theta); sin(theta) cos(theta)];
            
            max_iterations = params.max_iterations;
            convergence_threshold = params.convergence_threshold;
            W = linearAlgebra.matrixInverse(params.R);
            
            time_start = tic;
            
            gain = 10.^(power/10);
            los_candidate = obj.los_candidate;
            inliers = obj.inliers;
            xn = obj.xn;
            
            if ~isempty(los_candidate)
                inliers(los_candidate) = 0;
            end
            
            m_k = size(y,2);
            n_k = sum(inliers);
            m = zeros(2,n_k);
            
            [~,los_idx] = min(y(1,:));
            i = 1;
            
            % estimate landmark for every inlier
            for j = 1:m_k
                if inliers(j)
                    % initialize landmark
                    [mu,eta,~,H,~,~,ui,vi] = slam.compute_parameters(R,tx,y(:,j),xn(3),0,gain(j));
                    tau = y(1,j) - xn(4);
                    
                    eta2 = eta'*eta;
                    if j == los_idx && eta2 < params.eta2_threshold
                        gamma = 0.5;
                    else
                        gamma = eta'*(H*xn([1 2 4]) - mu) ...
                            /(tau*eta2);
                    end
                    p1 = tx(1:2)+tau*gamma*ui;
                    p2 = xn(1:2)+tau*(1-gamma)*vi;
                    xl = (p1 + p2)/2;
                    
                    % estimate landmark using the Gauss-Newton algorithm
                    [A,b,L] = slam.landmark_cost(tx,y(:,j),xn,xl,W);
                    for n = 2:max_iterations
                        invA = linearAlgebra.matrixInverse(A);
                        delta = invA*b;
                        
                        [A_up,b_up,L_up] = slam.landmark_cost(tx,y(:,j),xn,xl+delta,W);
                        
                        if L_up >= L
                            xl = xl+delta;
                        else
                            break;
                        end
                        
                        if abs((L - L_up)/L_up) < convergence_threshold
                            break
                        else
                            A = A_up;
                            b = b_up;
                            L = L_up;
                        end
                    end
                    m(:,i) = xl;
                    i = i + 1;
                end
            end
            obj.xl = m;
            dt = toc(time_start);
        end
        
        
        function [A,b,cost] = landmark_cost(tx,y,xn,xl,W)
            % compute cost of landmark estimate
            [mu, ~, Hl] = slam.h_func(tx,xn,xl);
            
            nu = y - mu;
            nu(2) = mod(nu(2) + pi,2*pi)-pi;
            nu(3) = mod(nu(3) + pi,2*pi)-pi;
            
            A = Hl'*W*Hl;
            b = Hl'*W*nu;
            cost = -0.5*nu'*W*nu;
        end
        
        function [mu,Hn,Hl] = h_func(tx,rx,p)
            % geometric measurement model for scattering
            
            % help variables
            sp = p;
            dx1 = tx(1) - sp(1);
            dy1 = tx(2) - sp(2);
            dx2 = sp(1) - rx(1);
            dy2 = sp(2) - rx(2);
            d12 = dx1^2+dy1^2;
            d22 = dx2^2+dy2^2;
            d1 = sqrt(d12);
            d2 = sqrt(d22);
            
            % compute measurement
            mu = [d1+d2 + rx(4); ...                    % TOA
                atan2(-dy1,-dx1) - tx(3); ...   % AOD azimuth
                atan2(dy2,dx2) - rx(3)];        % AOA azimuth
            
            % Compute Jacobians w.r.t. to UE and landmark
            if nargout > 1
                % Compute w.r.t. to UE
                Hn = [-dx2/d2 -dy2/d2 0 1; ...
                    0 0 0 0; ...
                    dy2/d22 -dx2/d22 -1 0];
                
                % Compute w.r.t. to landmark
                Hl = [dx2/d2-dx1/d1 dy2/d2-dy1/d1; ...
                    dy1/d12 -dx1/d12; ...
                    -dy2/d22 dx2/d22];
            end
        end
    end
end