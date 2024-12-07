function [sim,params,obj] = simulation_setup

    load('measurements/Kampusareena_slam_measurements_svd_p999.mat','sim');

    params.xn_dim = 4;                          % vehicle dimension
    params.xl_dim = 2;                          % landmark dimension
    params.h_dim = 3;                           % measurement dimension
    params.R = ...
        diag([1e-9*physconst('LightSpeed') 3*pi/180 3*pi/180].^2);                      % Covariance of measurements  (range, angle)
    params.T = 45;                              % simulation length
    params.bias = 10;                           % clock bias        
    params.sigma2_bias = 1;                     % process noise of clock
    params.convergence_threshold = 1e-1;        % convergence threshold (Gauss-Newton algorithm)
    params.max_iterations = 10;                 % maximum number of iterations (Gauss-Newton algorithm)
    params.pathloss.theta = [-13 -1.7 1.8^2];   % parameters of log-distance pathloss model
    params.pathloss.threshold = -10.8;          % threshold of LoS detection
    params.PLOT_ON = 1;                         % plot flag
    params.MEX = 1;                             % MEX flag
    params.N = 361;                             % Number of orientation grid points
    params.tx_offset = -1.6*pi/180;             % TX calibration value
    params.epsilon = 0.1;                       % threshold paramter
    params.eta2_threshold = 0.1;                % threshold paramter

    % add clock bias to UE state and calibrate BS
    bias = params.bias;
    sim.rx = vertcat(sim.rx,zeros(1,params.T));
    for k = 1:params.T
        if k < params.T
            bias = bias + sqrt(params.sigma2_bias)*randn;
        end
        sim.rx(4,k) = bias;

        sim.tx(3,k) = sim.tx(3,k) + params.tx_offset;

        sim.y{k}(1,:) = sim.y{k}(1,:) + sim.rx(4,k);
        sim.y{k}(2,:) = mod(sim.y{k}(2,:) + pi,2*pi)-pi;
        sim.y{k}(3,:) = mod(sim.y{k}(3,:) + pi,2*pi)-pi;
    end
    
    obj = cell(1,params.T);