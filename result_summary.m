function [RMSE,TIME_mapping,se] = result_summary(sim,obj,params)

tx = sim.tx;
rx = sim.rx;

% process SLAM results
X = nan(4,params.T);
M = cell(1,params.T);
LOS = -ones(1,params.T);
T = zeros(2,params.T);
for k = 1:params.T
    X(:,k) = obj{k}.xn;
    M{1,k} = obj{k}.xl;
    T(:,k) = [obj{k}.dt(1)+obj{k}.dt(2) obj{k}.dt(3)];
    if isempty(obj{k}.los_candidate)
        LOS(k) = 0;
    else
        LOS(k) = 1;
    end
end

% compute performance metrics
K = sum(LOS ~= -1);
N  = [sum(LOS==1) sum(LOS==0) K];
se = [sum((rx(1:2,:) - X(1:2,:)).^2,1); (rx(3:4,:) - X(3:4,:)).^2];

RMSE = zeros(3,3);
for i = 1:3
    RMSE(i,:) = [sqrt(mean(se(i,LOS==1),'omitnan')) sqrt(mean(se(i,LOS==0),'omitnan')) sqrt(mean(se(i,:),'omitnan'))];
end

TIME_localization = [mean(T(1,LOS==1),'omitnan') mean(T(1,LOS==0),'omitnan') mean(T(1,LOS~=-1),'omitnan')];
TIME_mapping = [mean(T(2,LOS==1),'omitnan') mean(T(2,LOS==0),'omitnan') mean(T(2,LOS~=-1),'omitnan')];

c = physconst('LightSpeed');
str = {'LOS','NLOS','ALL'};
for i = 1:3
    fprintf("%4s, estimates: %d/%d, rmse: %.4f [m] %.4f [deg] %.4f [ns], time: %.4f/%.4f [ms]\n",...
        str{i}, N(i),K,RMSE(1,i),RMSE(2,i)*180/pi,RMSE(3,i)*1e9/c,TIME_localization(i)*1000,TIME_mapping(i)*1000);
end


% illustrate SLAM performance
if params.PLOT_ON
    x_lim = [-7 11];
    y_lim = [-12 4];
    
    figure(100); clf; hold on; box on;
    set(gca,'ticklabelinterpreter','latex','fontsize',16)
    axis(gca,[x_lim y_lim])
    xlabel('$x$ [m]','interpreter','latex')
    ylabel('$y$ [m]','interpreter','latex')

    layout = imread('measurements/kampusarena_depleted_background_medium.png');
    [Ny,Nx,~] = size(layout);
    pixel_size_y = 16/Ny;
    pixel_size_x = 18/Nx;
    layout = layout(20:Ny-20,20:Nx-20,:);
    
    x_ = (1:Nx)*pixel_size_x;
    y_ = (1:Ny)*pixel_size_y;
    x_ = x_ - 7;
    y_ = y_ - 12;
    y_ = fliplr(y_);
    im = imagesc(x_,y_,layout);
    im.AlphaData = 0.75;

    plot(tx(1,1),tx(2,1),'bv','markersize',12,'linewidth',2);
    plot(rx(1,:),rx(2,:),'b-','linewidth',2);
    
    mu = cell2mat(M);
    plot(mu(1,:),mu(2,:),'r+','markersize',8,'linewidth',1);
    plot(X(1,:),X(2,:),'kx','markersize',8,'linewidth',2);
    for k = 1:K
        plot([rx(1,k) X(1,k)],[rx(2,k) X(2,k)],'k:','linewidth',1);
    end
end
