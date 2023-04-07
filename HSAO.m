% 天鹰算法——镜像单纯形法+
function [Leader_score,Leader_pos,Convergence_curve,pc,pn] = HSAO(SearchAgents_no,Max_iter,lb,ub,fobj,dim,alpha,delta,omega,u,r0)
N = dim/2;
Fit = zeros(1,SearchAgents_no);% 临时变量
Convergence_curve = zeros(1,Max_iter);
% 次优值
Suboptimal_score = inf;
Suboptimal_pos = zeros(1,dim);
Leader_score = inf;
Leader_pos = zeros(1,dim);
pc = 0;
pn = 0;
% 初始化
Num =1;    % 默认初始化
Positions = initialization(Num,SearchAgents_no, dim, ub, lb);
% 初始化函数值
for i = 1 : SearchAgents_no
    [Fit(i),~,~] = fobj(Positions(i,:));
end
% 找到最优值
for i = 1 : SearchAgents_no 
    % 次优值
    if(Fit(i) < Suboptimal_score)
       Suboptimal_score= Fit(i);
       Suboptimal_pos = Positions(i,:); 
    end
    if(Suboptimal_score < Leader_score)
       Leader_score= Suboptimal_score;
       Leader_pos = Suboptimal_pos; 
    end
end
% 主循环
for iter=1:Max_iter
%     fprintf('Star the 进化 no. %d.\n', iter);
    % 阶段1 最优附近围绕
    if iter< Max_iter*1/3
        x = (r0+u*(1:dim)).*sin(-omega*(1:dim)+3*pi/2);% 公式4
        y = (r0+u*(1:dim)).*cos(-omega*(1:dim)+3*pi/2);% 公式3
        for i = 1:SearchAgents_no
            if rand < 0.5
                Pos_mean = get_mean_pos(dim,SearchAgents_no,Positions);% 当前所有个体的平均位置
                Positions(i,:) = Leader_pos*(1-iter/Max_iter) + rand*(Pos_mean-Leader_pos); % 公式1——信息交流机制
            else
                r_id = randperm(SearchAgents_no,1);
                Positions(i,:) = Leader_pos.*Levy(dim)+Positions(r_id,:)+(y-x)*rand;% 公式2
            end
            % 越界检查，越界后再解空间随机
            Positions(i,:) = Bounds(Positions(i,:), lb, ub,N);
            [Fit(i),~,~] = fobj(Positions(i,:)); 
        end   
    else
        QF=iter^((2*rand()-1)/(1-Max_iter)^2);
        G1 = 2*(1-iter/Max_iter);
        G2 = 2*rand-1;
        for i = 1:SearchAgents_no
            % 将平均扩展为随机步长机制和螺旋交流机制
            % 螺旋因子
            a2 = -1+iter*((-1)/Max_iter);
            b = 1;               	
            l = (a2-1)*rand + 1;   	
            lx = exp(b.*l).*cos(l.*2*pi);
            if rand < 0.5
%                 Pos_mean = get_mean_pos(dim,SearchAgents_no,Positions);
%                 Positions(i,:) = (Leader_pos-Pos_mean)*alpha - rand+delta*(lb+(ub-lb)*rand);% 公式5 
                r_id = randperm(SearchAgents_no,1);
                Positions(i,:) =Positions(i,:)+lx*(Leader_pos-Positions(i,:))+rand*(Positions(r_id,:)-Positions(i,:));
            else
                % 将莱维飞行变为
                Positions(i,:) = QF*Leader_pos-(G2*Positions(i,:)*rand) - G1.*Levy(dim)+rand*G2;% 公式6
            end
            % 越界检查，越界后再解空间随机
            Positions(i,:) = Bounds(Positions(i,:), lb, ub,N);
            [Fit(i),~,~] = fobj(Positions(i,:)); 
        end
    end
    % 找到最优值和次优值
    for i = 1 : SearchAgents_no 
        % 次优值
        if(Fit(i) < Suboptimal_score)
           Suboptimal_score= Fit(i);
           Suboptimal_pos = Positions(i,:); 
        end
        if(Suboptimal_score < Leader_score)
           Leader_score= Suboptimal_score;
           Leader_pos = Suboptimal_pos; 
        end
    end
    % ===============镜像单纯形法=================
    % 反射系数fs,扩张系数kz，收缩系数ss,内收缩系数nss,缩放因子sf
    fs = 1;
    kz = 2;
    ss = 0.5;
    nss = 0.5;
    sf = 1;
    for i = 1 : SearchAgents_no 
        x_c = Suboptimal_pos+Leader_pos;
        x_r = x_c+fs*(x_c-Positions(i,:));
        x_r = Bounds(x_r, lb, ub,N);
        [Fx_r,~,~] = fobj(x_r);
        if Fx_r<Leader_score  % 反向操作成功
            x_e = x_c+kz*(x_r-x_c);
            x_e = Bounds(x_e, lb, ub,N);
            [Fx_e,~,~] = fobj(x_e);
            if Fx_e<Leader_score    % 扩张成功
                Fit(i)=Fx_e;
                Positions(i,:) = x_e;
            else
                x_t = x_c+ss*(Positions(i,:)-x_c);
                x_t = Bounds(x_t, lb, ub,N);
                [Fx_t,~,~] = fobj(x_t);
                if Fx_t<Leader_score    % 外收缩成功
                    Fit(i)=Fx_t;
                    Positions(i,:) = x_t;
                end
            end
        end
        if Fx_r<Fit(i) && Fx_r>Leader_score  % 内缩操作
            x_w = x_c + nss*(Positions(i,:)-x_c);
            x_w = Bounds(x_w, lb, ub,N);
            [Fx_w,~,~] = fobj(x_w);
            if Fx_w<Leader_score    
                Fit(i)=Fx_w;
                Positions(i,:) = x_w;
            end
        end
        x_j = (ub+lb)/2+(ub+lb)/(2*sf)-Positions(i,:)./sf;
        x_j = Bounds(x_j, lb, ub,N);
        [Fx_j,~,~] = fobj(x_j); 
        if Fx_j<Leader_score    
            Fit(i)=Fx_j;
            Positions(i,:) = x_j;
        end
    end
    
    
    % 再次寻找最优值
    for i = 1 : SearchAgents_no 
        if(Fit(i) < Leader_score)
           Leader_score= Fit(i);
           Leader_pos = Positions(i,:); 
        end
    end
    
    [Leader_score,pc,pn] = fobj(Leader_pos);
    Convergence_curve(iter) = Leader_score;
end
% 获取种群平均位置
end

function pos_mean = get_mean_pos(dim,SearchAgents_no,Positions)
    pos_mean = zeros(1,dim);
    for i=1:SearchAgents_no
        pos_mean = pos_mean + Positions(i,:);
    end
    pos_mean = pos_mean/SearchAgents_no;
end

function o=Levy(d)
    beta=1.5;
    sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    u=randn(1,d)*sigma;v=randn(1,d);
    step=u./abs(v).^(1/beta);
    o=step;
end

