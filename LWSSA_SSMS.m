% ===========麻雀算法===========
function [Leader_score,Leader_pos,Convergence_curve,pc,pn] = LWSSA_SSMS(SearchAgents_no,Max_iter,lb,ub,fobj,dim)

% 划分发现者和跟随者
% 发现者比例为0.2；发现者数量
ST = 0.8;% cos减少[1,0.5]  
P_percent = 0.2; % cos减少[1,0]
pNum = round(SearchAgents_no *  P_percent);  

% 位置，适应度，收敛曲线
Positions = zeros(SearchAgents_no,dim);
Fit = zeros(1,SearchAgents_no);
Convergence_curve = zeros(1,Max_iter);

% 初始化种群
for i = 1 : SearchAgents_no
    Positions(i,:) = lb + (ub - lb) .* rand(1,dim); % x为代理数*dim的矩阵
    [Fit(i),~,~] = fobj(Positions(i,:));                 % 保存函数值       
end
pFit = Fit;                      
pX = Positions;                              % 与pFit对应的个人最佳位置
[Leader_score, bestI] = min(Fit);            % Leader_score表示全局最佳适应值，和当前最佳代理数bestI
Leader_pos = Positions(bestI,:);             % Leader_pos表示Leader_score对应的全局最佳位置

% 主循环
for iter=1:Max_iter
    
% ==================第一步：更新发现者位置_螺旋====================
    [~,sortIndex] = sort(pFit);             % sortIndex:代理按升序排序
    
    [Worst_score,B]=max(pFit);               % fmax表示全局最差适应值，和当前最差代理数B
    Worst_pos= Bounds(Positions(B,:),lb,ub);               % worse表示全局最差位置

    r2=rand(1);
    % 螺旋因子
    a2 = -1+iter*((-1)/Max_iter);
    b = 1;               	
    l = (a2-1)*rand + 1;   	
    z = exp(b.*l).*cos(l.*2*pi);

    % 当没有危险时，麻雀每一维都在缩小
    if(r2<ST)                          % 警戒值为[0.5 1],此刻为0.8
        for i = 1 : pNum                % 方程式（3）
            r1=rand(1);
            Positions(sortIndex(i),:) = Bounds(z*pX(sortIndex(i),:)*exp(-(i)/(r1*Max_iter)),lb,ub);
            [Fit(sortIndex(i)),~,~] = fobj(Positions(sortIndex(i),:));   
        end
    % 当大于警戒值时，发现危险；该发现者将按正太分布随机转移到当前位置附近。（其值收敛于最优位置）    
    else
        for i = 1 : pNum   
            Positions(sortIndex(i),:) = Bounds(pX(sortIndex(i),:)+z*randn(1)*ones(1,dim),lb,ub);
            [Fit(sortIndex(i)),~,~] = fobj(Positions(sortIndex(i),:));
        end   
    end
% ==================第二步：更新跟随者位置_莱维飞行====================
    [~, bestII] = min(Fit);      % fMMin表示全局最佳适应值，和当前最佳代理数bestII
    bestXX = Positions(bestII,:);            

    for i = (pNum+1) : SearchAgents_no                     % 方程式（4）
        A=floor(rand(1,dim)*2)*2-1;
        % 种群收敛时，符合标准正态分布随机数
        if(i>(SearchAgents_no/4))
            Positions(sortIndex(i),:)=randn(1)*exp((Worst_pos-pX(sortIndex( i ),:))/(i)^2);
        % 种群收敛时，其值通过莱维飞行的扰动收敛于最优位置
        else
            Positions(sortIndex(i),:)=bestXX+(abs(( pX( sortIndex(i),:)-bestXX)))*(A'/dim)*ones(1,dim).*Levy(dim);  
        end  
        Positions( sortIndex(i),:) = Bounds(Positions(sortIndex(i),:), lb, ub);
        [Fit(sortIndex(i)),~,~] = fobj( Positions(sortIndex(i),:));
    end

% 当发现危险时，所有麻雀放弃当前食物而转移到新的位置
    c=randperm(numel(sortIndex));
    b=sortIndex(c(1:20));
    for j =  1  : length(b)      % 方程式（5）
        % 如果不是最优，那么其值收敛于最优
        if(pFit( sortIndex( b(j) ) )>(Leader_score))
            Positions( sortIndex( b(j) ), : )=Leader_pos+(randn(1,dim)).*(abs(( pX( sortIndex( b(j) ), : ) -Leader_pos)));
        % 如果是最优，那么其值会远离最优
        else
            Positions( sortIndex( b(j) ), : ) =pX( sortIndex( b(j) ), : )+(2*rand(1)-1)*(abs(pX( sortIndex( b(j) ), : )-Worst_pos))/ ( pFit( sortIndex( b(j) ) )-Worst_score+1e-50);
        end
        Positions( sortIndex(b(j) ), : ) = Bounds(Positions( sortIndex(b(j) ), : ),lb,ub);
        [Fit( sortIndex(b(j))),~,~]= fobj( Positions( sortIndex(b(j)),:));
    end
    % 找到最优值
    for i = 1 : SearchAgents_no 
        if ( Fit(i) < pFit(i))
            pFit(i) = Fit(i);
            pX(i,:) = Positions(i,:);
        end
        if(pFit(i) < Leader_score)
           Leader_score= pFit(i);
           Leader_pos = pX(i,:); 
        end
    end
    [~,pc,pn]= fobj(Leader_pos);
    Convergence_curve(iter) = Leader_score;
end
end

function o=Levy(d)
    beta=1.5;% 指数1.5
    sigma=(gamma(1+beta)*sin(pi*beta/2)/(gamma((1+beta)/2)*beta*2^((beta-1)/2)))^(1/beta);
    u=randn(1,d)*sigma;v=randn(1,d);%服从正太分布
    step=u./abs(v).^(1/beta);% 飞行步长
    o=step;
end
% 麻雀算法会收敛于原点和当前最优位置
function s = Bounds(s, lb, ub)
    % 返回超出搜索空间边界的搜索代理
        Flag4ub = s > ub;
        Flag4lb = s < lb;
        s = (s.*(~(Flag4ub+Flag4lb))) + (ub-0.25.*(ub-lb).*rand).*Flag4ub + (lb+0.25.*(ub-lb).*rand).*Flag4lb;
end
% 个人建议：1.去除向原点收敛的操作；2.减少向最优位置的跳跃，换成向最优位置移动。
