% ===========================GAWOA===========================
function [Leader_score,Leader_pos,Convergence_curve,pc,pn] = GLAWOA(SearchAgents_no,Max_iter,lb, ub,fobj,dim)

N = dim;
% 初始化
Leader_pos = zeros(1,dim);
Leader_score = inf; 
Leader_score_pre = Leader_score;
Fit = zeros(1,SearchAgents_no);

% 容差
delta = 1e-6;
todoTol = 0;
Flag = 0;
% 初始化
Num = 0;    % 默认初始化
Positions = initialization(Num,SearchAgents_no, dim, ub, lb);
Convergence_curve = zeros(1,Max_iter);
% 循环计数
iter = 0;
% 主循环
while iter < Max_iter && Flag <= 3
%     w = exp((-10*iter/Max_iter).^2);
    xGroupA = [15 10 6 5 3 2 ];
    rand_lxGroup_index = floor(5*iter/Max_iter+1);
    xGroup = xGroupA(rand_lxGroup_index);
    for i = 1:size(Positions,1)
        Positions(i,:) = Bounds(Positions(i,:), lb, ub,N);
        % 计算每个搜索代理的目标函数
        [Fit(i),~,~] = fobj(Positions(i,:));

        % 更新最优鲸鱼
        if Fit(i) < Leader_score % Change this to > for maximization problem
            Leader_score = Fit(i);         % Update alpha
            Leader_pos = Positions(i,:);
            [~,pc,pn] = fobj(Leader_pos);
        end       
    end
    pFit = Fit; 
    [~,sortIndex] = sort(pFit);             % sortIndex:代理按升序排序
    
    a = 2*cos(pi/2*iter/Max_iter); % 1.非线性收敛因子 
    
    % a2 [-1,-2]线性增加
    a2 = -1+iter*((-1)/Max_iter);
    
    % 更新代理
    for i = 1:size(Positions,1)
        r1 = rand();    % r1 is a random number in [0,1]
        r2 = rand();    % r2 is a random number in [0,1]
        
        A = 2*a*r1-a;   % Eq. (2.3) in the paper
        C = 2*r2;       % Eq. (2.4) in the paper
        
        % parameters for spiral updating position
        b = 1;               	%  parameters in Eq. (2.5)
        l = (a2-1)*rand + 1;   	%  parameters in Eq. (2.5)
        
        p = rand();        		% p in Eq. (2.6)
        %===========================分组==================================
        for j = 1:size(Positions,2)
            % follow the shrinking encircling mechanism or prey search
            if p < 0.5   
                % 搜索
                if abs(A) >= 1
                    rand_leader_index = floor(SearchAgents_no*rand()+1);
                    X_rand = Positions(rand_leader_index, :);
                    D_X_rand = abs(C*X_rand(j) - Positions(i,j)); % Eq. (2.7)
                    Positions(i,j) = X_rand(j) - A*D_X_rand;      % Eq. (2.8)
                % 缩小包围机制  
                elseif abs(A) < 1
                    groupBestPos = groupBest(SearchAgents_no,Positions,sortIndex,ub(j), lb(j),Leader_pos(j),Max_iter,iter,j,xGroup);
%                     D_Leader = abs(C*Leader_pos(j) - Positions(i,j)); % Eq. (2.1)
%                     Positions(i,j) = Leader_pos(j) - A*D_Leader;      % Eq. (2.2)
                    D_Leader = abs(C*groupBestPos - Positions(i,j)); % Eq. (2.1)
                    Positions(i,j) = groupBestPos - A*D_Leader;      % Eq. (2.2)


                end
            % follow the spiral-shaped path
            elseif p>=0.5
                % 螺旋
                distance2Leader = abs(Leader_pos(j)-Positions(i,j));
                Positions(i,j) = distance2Leader*exp(b.*l).*cos(l.*2*pi) + Leader_pos(j);
            end
        end
    end
    % increase the iteration index by 1
    iter = iter + 1;
    Convergence_curve(iter) = Leader_score;
%     [iter Leader_score]
    
    if todoTol == 1 && abs(Leader_score - Leader_score_pre) < delta
        Flag = Flag + 1;
		Convergence_curve = Convergence_curve(1,1:iter);
    end
    Leader_score_pre = Leader_score;
end

end

function groupBestPos = groupBest(SearchAgents_no,Positions,sortIndex,ub, lb,Leader_pos,Max_iter,iter,j,xGroup)
    
    beta=exp((1-iter/Max_iter)*log(1/rand));
    xNum = SearchAgents_no/xGroup;
    X = zeros(1,xGroup);
    for i = 1:xGroup
        rand_leader_index = floor(unifrnd (1+((i-1)*xNum),i*xNum)+1);
        X(i) = Positions(sortIndex(rand_leader_index),j);
    end
    groupBestPos = beta*sum(abs(Leader_pos-X));
    Flag4ub = groupBestPos > ub;
    Flag4lb = groupBestPos < lb;
    groupBestPos = (groupBestPos.*(~(Flag4ub+Flag4lb))) + ub.*Flag4ub + lb.*Flag4lb;
    groupBestPos = floor(groupBestPos);
end



