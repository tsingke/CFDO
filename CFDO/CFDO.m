

function [gbestx,gbestfitness,BestCost1] = CFDO(~,popsize,dimension,xmax,xmin,~,~,maxiter,Func,FuncId,~)
%% Problem Definition
CostFunction=Func;
nVar=dimension;                 % Number of Decision Variables
VarSize=[1 nVar];       % Decision Variables Matrix Size
VarMin=xmin;             % Decision Variables Lower Bound
VarMax= xmax;         % Decision Variables Upper Bound

%%  Parameters
FEs = 0;
MaxFEs = nVar*10000;
rank1 = round(0.2*popsize);
rank2 = round(0.4*popsize);
rank3 = popsize-rank1-rank2;
pxmin=xmin.*rand(1,dimension);
pxmax=xmax.*rand(1,dimension);
dist2 = sum(abs(pxmin - pxmax));% 计算曼哈顿距离

nPop = popsize;         % Population Size
%% Initialization
empty_plant.Position = [];
empty_plant.Cost = [];
empty_plant.FitDis  = [];

pop = repmat(empty_plant, nPop, 1);    % Initial Population Array
% Initialize Best Solution Ever Found
BestSol.Cost=inf;

for i = 1:numel(pop)

    % Initialize Position
    pop(i).Position = unifrnd(VarMin, VarMax, VarSize);
    % Evaluation
    pop(i).Cost = CostFunction(pop(i).Position',FuncId);
    FEs = FEs+1;


    if pop(i).Cost<=BestSol.Cost
        BestSol.Position=pop(i).Position;
        BestSol.Cost=pop(i).Cost;
    end
    BestCost1(FEs)=BestSol.Cost;
end

%%  Main Loop

% FEs=nPop;
% % for FEs=1:MaxFEs
while FEs<=MaxFEs
    % Get Best and Worst Cost Value

    [~, SortOrder]=sort([pop.Cost]);%排序
    pop = pop(SortOrder);
    Costs = [pop.Cost];

    BestCost = min(Costs);
    WorstCost = max(Costs);

    rank3= (0.5 - 0.005 * 10 ^ (2 * FEs / MaxFEs))*popsize;
    rank2=popsize-rank3;

    bestp=pop(1).Position;
    medinp=pop(floor(popsize/2)).Position;

    dist1  = sum(abs(bestp - medinp)); % 计算曼哈顿距离
    d=dist1/dist2;

    for i = 1 : popsize
        value = 0;
        for j = 1 : dimension
            value = value + abs(bestp(j) -pop(i).Position(j));
        end
        distances(i) = value;

    end

    MaxMinFitness = WorstCost - BestCost;
    MinDistance = min(distances);
    medDistance =median(distances);
    MaxMinDistance = max(distances) - MinDistance;
    if d > rand
        d=1-d;
    end

    for i = 1 : popsize
        Fitness(i) = 1 - ((Costs(i) - BestCost) / MaxMinFitness);
        Distances(i) = (distances(i) - MinDistance) / MaxMinDistance;

        divDistances(i) = d*Fitness(i) +(1-d)*Distances(i);
        pop(i).FitDis =  divDistances(i);
    end
    FitDiss = [pop.FitDis];

    rlist = randperm(popsize);
    rpairs = [rlist(1:ceil(popsize/2)); rlist(floor(popsize/2) + 1:popsize)]';
    % 粒子对竞争
    mask = (Costs(rpairs(:,1))>Costs(rpairs(:,2)));
    for k = 1:ceil(popsize/2)
        newsol = empty_plant;
        if mask(k)==0
            los=rpairs(k,2);
            win=rpairs(k,1);
        else
            los=rpairs(k,1);
            win=rpairs(k,2);
        end
    
        selected_indices=selection(Costs,popsize,los);
        rp=selected_indices(popsize-2:popsize);
        a3=rp(1);
        a2=rp(2);
        a1=rp(3);
%         AA1=((pop(a2).Cost-pop(a3).Cost)/abs(pop(a3).Cost-pop(a2).Cost));
%         if (pop(a2).Cost-pop(a3).Cost)==0
%             AA1=1;
%         end
%         A2=((pop(a1).Cost-pop(los).Cost)/abs(pop(a1).Cost-pop(los).Cost));
%         if (pop(a1).Cost-pop(los).Cost)==0
%             AA2=1;
%         end

            if rand<0.5 %% 
                A1=(pop(a3).FitDis-pop(a1).FitDis)/(pop(a2).FitDis-pop(a1).FitDis);
                A2=-A1;




            else
                A2=(pop(a3).FitDis-pop(a1).FitDis)/(pop(a2).FitDis-pop(a1).FitDis);

                A1=-A2;
            end

        Pos= pop(los).Position+(A1*pop(win).FitDis).*(pop(los).Position-pop(a1).Position)+(A2*pop(win).FitDis).*(pop(a3).Position-pop(a2).Position);

        FEs = FEs+1;

        newsol.Position =Pos;
        newsol.Position = max(newsol.Position, VarMin);
        newsol.Position = min(newsol.Position, VarMax);
        % Evaluate
        newsol.Cost = CostFunction(newsol.Position',FuncId);

        if newsol.Cost<pop(los).Cost
            pop(los).Position=newsol.Position;
            pop(los).Cost=newsol.Cost;

            if pop(los).Cost<=BestSol.Cost
                BestSol=pop(los);
            end

        end

        BestCost1(FEs)=BestSol.Cost;
%      ratioi= (pop(win).Cost - WorstCost)/(BestCost - WorstCost);

     r=(pop(win).Cost- BestCost) / MaxMinFitness;
     M=rand.*(1-r);  
     ratioi=M*d;
   
        newsol = empty_plant;
        
        selected_indices=selection(Costs,popsize,i);
        rp=selected_indices(popsize-2:popsize);
        a3=rp(1);
        a2=rp(2);
        a1=rp(3);

        sigma1=ratioi;
        newsol.Position =pop(a1).Position+(sigma1).*(pop(a3).Position-pop(a2).Position);
        see=pop(win).Position;
        
        for j =1:nVar
            if rand>0.9
                Pos1(1,j) =newsol.Position(1,j);
            else
                Pos1(1,j)= see(1,j);
            end
        end

        newsol.Position =Pos1;
        % Apply Lower/Upper Bounds
        newsol.Position = max(newsol.Position, VarMin);
        newsol.Position = min(newsol.Position, VarMax);
        
        % Evaluate
        newsol.Cost = CostFunction(newsol.Position',FuncId);  
        FEs=FEs+1;
        if newsol.Cost<pop(win).Cost
            pop(win).Position=newsol.Position;
            pop(win).Cost=newsol.Cost;

            if pop(win).Cost<=BestSol.Cost
                BestSol=pop(win);
            end
        end

     
        BestCost1(FEs)=BestSol.Cost;








        if mod(FEs,MaxFEs/10) == 0 && FEs <= MaxFEs
            fprintf("BCB 第%d次评价，最佳适应度 = %e\n",FEs,BestSol.Cost);
        end
    end
end

gbestx = BestSol.Position;
gbestfitness = BestSol.Cost;



if FEs<MaxFEs
    BestCost1(FEs+1:MaxFEs)=gbestfitness;
else
    if FEs>MaxFEs
        BestCost1(MaxFEs+1:end)=[];
    end
end
end
function selected_indices = selection(fitness, popsize,i)
fitness(i)=0;
    % 计算适应度总和
    total_fitness = sum(fitness);
    
    % 计算每个个体的相对适应度
    relative_fitness = fitness / total_fitness;
    
    % 计算轮盘上的累积概率
    cumulative_probability = cumsum(relative_fitness);
    
    % 初始化所选个体的索引数组
    selected_indices = zeros(1, popsize);
    
    % 执行轮盘赌选择操作
    for i = 1:popsize
        % 生成一个随机数，范围在0到1之间
        random_value = rand();
        
        % 找到第一个累积概率大于随机数的个体
        j = 1;
        while cumulative_probability(j) < random_value
            j = j + 1;
        end
        
        % 选中该个体
        selected_indices(i) = j;
    end
end

