%选择操作，使用锦标赛方法加精英保留
function [newgroup,bestindividual,bestfitness,bestfitnessarray,averagefitnessarray,historybestindividual,historybestfitness]=selection(newgroup,group,fitness,populationsize,bestfitnessarray,averagefitnessarray,historybestindividual,historybestfitness,generation,participatenum)
bestfitness=fitness(1);
bestindividual=group(1,:);%暂时将最佳染色体设置为第一个染色体
bestposition=1;%暂时将最佳染色体的号码设置为1
for index22=1:populationsize
    if fitness(index22)>bestfitness
        bestposition=index22;
        bestfitness=fitness(index22);%记录当前代最佳染色体的适应度
    end
end
bestindividual=group(bestposition,:);%记录当前代最佳染色体
if bestfitness>=historybestfitness
    historybestindividual=bestindividual;
    historybestfitness=bestfitness;
end
historybestfitness
newgroup(1,:)=historybestindividual;
bestfitnessarray(generation)=historybestfitness;
averagefitnessarray(generation)=mean(fitness);
%采用锦标赛法进行选择(选择次数少于种群规模数1，因为精英保留)
for index21=2:populationsize
    sortedindividual=randperm(populationsize);
    tempindividual=sortedindividual(1:participatenum);
    comparemax=fitness(tempindividual(1));
    maxposition=tempindividual(1);
   for index23=1:numel(tempindividual)
     if   fitness(tempindividual(index23))>comparemax
         maxposition=tempindividual(index23); 
     end
   end
   newgroup(index21,:)=group(maxposition,:);
end