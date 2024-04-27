%计算适应度
function fitness=figurefitness(fitness,observationresult,computationresult,populationsize)
SSE=[];
SST=[];
average=mean(observationresult);
for index12=1:populationsize
    SSE(index12)=sum((observationresult(:)-computationresult(:,index12)).^2);
    SST(index12)=sum((computationresult(:,index12)-average).^2);
    fitness(index12)=1-SSE(index12)/SST(index12); 
end
for index20=1:populationsize%为不合法的染色体赋予一个很低的适应度（-Inf）
    if isreal(fitness(index20))==0
        fitness(index20)=-Inf;
    end
   if isnan(fitness(index20))==1
       fitness(index20)=-Inf;
   end
    if isinf(fitness(index20))==1
       fitness(index20)=-Inf;
   end
end