%ѡ�������ʹ�ý����������Ӿ�Ӣ����
function [newgroup,bestindividual,bestfitness,bestfitnessarray,averagefitnessarray,historybestindividual,historybestfitness]=selection(newgroup,group,fitness,populationsize,bestfitnessarray,averagefitnessarray,historybestindividual,historybestfitness,generation,participatenum)
bestfitness=fitness(1);
bestindividual=group(1,:);%��ʱ�����Ⱦɫ������Ϊ��һ��Ⱦɫ��
bestposition=1;%��ʱ�����Ⱦɫ��ĺ�������Ϊ1
for index22=1:populationsize
    if fitness(index22)>bestfitness
        bestposition=index22;
        bestfitness=fitness(index22);%��¼��ǰ�����Ⱦɫ�����Ӧ��
    end
end
bestindividual=group(bestposition,:);%��¼��ǰ�����Ⱦɫ��
if bestfitness>=historybestfitness
    historybestindividual=bestindividual;
    historybestfitness=bestfitness;
end
historybestfitness
newgroup(1,:)=historybestindividual;
bestfitnessarray(generation)=historybestfitness;
averagefitnessarray(generation)=mean(fitness);
%���ý�����������ѡ��(ѡ�����������Ⱥ��ģ��1����Ϊ��Ӣ����)
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