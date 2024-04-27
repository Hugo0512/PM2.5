warning off
format long
clear all
clc
result1=[];%存训练适应度
result2=[];%存储测试适应度
result3={};%存储最终函数
result4=[];%存储计算数值和观测值
for qqq=1:10
    qqq
[num,txt,raw]=xlsread('北京');
data=num(:,[5:9,4]);

[row,col]=size(data);
variablenumber=col-1;
selectnum=randperm(size(data,1));
traindata=data(selectnum(1:round(numel(selectnum)*(3/4))),1:variablenumber);
testdata=data(selectnum(round(numel(selectnum)*(3/4))+1:numel(selectnum)),1:variablenumber);
data1=traindata;
observationresult=data(selectnum(1:round(numel(selectnum)*(3/4))),col);%最后一列是目标值
data2=testdata;
observationresult1=data(selectnum(round(numel(selectnum)*(3/4))+1:numel(selectnum)),col);

maxgeneration=100;%进化代数
populationsize=100;%种群规模
participatenum=30;%参加锦标赛的个体数量
genenumber=6;%基因个数
genesize=13;%基因长度
linkfunction='+';%连接函数
mutationrate=0.2;%变异率
oneprecombinationrate=0.3;%单点重组率
twoprecombinationrate=0.3; %双点重组率
generecombinationrate=0.3; %基因重组率
ISelementlength=[1,2,3,4,5];%IS元素长度
IStranspositionrate=0.2;%IS变异率
RISelementlength=[1,2,3,4,5];%RIS元素长度
RIStranspositionrate=0.2;%RIS变异率
genetranspositionrate=0.2;%基因变异率
functionset={'+','-','*','/','S','L','C','I','X','G','A','Q','E','F','D','R','P','~','M','T','O','N','J'};%函数符集合
functionparameter=[2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,1,1];%函数参数个数
for index00=1:variablenumber%终止符集合
    terminalset{index00}=strcat('x',num2str(index00));
end
group={};
headsize=6;%头部长度
tailsize=7;%尾部长度
historybestfitness=-Inf;%记录历史最佳染色体的适应度
historybestindividual=[];%记录历史最佳染色体
bestfitnessarray=[];
averagefitnessarray=[];
group=initial(group,populationsize,genenumber,headsize,genesize,functionset,terminalset,tailsize);%初始化种群
for generation=1:maxgeneration
lengthofORF=[];%存储每个基因的有效长度（ORF）
computationresult=[];
computationresult=computevalue(populationsize,headsize,tailsize,genenumber,group,functionset,lengthofORF,computationresult,genesize,functionparameter,data1,terminalset);%计算数值
fitness=[];
fitness=figurefitness(fitness,observationresult,computationresult,populationsize);%计算适应度
newgroup={};
[newgroup,bestindividual,bestfitness,bestfitnessarray,averagefitnessarray,historybestindividual,historybestfitness]=selection(newgroup,group,fitness,populationsize,bestfitnessarray,averagefitnessarray,historybestindividual,historybestfitness,generation,participatenum);%选择操作
newgroup=geneticoperation(headsize,newgroup,mutationrate,populationsize,genesize,genenumber,ISelementlength,IStranspositionrate,RIStranspositionrate,RISelementlength,genetranspositionrate,generecombinationrate,oneprecombinationrate,twoprecombinationrate,terminalset,functionset,tailsize);%遗传操作
group=newgroup;
end
char(historybestindividual)
computationresult1=[];
fitness1=[];
computationresult1=computevalue(1,headsize,tailsize,genenumber,historybestindividual,functionset,lengthofORF,computationresult1,genesize,functionparameter,data2,terminalset);%计算数值
fitness1=figurefitness(fitness1,observationresult1,computationresult1,1);%计算适应度
a=1:size(observationresult1,1);
plot(a,computationresult1,'-r+',a,observationresult1,'-bx');
xlabel('NO.');
ylabel('value');
legend('计算数值','观测值');
result1(qqq)=historybestfitness;
result2(qqq)=fitness1;
result3(qqq,:)=historybestindividual;
result4(1,:,qqq)=computationresult1;
result4(2,:,qqq)=observationresult1;
end