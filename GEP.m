warning off
format long
clear all
clc
result1=[];%��ѵ����Ӧ��
result2=[];%�洢������Ӧ��
result3={};%�洢���պ���
result4=[];%�洢������ֵ�͹۲�ֵ
for qqq=1:10
    qqq
[num,txt,raw]=xlsread('����');
data=num(:,[5:9,4]);

[row,col]=size(data);
variablenumber=col-1;
selectnum=randperm(size(data,1));
traindata=data(selectnum(1:round(numel(selectnum)*(3/4))),1:variablenumber);
testdata=data(selectnum(round(numel(selectnum)*(3/4))+1:numel(selectnum)),1:variablenumber);
data1=traindata;
observationresult=data(selectnum(1:round(numel(selectnum)*(3/4))),col);%���һ����Ŀ��ֵ
data2=testdata;
observationresult1=data(selectnum(round(numel(selectnum)*(3/4))+1:numel(selectnum)),col);

maxgeneration=100;%��������
populationsize=100;%��Ⱥ��ģ
participatenum=30;%�μӽ������ĸ�������
genenumber=6;%�������
genesize=13;%���򳤶�
linkfunction='+';%���Ӻ���
mutationrate=0.2;%������
oneprecombinationrate=0.3;%����������
twoprecombinationrate=0.3; %˫��������
generecombinationrate=0.3; %����������
ISelementlength=[1,2,3,4,5];%ISԪ�س���
IStranspositionrate=0.2;%IS������
RISelementlength=[1,2,3,4,5];%RISԪ�س���
RIStranspositionrate=0.2;%RIS������
genetranspositionrate=0.2;%���������
functionset={'+','-','*','/','S','L','C','I','X','G','A','Q','E','F','D','R','P','~','M','T','O','N','J'};%����������
functionparameter=[2,2,2,2,1,1,1,1,1,1,1,1,1,1,1,1,1,1,2,2,2,1,1];%������������
for index00=1:variablenumber%��ֹ������
    terminalset{index00}=strcat('x',num2str(index00));
end
group={};
headsize=6;%ͷ������
tailsize=7;%β������
historybestfitness=-Inf;%��¼��ʷ���Ⱦɫ�����Ӧ��
historybestindividual=[];%��¼��ʷ���Ⱦɫ��
bestfitnessarray=[];
averagefitnessarray=[];
group=initial(group,populationsize,genenumber,headsize,genesize,functionset,terminalset,tailsize);%��ʼ����Ⱥ
for generation=1:maxgeneration
lengthofORF=[];%�洢ÿ���������Ч���ȣ�ORF��
computationresult=[];
computationresult=computevalue(populationsize,headsize,tailsize,genenumber,group,functionset,lengthofORF,computationresult,genesize,functionparameter,data1,terminalset);%������ֵ
fitness=[];
fitness=figurefitness(fitness,observationresult,computationresult,populationsize);%������Ӧ��
newgroup={};
[newgroup,bestindividual,bestfitness,bestfitnessarray,averagefitnessarray,historybestindividual,historybestfitness]=selection(newgroup,group,fitness,populationsize,bestfitnessarray,averagefitnessarray,historybestindividual,historybestfitness,generation,participatenum);%ѡ�����
newgroup=geneticoperation(headsize,newgroup,mutationrate,populationsize,genesize,genenumber,ISelementlength,IStranspositionrate,RIStranspositionrate,RISelementlength,genetranspositionrate,generecombinationrate,oneprecombinationrate,twoprecombinationrate,terminalset,functionset,tailsize);%�Ŵ�����
group=newgroup;
end
char(historybestindividual)
computationresult1=[];
fitness1=[];
computationresult1=computevalue(1,headsize,tailsize,genenumber,historybestindividual,functionset,lengthofORF,computationresult1,genesize,functionparameter,data2,terminalset);%������ֵ
fitness1=figurefitness(fitness1,observationresult1,computationresult1,1);%������Ӧ��
a=1:size(observationresult1,1);
plot(a,computationresult1,'-r+',a,observationresult1,'-bx');
xlabel('NO.');
ylabel('value');
legend('������ֵ','�۲�ֵ');
result1(qqq)=historybestfitness;
result2(qqq)=fitness1;
result3(qqq,:)=historybestindividual;
result4(1,:,qqq)=computationresult1;
result4(2,:,qqq)=observationresult1;
end