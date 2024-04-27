%遗传操作，包括点变异，串变异和重组
function newgroup=geneticoperation(headsize,newgroup,mutationrate,populationsize,genesize,genenumber,ISelementlength,IStranspositionrate,RIStranspositionrate,RISelementlength,genetranspositionrate,generecombinationrate,oneprecombinationrate,twoprecombinationrate,terminalset,functionset,tailsize)
for index13=1:populationsize*mutationrate%变异操作
    mutationpoint=[];
    tempindividual=[];
    selectedchrom=[];
   selectedchrom=randint(1,1,[1,populationsize]);
   tempindividual=newgroup(selectedchrom,:);%随机选择一个染色体，记录染色体的编号，并且将其赋值给tempindividual
    mutationpoint=randint(1,1,[1,genesize*genenumber]);%随机产生变异位置
    if mod(mutationpoint,genesize)==0
        tempindividual(mutationpoint)=terminalset(randint(1,1,[1,numel(terminalset)]));
    else
        if mod(mutationpoint,genesize)>headsize
            tempindividual(mutationpoint)=terminalset(randint(1,1,[1,numel(terminalset)]));
        else
            if rand()>0.5
                tempindividual(mutationpoint)=functionset(randint(1,1,[1,numel(functionset)]));
            else 
                tempindividual(mutationpoint)=terminalset(randint(1,1,[1,numel(terminalset)]));
            end
        end
    end
    newgroup(selectedchrom,:)=tempindividual;
end

for index14=1:populationsize*IStranspositionrate%IS转座
     tempindividual=[];
     selectedchromosome=[];
     selectedpoint=[];
     temphead={};%存放被修改基因的临时头部
      tempISelement=[];
      selectedgene=[];
      ISelementstartpoint=[];
      selectedchromosome=randint(1,1,[1,populationsize]);
      tempindividual=newgroup(selectedchromosome,:);%随机选择一个染色体，记录染色体的编号，并且将其赋值给tempindividual
      chosenISelementlength=ISelementlength(randint(1,1,[1,numel(ISelementlength)]));%随机选择IS元素的长度
      ISelementstartpoint=randint(1,1,[1,(genesize*genenumber-chosenISelementlength+1)]);%随机选择IS元素的起始位置
      tempISelement=tempindividual(ISelementstartpoint:ISelementstartpoint+(chosenISelementlength-1));%存储IS元素的副本
      selectedgene=randint(1,1,[1,genenumber]);%随机选择被IS元素插入的基因号码
      selectedpoint=randint(1,1,[2,headsize-chosenISelementlength+1]);%选择被插入元素在头部的位置,插入接缝从2开始计数(注意头部的插入位置)
      temphead(1:(selectedpoint-1))=tempindividual(((selectedgene-1)*genesize+1):((selectedgene-1)*genesize+(selectedpoint-1)));
      temphead(selectedpoint:(selectedpoint+chosenISelementlength-1))=tempISelement;
      temphead((selectedpoint+chosenISelementlength):headsize)=tempindividual(((selectedgene-1)*genesize+selectedpoint):((selectedgene-1)*genesize+headsize-chosenISelementlength));
      tempindividual(((selectedgene-1)*genesize+1):((selectedgene-1)*genesize+headsize))=temphead;
      newgroup(selectedchromosome,:)=tempindividual;
end
for index15=1:populationsize*RIStranspositionrate%RIS转座
     tempindividual=[];
     selectedchromosome=[];
     temphead={};%存放被修改基因的临时头部
     tempRISelement=[];
      selectedgene=[];
      tempgene=[];
      randscanpoint=[];
      RISelementstartpoint=[];
      selectedchromosome=randint(1,1,[1,populationsize]);%随机选择需要修改的染色体
      tempindividual=newgroup(selectedchromosome,:);%赋值给tempindividual
      selectedgene=randint(1,1,[1,genenumber]);%随机选择被RIS元素插入的基因号码
      tempgene=tempindividual(((selectedgene-1)*genesize+1):((selectedgene-1)*genesize+genesize));
      randscanpoint=randint(1,1,[1,(genenumber*genesize-tailsize)]);%随机选择RIS元素的扫描起始点
      for index45=randscanpoint:genenumber*genesize%从扫描起始位置开始寻找第一个函数
          if ismember(tempindividual(index45),functionset)==1
              RISelementstartpoint=index45;
              break;
          end
      end
      if isempty(RISelementstartpoint)==0 
             chosenRISelementlength=RISelementlength(randint(1,1,[1,numel(RISelementlength)]));%随机选择RIS元素的长度
            tempRISelement=tempindividual(RISelementstartpoint:RISelementstartpoint+chosenRISelementlength-1);
            temphead(1:chosenRISelementlength)=tempRISelement;
            temphead((chosenRISelementlength+1):headsize)=tempgene(1:(headsize-chosenRISelementlength));
            tempindividual(((selectedgene-1)*genesize+1):((selectedgene-1)*genesize+headsize))=temphead;
            newgroup(selectedchromosome,:)=tempindividual;
      end   
end
for index010=1:populationsize*genetranspositionrate%gene转座
     tempindividual=[];
     tempindividual1={};
     selectedchromosome=[];
     %temphead={};%存放被修改基因的临时头部
     tempRISelement=[];
      selectedgene=[];
      tempgene=[];
      randscanpoint=[];
      RISelementstartpoint=[];
      selectedchromosome=randint(1,1,[1,populationsize]);%随机选择需要修改的染色体
      tempindividual=newgroup(selectedchromosome,:);%赋值给tempindividual
      selectedgene=randint(1,1,[1,genenumber]);%随机选择基因号码
      tempgene=tempindividual(((selectedgene-1)*genesize+1):((selectedgene-1)*genesize+genesize));
      tempindividual1(1:genesize)=tempgene;
      tempindividual1((genesize+1):(genesize*selectedgene))=tempindividual(1:(genesize*(selectedgene-1)));
      tempindividual1((genesize*selectedgene+1):(genesize*genenumber))=tempindividual((genesize*selectedgene+1):(genenumber*genesize));
      newgroup(selectedchromosome,:)=tempindividual1;
end
for index17=1:populationsize*twoprecombinationrate%两点交叉
      cross=[];%记录两个染色体的编号
      crosspoint=[];
      tempindividual1=[];
      tempindividual2=[];
      temp=[];
      numtemp=[];
      tempcross=randperm(populationsize);
     cross=tempcross(1:2);%随机选择两个参加重组的染色体
     tempcrosspoint=randperm(genenumber*genesize);
     crosspoint=tempcrosspoint(1:2);%随机选择两个重组点
     if crosspoint(1)>crosspoint(2)
         numtemp=crosspoint(1);
         crosspoint(1)=crosspoint(2);
         crosspoint(2)=numtemp;
     end
     tempindividual1=newgroup(cross(1),:);
     tempindividual2=newgroup(cross(2),:);
     temp=tempindividual1(crosspoint(1):crosspoint(2));
     tempindividual1(crosspoint(1):crosspoint(2))=tempindividual2(crosspoint(1):crosspoint(2));
     tempindividual2(crosspoint(1):crosspoint(2))=temp;
     newgroup(cross(1),:)=tempindividual1;
     newgroup(cross(2),:)=tempindividual2;
end
for index16=1:populationsize*oneprecombinationrate%单点交叉
      cross=[];
      crosspoint=[];
      tempindividual1=[];
      tempindividual2=[];
      temp=[];
     tempcross=randperm(populationsize);
     cross=tempcross(1:2);%随机选择两个参加重组的染色体
     crosspoint=randint(1,1,[1,genenumber*genesize]);
     tempindividual1=newgroup(cross(1),:);
     tempindividual2=newgroup(cross(2),:);
     temp=tempindividual1(crosspoint+1:genenumber*genesize);
     tempindividual1(crosspoint+1:genenumber*genesize)=tempindividual2(crosspoint+1:genenumber*genesize);
     tempindividual2(crosspoint+1:genenumber*genesize)=temp;
    newgroup(cross(1),:)=tempindividual1;
     newgroup(cross(2),:)=tempindividual2;
end
for index33=1:populationsize*generecombinationrate%基因重组
      cross=[];
      tempindividual1=[];
      tempindividual2=[];
      temp=[];
      selectedgene=[];
     tempcross=randperm(populationsize);
     cross=tempcross(1:2);%随机选择两个参加重组的染色体
     tempindividual1=newgroup(cross(1),:);
     tempindividual2=newgroup(cross(2),:);
     selectedgene=randint(1,1,[1,genenumber]);%随机选择交换基因的号码
     temp=tempindividual1((selectedgene-1)*genesize+1:(selectedgene-1)*genesize+genesize);
     tempindividual1((selectedgene-1)*genesize+1:(selectedgene-1)*genesize+genesize)=tempindividual2((selectedgene-1)*genesize+1:(selectedgene-1)*genesize+genesize);
     tempindividual2((selectedgene-1)*genesize+1:(selectedgene-1)*genesize+genesize)=temp;
      newgroup(cross(1),:)=tempindividual1;
     newgroup(cross(2),:)=tempindividual2;
end