%�Ŵ���������������죬�����������
function newgroup=geneticoperation(headsize,newgroup,mutationrate,populationsize,genesize,genenumber,ISelementlength,IStranspositionrate,RIStranspositionrate,RISelementlength,genetranspositionrate,generecombinationrate,oneprecombinationrate,twoprecombinationrate,terminalset,functionset,tailsize)
for index13=1:populationsize*mutationrate%�������
    mutationpoint=[];
    tempindividual=[];
    selectedchrom=[];
   selectedchrom=randint(1,1,[1,populationsize]);
   tempindividual=newgroup(selectedchrom,:);%���ѡ��һ��Ⱦɫ�壬��¼Ⱦɫ��ı�ţ����ҽ��丳ֵ��tempindividual
    mutationpoint=randint(1,1,[1,genesize*genenumber]);%�����������λ��
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

for index14=1:populationsize*IStranspositionrate%ISת��
     tempindividual=[];
     selectedchromosome=[];
     selectedpoint=[];
     temphead={};%��ű��޸Ļ������ʱͷ��
      tempISelement=[];
      selectedgene=[];
      ISelementstartpoint=[];
      selectedchromosome=randint(1,1,[1,populationsize]);
      tempindividual=newgroup(selectedchromosome,:);%���ѡ��һ��Ⱦɫ�壬��¼Ⱦɫ��ı�ţ����ҽ��丳ֵ��tempindividual
      chosenISelementlength=ISelementlength(randint(1,1,[1,numel(ISelementlength)]));%���ѡ��ISԪ�صĳ���
      ISelementstartpoint=randint(1,1,[1,(genesize*genenumber-chosenISelementlength+1)]);%���ѡ��ISԪ�ص���ʼλ��
      tempISelement=tempindividual(ISelementstartpoint:ISelementstartpoint+(chosenISelementlength-1));%�洢ISԪ�صĸ���
      selectedgene=randint(1,1,[1,genenumber]);%���ѡ��ISԪ�ز���Ļ������
      selectedpoint=randint(1,1,[2,headsize-chosenISelementlength+1]);%ѡ�񱻲���Ԫ����ͷ����λ��,����ӷ��2��ʼ����(ע��ͷ���Ĳ���λ��)
      temphead(1:(selectedpoint-1))=tempindividual(((selectedgene-1)*genesize+1):((selectedgene-1)*genesize+(selectedpoint-1)));
      temphead(selectedpoint:(selectedpoint+chosenISelementlength-1))=tempISelement;
      temphead((selectedpoint+chosenISelementlength):headsize)=tempindividual(((selectedgene-1)*genesize+selectedpoint):((selectedgene-1)*genesize+headsize-chosenISelementlength));
      tempindividual(((selectedgene-1)*genesize+1):((selectedgene-1)*genesize+headsize))=temphead;
      newgroup(selectedchromosome,:)=tempindividual;
end
for index15=1:populationsize*RIStranspositionrate%RISת��
     tempindividual=[];
     selectedchromosome=[];
     temphead={};%��ű��޸Ļ������ʱͷ��
     tempRISelement=[];
      selectedgene=[];
      tempgene=[];
      randscanpoint=[];
      RISelementstartpoint=[];
      selectedchromosome=randint(1,1,[1,populationsize]);%���ѡ����Ҫ�޸ĵ�Ⱦɫ��
      tempindividual=newgroup(selectedchromosome,:);%��ֵ��tempindividual
      selectedgene=randint(1,1,[1,genenumber]);%���ѡ��RISԪ�ز���Ļ������
      tempgene=tempindividual(((selectedgene-1)*genesize+1):((selectedgene-1)*genesize+genesize));
      randscanpoint=randint(1,1,[1,(genenumber*genesize-tailsize)]);%���ѡ��RISԪ�ص�ɨ����ʼ��
      for index45=randscanpoint:genenumber*genesize%��ɨ����ʼλ�ÿ�ʼѰ�ҵ�һ������
          if ismember(tempindividual(index45),functionset)==1
              RISelementstartpoint=index45;
              break;
          end
      end
      if isempty(RISelementstartpoint)==0 
             chosenRISelementlength=RISelementlength(randint(1,1,[1,numel(RISelementlength)]));%���ѡ��RISԪ�صĳ���
            tempRISelement=tempindividual(RISelementstartpoint:RISelementstartpoint+chosenRISelementlength-1);
            temphead(1:chosenRISelementlength)=tempRISelement;
            temphead((chosenRISelementlength+1):headsize)=tempgene(1:(headsize-chosenRISelementlength));
            tempindividual(((selectedgene-1)*genesize+1):((selectedgene-1)*genesize+headsize))=temphead;
            newgroup(selectedchromosome,:)=tempindividual;
      end   
end
for index010=1:populationsize*genetranspositionrate%geneת��
     tempindividual=[];
     tempindividual1={};
     selectedchromosome=[];
     %temphead={};%��ű��޸Ļ������ʱͷ��
     tempRISelement=[];
      selectedgene=[];
      tempgene=[];
      randscanpoint=[];
      RISelementstartpoint=[];
      selectedchromosome=randint(1,1,[1,populationsize]);%���ѡ����Ҫ�޸ĵ�Ⱦɫ��
      tempindividual=newgroup(selectedchromosome,:);%��ֵ��tempindividual
      selectedgene=randint(1,1,[1,genenumber]);%���ѡ��������
      tempgene=tempindividual(((selectedgene-1)*genesize+1):((selectedgene-1)*genesize+genesize));
      tempindividual1(1:genesize)=tempgene;
      tempindividual1((genesize+1):(genesize*selectedgene))=tempindividual(1:(genesize*(selectedgene-1)));
      tempindividual1((genesize*selectedgene+1):(genesize*genenumber))=tempindividual((genesize*selectedgene+1):(genenumber*genesize));
      newgroup(selectedchromosome,:)=tempindividual1;
end
for index17=1:populationsize*twoprecombinationrate%���㽻��
      cross=[];%��¼����Ⱦɫ��ı��
      crosspoint=[];
      tempindividual1=[];
      tempindividual2=[];
      temp=[];
      numtemp=[];
      tempcross=randperm(populationsize);
     cross=tempcross(1:2);%���ѡ�������μ������Ⱦɫ��
     tempcrosspoint=randperm(genenumber*genesize);
     crosspoint=tempcrosspoint(1:2);%���ѡ�����������
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
for index16=1:populationsize*oneprecombinationrate%���㽻��
      cross=[];
      crosspoint=[];
      tempindividual1=[];
      tempindividual2=[];
      temp=[];
     tempcross=randperm(populationsize);
     cross=tempcross(1:2);%���ѡ�������μ������Ⱦɫ��
     crosspoint=randint(1,1,[1,genenumber*genesize]);
     tempindividual1=newgroup(cross(1),:);
     tempindividual2=newgroup(cross(2),:);
     temp=tempindividual1(crosspoint+1:genenumber*genesize);
     tempindividual1(crosspoint+1:genenumber*genesize)=tempindividual2(crosspoint+1:genenumber*genesize);
     tempindividual2(crosspoint+1:genenumber*genesize)=temp;
    newgroup(cross(1),:)=tempindividual1;
     newgroup(cross(2),:)=tempindividual2;
end
for index33=1:populationsize*generecombinationrate%��������
      cross=[];
      tempindividual1=[];
      tempindividual2=[];
      temp=[];
      selectedgene=[];
     tempcross=randperm(populationsize);
     cross=tempcross(1:2);%���ѡ�������μ������Ⱦɫ��
     tempindividual1=newgroup(cross(1),:);
     tempindividual2=newgroup(cross(2),:);
     selectedgene=randint(1,1,[1,genenumber]);%���ѡ�񽻻�����ĺ���
     temp=tempindividual1((selectedgene-1)*genesize+1:(selectedgene-1)*genesize+genesize);
     tempindividual1((selectedgene-1)*genesize+1:(selectedgene-1)*genesize+genesize)=tempindividual2((selectedgene-1)*genesize+1:(selectedgene-1)*genesize+genesize);
     tempindividual2((selectedgene-1)*genesize+1:(selectedgene-1)*genesize+genesize)=temp;
      newgroup(cross(1),:)=tempindividual1;
     newgroup(cross(2),:)=tempindividual2;
end