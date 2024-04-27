%计算数值（包括计算ORF长度和计算数值）
function computationresult=computevalue(populationsize,headsize,tailsize,genenumber,group,functionset,lengthofORF,computationresult,genesize,functionparameter,data1,terminalset)
for index5=1:populationsize%解码染色体，计算每个基因的有效长度
    for index6=0:genenumber-1 
            token=0;%此变量代表标志（q）
            position=0;%指示当前搜索位置（i）
        for index7=1:headsize+tailsize 
               if  ismember(group(index5,index6*genesize+index7),functionset)==1
                   token=token+functionparameter(ismember(functionset,group(index5,index6*genesize+index7)));
                   position=position+1;
               else
                   position=position+1;
               end
               if position>token
                   break;
               end
        end
       lengthofORF(index5,index6+1)=position;
    end
end

temp=[];%临时存放有效的基因
tempresult=[];%存放数值计算的中间结果
computationresult=[];%存放由公式计算的结果
[m,n]=size(data1);
tempterminal=[];
tempfunction=[];

for dataindex=1:m    
    for index8=1:populationsize
        for index9=1:genenumber
           temp=group(index8,((index9-1)*genesize+1):((index9-1)*genesize+lengthofORF(index8,index9)));
           tempcopy=[];
           templength=lengthofORF(index8,index9);
           functionnum=0;%记录有效基因中函数符的个数，函数符的数量就是从后向前计算的次数
           for index78=1:lengthofORF(index8,index9)
               if ismember(temp(index78),functionset)==1
                   functionnum=functionnum+1;
               end
           end
     for index10=1:numel(temp)  
        if ismember(temp(index10),terminalset)==1
             tempcopy(index10)=data1(dataindex,ismember(terminalset,temp(index10)));
        end
     end
                if numel(tempcopy)==1
                    tempresult(1)=tempcopy(1);
                else   
                    for inde56=1:functionnum   
                     for index11=numel(temp):-1:1
                       if ismember(temp(index11),functionset)==1
                        tempfunction=temp{index11};
                        switch tempfunction
                            case '+'
                                tempcopy(index11)=tempcopy(templength)+tempcopy(templength-1);
                                templength=templength-functionparameter(ismember(functionset,tempfunction));
                                temp{index11}='@';%将计算过的函数符号换成其他符号
                            case '-'
                                tempcopy(index11)=tempcopy(templength)-tempcopy(templength-1);
                                templength=templength-functionparameter(ismember(functionset,tempfunction));
                                temp{index11}='@';%将计算过的函数符号换成其他符号
                            case '*'
                                tempcopy(index11)=tempcopy(templength)*tempcopy(templength-1);
                                 templength=templength-functionparameter(ismember(functionset,tempfunction));
                                 temp{index11}='@';%将计算过的函数符号换成其他符号
                            case '/'
                                tempcopy(index11)=tempcopy(templength)/tempcopy(templength-1);
                                templength=templength-functionparameter(ismember(functionset,tempfunction));
                                temp{index11}='@';%将计算过的函数符号换成其他符号
                            case 'S'
                                tempcopy(index11)=sinh(tempcopy(templength));
                                templength=templength-functionparameter(ismember(functionset,tempfunction));
                                temp{index11}='@';%将计算过的函数符号换成其他符号
                           case 'C'
                                tempcopy(index11)=cosh(tempcopy(templength));
                                templength=templength-functionparameter(ismember(functionset,tempfunction));
                                temp{index11}='@';%将计算过的函数符号换成其他符号
                           case 'T'
                                tempcopy(index11)=tan(tempcopy(templength));
                                templength=templength-functionparameter(ismember(functionset,tempfunction));
                                 temp{index11}='@';%将计算过的函数符号换成其他符号
                           case 'L'
                                tempcopy(index11)=log10(tempcopy(templength));
                                templength=templength-functionparameter(ismember(functionset,tempfunction));
                                temp{index11}='@';%将计算过的函数符号换成其他符号
                           case 'G'
                                tempcopy(index11)=log2(tempcopy(templength));
                                templength=templength-functionparameter(ismember(functionset,tempfunction));
                                temp{index11}='@';%将计算过的函数符号换成其他符号
                          case 'X'
                                tempcopy(index11)=log(tempcopy(templength));
                                templength=templength-functionparameter(ismember(functionset,tempfunction));
                                temp{index11}='@';%将计算过的函数符号换成其他符号
                          case 'Q'
                                tempcopy(index11)=sqrt(tempcopy(templength));
                                templength=templength-functionparameter(ismember(functionset,tempfunction));
                                temp{index11}='@';%将计算过的函数符号换成其他符号
                          case 'A'
                                tempcopy(index11)=abs(tempcopy(templength));
                                templength=templength-functionparameter(ismember(functionset,tempfunction));
                                temp{index11}='@';%将计算过的函数符号换成其他符号
                          case 'E'
                                tempcopy(index11)=exp(tempcopy(templength));
                                templength=templength-functionparameter(ismember(functionset,tempfunction));
                                temp{index11}='@';%将计算过的函数符号换成其他符号
                          case 'F'
                                tempcopy(index11)=power(tempcopy(templength),2);
                                templength=templength-functionparameter(ismember(functionset,tempfunction));
                                temp{index11}='@';%将计算过的函数符号换成其他符号
                           case 'M'
                                tempcopy(index11)=max(tempcopy(templength),tempcopy(templength-1));
                                templength=templength-functionparameter(ismember(functionset,tempfunction));
                                temp{index11}='@';%将计算过的函数符号换成其他符号
                           case 'I'
                                tempcopy(index11)=min(tempcopy(templength),tempcopy(templength-1));
                                templength=templength-functionparameter(ismember(functionset,tempfunction));
                                temp{index11}='@';%将计算过的函数符号换成其他符号
                           case 'O'
                                tempcopy(index11)=(tempcopy(templength)+tempcopy(templength-1))/2;
                                templength=templength-functionparameter(ismember(functionset,tempfunction));
                                temp{index11}='@';%将计算过的函数符号换成其他符号
                             case 'R'
                                tempcopy(index11)=1/tempcopy(templength);
                                templength=templength-functionparameter(ismember(functionset,tempfunction));
                                temp{index11}='@';%将计算过的函数符号换成其他符号      
                           case 'D'
                                tempcopy(index11)=tempcopy(templength)*2;
                                templength=templength-functionparameter(ismember(functionset,tempfunction));
                                temp{index11}='@';%将计算过的函数符号换成其他符号 
                           case 'P'
                                tempcopy(index11)=tempcopy(templength)*3;
                                templength=templength-functionparameter(ismember(functionset,tempfunction));
                                temp{index11}='@';%将计算过的函数符号换成其他符号 
                             case '~'
                                tempcopy(index11)=-tempcopy(templength);
                                templength=templength-functionparameter(ismember(functionset,tempfunction));
                                temp{index11}='@';%将计算过的函数符号换成其他符号  
                             case 'N'
                                tempcopy(index11)=sin(templength);
                                templength=templength-functionparameter(ismember(functionset,tempfunction));
                                temp{index11}='@';%将计算过的函数符号换成其他符号  
                             case 'J'
                                tempcopy(index11)=cos(templength);
                                templength=templength-functionparameter(ismember(functionset,tempfunction));
                                temp{index11}='@';%将计算过的函数符号换成其他符号  
                        end%与switch配对
                        break;
                       end%与if配对
                     end%与内层for配对
                     end%与外层for配对 
                      if templength==1
                          tempresult(index9)=tempcopy(1);
                      end
                end%与if/else配对            
        end%与index9配对
     computationresult(dataindex,index8)=sum(tempresult);  
     %computationresult(dataindex,index8)=prod(tempresult);
    end%与index8配对
end%与dataindex配对