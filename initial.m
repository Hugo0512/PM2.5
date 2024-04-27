%种群初始化
function group=initial(group,populationsize,genenumber,headsize,genesize,functionset,terminalset,tailsize)
for index1=1:populationsize%产生初始种群
    for index2=0:genenumber-1
        for index3=1:headsize
            if rand()>0.5
            group(index1,index2*genesize+index3)=functionset(randint(1,1,[1,numel(functionset)]));
            else
            group(index1,index2*genesize+index3)=terminalset(randint(1,1,[1,numel(terminalset)]));
            end
        end
        for index4=headsize+1:headsize+tailsize
            group(index1,index2*genesize+index4)=terminalset(randint(1,1,[1,numel(terminalset)]));
        end
    end
end