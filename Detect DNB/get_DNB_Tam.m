clear;
clc;
close all;

fpi=fopen('processed_tongji_07262017.txt');
hline1 = textscan(fpi, '%s', 1, 'delimiter', '\n');
field1=textscan(hline1{1}{1},'%s');
format='%s';
for i=2:181
    format=[format,' %f'];
end
lines =textscan(fpi, format,100000,'delimiter', '\t');
pipi=lines{1};
profile = [];
for i = 2 : 181
    profile = [profile, lines{i}];
end
fclose(fpi);
psize=size(profile);             %22690*56
zprofile=zscore(profile');
zprofile=zprofile';

for i=1:12
    temppig(:,i,:)=zprofile(:,3*5*(i-1)+1:3*5*(i-1)+5);
    tempmimi(:,i,:)=zprofile(:,3*5*(i-1)+11:3*5*(i-1)+15);
    control(:,5*(i-1)+1:5*(i-1)+5)=zprofile(:,3*5*(i-1)+1:3*5*(i-1)+5);
    Tamoxifen(:,5*(i-1)+1:5*(i-1)+5)=zprofile(:,3*5*(i-1)+11:3*5*(i-1)+15);
end

pprofile=control;
mprofile=Tamoxifen;

psize=size(temppig);
qthresh=0.01;
for i=1:psize(1)
    for j=1:psize(2)
        [h(j),pvalue(i,j),ci(j,:)]=ttest2(temppig(i,j,:),tempmimi(i,j,:),qthresh);
        mysd(i,j)=std(tempmimi(i,j,:));
    end
    if mod(i,100)==0
        i
    end
end
for j=1:psize(2)
    [pv(:,j),idx(:,j)]=sort(pvalue(:,j));
end


j=zeros(psize(2),1);
k=ones(psize(2),1);
clear index;
for t=1:psize(2)
    %filename=['time',int2str(t),'.txt'];
    %fw=fopen(filename,'w');
    while pv(k(t),t)<(k(t)/psize(1))*qthresh
        if std(tempmimi(idx(k(t),t),t,:))/std(temppig(idx(k(t),t),t,:))<4
            k(t)=k(t)+1;
            continue
        end
        j(t)=j(t)+1;
        usefulPvalue(j(t),t,:)=pv(k(t),:);
        index(j(t),t)=idx(k(t),t);
        usefulRNA_id(j(t),t)=pipi(idx(k(t),t));
        usefulRNA_pig(j(t),t,:)=temppig(idx(k(t),t),t,:);
        usefulRNA_mimi(j(t),t,:)=tempmimi(idx(k(t),t),t,:);
        usefulRNA_healthy(j(t),t,:)=pprofile(idx(k(t),t),:);
        usefulRNA_sick(j(t),t,:)=mprofile(idx(k(t),t),:);
        k(t)=k(t)+1;
    end
    %fw.close();
end
jk=k'
jy=j'

clear Y4 Z4 T4
if jy(4)>1
    Y4=pdist(usefulRNA_mimi(1: jy(4),4,:),@cofun);
    Z4=linkage(Y4);
    T4=cluster(Z4,'maxclust',180);
    num4=max(T4);
    groupnum4=zeros(1,num4);
    for i=1:length(T4)
        groupnum4(T4(i))=groupnum4(T4(i))+1;
    end
    max(groupnum4)
end

mkdir('cluster_gene_group_Tam');
path=cd;
newpath=strcat(path,'/','cluster_gene_group_Tam','/');

clear A;
fla=0;
validgroupnum=0;

if jy(4)>1
    filename=strcat(newpath,'time',int2str(4),'.txt');
    fpi=fopen(filename,'wt');
    candidategroup4=0;
    for i=1:num4
        fprintf(fpi,'group num: %d:\n',i);
        clear tm;
        clear k; clear tmm; clear tp; clear zjj;
        stdRNA=zeros(psize(2),1);
        mixpcc=zeros(psize(2),1);
        k=0;
        
        for j=1:length(T4)
            if i==T4(j)
                fprintf(fpi,'%d\t%s\n',j,usefulRNA_id{j,4});
                k=k+1;
                %             tm(k,:)=usefulRNA2_z(j,:);
                tm(k,:)=usefulRNA_sick(j,4,:);
                tp(k,:)=usefulRNA_healthy(j,4,:);
                
                for x=1:psize(2)
                    stdRNA(x)=stdRNA(x)+std(tm(k,5*(x-1)+1:5*x));  %/(pvl(k,x)+0.000001)^0.5;
                end
            end
        end
        if k<2
            fprintf(fpi,'\n');
            continue
        end
        
        if k<80
            continue;
        end
        tmm=tm;
        
        for x=1:psize(2)
            pccall=corr(tm(:,5*(x-1)+1:5*x)');
            pccrna(x)=(sum(sum(abs(pccall)))-k)/2;
            avestd(x)=stdRNA(x)/k;
        end
        
        %line=line+1;
        fprintf(fpi,'\n');
        validgroupnum=validgroupnum+1;
        stdRNA=stdRNA/k;
        pccrna=pccrna/(k*(k-1))*2;
        
        othergenes_mimi=setdiff(mprofile,tm,'rows');
        mimisize=size(othergenes_mimi);
        for rand_time=1:10
            rand_idx=randsample([1:mimisize(1)],200);
            rand_mimi=othergenes_mimi(rand_idx,:);
            for ii=1:k
                for jj=1:200
                    zscoretm=rand_mimi(jj,:);
                    for x=1:psize(2)
                        mixpcc(x)=mixpcc(x)+abs(corr(tm(ii,5*(x-1)+1:5*x)',...
                            zscoretm(5*(x-1)+1:5*x)'));
                    end
                end
            end
            mixpcc=mixpcc/k/200/10;
        end
        
        tag=0;
        for time=2:psize(2)-1
            if (stdRNA(time)>stdRNA(time-1))&&(stdRNA(time)>stdRNA(time+1))
                %  &&(pccrna(time)>pccrna(time-1))&&(pccrna(time)>pccrna(time+1))
                % &&(mixpcc(time)<mixpcc(time-1))&&(mixpcc(time)<mixpcc(time+1))
                time
                tag=1;
                break;
            end
        end
        if tag==0
            continue;
        end
        
        candidategroup4=candidategroup4+1;
        
        compositindex=zeros(psize(2),1);
        for s=1:psize(2)
            compositindex(s)=stdRNA(s)*pccrna(s);  %/mixpcc(s);
            
        end
        
        figure(validgroupnum);
        t=[1:psize(2)];
        
        subplot(1,4,1);
        plot(t,stdRNA,'r');
        tmp_name=[int2str(i),'-std4'];
        title(tmp_name);
        
        subplot(1,4,2);
        plot(t,pccrna,'r');
        tmp_name=[int2str(i),'-inpcc4'];
        title(tmp_name);
        
        subplot(1,4,3);
        plot(t,mixpcc,'r');
        tmp_name=[int2str(i),'-outpcc4'];
        title(tmp_name);
        
        subplot(1,4,4);
        plot(t,compositindex,'r');
        tmp_name=[int2str(i),'-composite index4'];
        title(tmp_name);
    end
    fclose(fpi);
end
