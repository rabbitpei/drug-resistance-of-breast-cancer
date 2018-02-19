clear;
clc;
close all;

fpi=fopen('processed_tongji_07282016.txt');
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

psize=size(tempmimi);
qthresh=0.01;
clear pvalue;
for i=1:psize(1)
    presample=[reshape(tempmimi(i,2,:),1,psize(3)),reshape(tempmimi(i,3,:),1,psize(3))];
    postsample=[reshape(tempmimi(i,6,:),1,psize(3)),reshape(tempmimi(i,7,:),1,psize(3))];  
    [h,pvalue(i),ci]=ttest2(presample,postsample,qthresh);    
    if mod(i,1000)==0
        i
    end
end
[pv,idx]=sort(pvalue);
fw=fopen('differential_genes23&67.txt','w');
for i=1:200
    fprintf(fw,'%s\n',pipi{idx(i)});
end
fclose(fw);
fw=fopen('differential_genes_23&67_005.txt','w');
for i=1:psize(1)
    if pv(i)>0.05
        break
    end
    fprintf(fw,'%s\n',pipi{idx(i)});
end
fclose(fw);

psize=size(tempmimi);
qthresh=0.01;
clear pvalue;
for i=1:psize(1)
    presample=[reshape(tempmimi(i,9,:),1,psize(3)),reshape(tempmimi(i,10,:),1,psize(3))];
    postsample=[reshape(tempmimi(i,6,:),1,psize(3)),reshape(tempmimi(i,7,:),1,psize(3))];  
    [h,pvalue(i),ci]=ttest2(presample,postsample,qthresh);    
    if mod(i,1000)==0
        i
    end
end
[pv,idx]=sort(pvalue);
fw=fopen('differential_genes67&910.txt','w');
for i=1:200
    fprintf(fw,'%s\n',pipi{idx(i)});
end
fclose(fw);
fw=fopen('differential_genes_67&910_005.txt','w');
for i=1:psize(1)
    if pv(i)>0.05
        break
    end
    fprintf(fw,'%s\n',pipi{idx(i)});
end
fclose(fw);

% 
% psize=size(temppig);
% qthresh=0.01;
% for i=1:psize(1)
%     for j=1:psize(2)
%         [h(j),pvalue(i,j),ci(j,:)]=ttest2(temppig(i,j,:),tempmimi(i,j,:),qthresh);
%     end
%     if mod(i,100)==0
%         i
%     end
% end
% for j=1:psize(2)
%     [pv(:,j),idx(:,j)]=sort(pvalue(:,j));
% end
% 
% j=zeros(psize(2),1);
% k=ones(psize(2),1);
% clear index;
% for t=1:psize(2)
%     %filename=['time',int2str(t),'.txt'];
%     %fw=fopen(filename,'w');
%     while pv(k(t),t)<(k(t)/psize(1))*qthresh
%         if std(tempmimi(idx(k(t),t),t,:))/std(temppig(idx(k(t),t),t,:))<4
%             k(t)=k(t)+1;
%             continue
%         end
%         j(t)=j(t)+1;
%         usefulPvalue(j(t),t,:)=pv(k(t),:);
%         index(j(t),t)=idx(k(t),t);
%         usefulRNA_id(j(t),t)=pipi(idx(k(t),t));
%         usefulRNA_pig(j(t),t,:)=temppig(idx(k(t),t),t,:);
%         usefulRNA_mimi(j(t),t,:)=tempmimi(idx(k(t),t),t,:);
%         usefulRNA_healthy(j(t),t,:)=pprofile(idx(k(t),t),:);
%         usefulRNA_sick(j(t),t,:)=mprofile(idx(k(t),t),:);
%         k(t)=k(t)+1;
%     end
%     %fw.close();
% end
% jk=k'
% jy=j'
% a=b
