clear;
clc;
close all;

fpi=fopen('MCF7_RPKM_180.txt');
hline1 = textscan(fpi, '%s', 1, 'delimiter', '\n');
field1=textscan(hline1{1}{1},'%s');
format='%s';
for i=2:181
    format=[format,' %s'];
end
lines =textscan(fpi, format,100000,'delimiter', '\t');
pipi=lines{1};
profile = [];
for i = 2 : 181
    profile = [profile, lines{i}];
end
fclose(fpi);
psize=size(profile);             %22690*56

processed_id=pipi;
processed_profile=zeros(size(profile));

k=0;
for i=1:psize(1)
    flag=1;
    tmp=zeros(psize(2),1);
    for j=1:12
        tp=str2float(profile(i,3*5*(j-1)+1:3*5*(j-1)+5));
        tm=str2float(profile(i,3*5*(j-1)+11:3*5*(j-1)+15));
        if length(find(tp<0.001))>2||length(find(tm<0.001))>2
            flag=0;
            break
        else 
            if length(find(tp<0.001))>0
                tp(find(tp<0.001))=mean(tp(find(tp>0.001)));
            end
            tmp(3*5*(j-1)+1:3*5*(j-1)+5)=tp;
            if length(find(tm<0.001))>0
                tm(find(tm<0.001))=mean(tm(find(tm>0.001)));
            end
            tmp(3*5*(j-1)+11:3*5*(j-1)+15)=tm;          
        end
    end
    if flag==1
        k=k+1;
        processed_id{k}=pipi{i};
        processed_profile(k,:)=tmp;
    end  
    if mod(i,1000)==0
        i
    end
end

fw=fopen('processed_tongji_07262017.txt','w');
processed_profile=processed_profile(1:k,:);
processed_id=processed_id(1:k);
for i=1:k
    fprintf(fw,'%s',processed_id{i});
    for j=1:size(processed_profile,2)
        fprintf(fw,'\t%f',processed_profile(i,j));
    end
    fprintf(fw,'\n');
end
fclose(fw);
