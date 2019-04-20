Data=importdata('Results.txt');
rslKs=4;rslJs=10;rslJ0=4;
Data=Data(:,1:4);
Data(:,4)=Data(:,4)./sum(Data(:,4));
cdf=0;
len=length(Data(:,4));
Mcdf=zeros(1,len);
for pf=1:len
cdf=cdf+Data(pf,4);
Mcdf(pf)=cdf;
end
u=rand(2000000,1);
countu=histc(u,[0,Mcdf]);
cpf=0;
ReSamPara=zeros(length(u),3);
for pf=1:len
ReSamPara(cpf+1:cpf+countu(pf),:)=[rslKs*(rand(countu(pf),1)-0.5)+Data(pf,1),rslJs*(rand(countu(pf),1)-0.5)+Data(pf,2),rslJ0*(rand(countu(pf),1)-0.5)+Data(pf,3)];
cpf=cpf+countu(pf);
end
save('ReSamPara', 'ReSamPara', '-ascii');
