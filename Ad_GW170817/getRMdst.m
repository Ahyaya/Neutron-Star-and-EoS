Z=0;
Cent{1}=linspace(8.8,17.5,42);
Cent{2}=linspace(0.8,2.7,32);

fp=fopen("InteralData.txt","r");

for pf=1:221
D=fscanf(fp,"%lf",[4,500000]);
D=D';D=D(:,1:2);
tempZ=hist3(D,"Ctrs",Cent);
Z=tempZ+Z;
end

fclose(fp);
clear D,tempZ;

Z=Z';
save WPout_RMdst;
