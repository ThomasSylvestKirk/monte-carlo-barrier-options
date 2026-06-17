clc
clear all
seed=1001;
rng(seed) % save seed number
dim=250000; % dim is repetetions
n=254; %Monitoring points to maturity
T=0.5; %how long time to maturity in trading days 
h=T/n; %steplength
r=0.1; 
S0=100; 
K=95;
D=95;
sigma=0.1;
Z=normrnd(0,1,n+1,dim);

%% Check if function works
Settle='01-Jan-2015';
Maturity='01-Jul-2015';
RateSpec = intenvset('ValuationDate', Settle, 'StartDates', Settle, 'EndDates', Maturity, ...
'Rates', r, 'Compounding', -1, 'Basis', 1);
StockSpec=stockspec(sigma,S0);
Strike = K;
OptSpec = 'call';
Barrier = D;
BarrierSpec = 'DO';
Price = barrierbybls(RateSpec, StockSpec, OptSpec, Strike, Settle,...
Maturity,  BarrierSpec, Barrier)
Price= DOCbls(S0,T,K,r,sigma,D)

%% Naive comparasion of continuous value and the discrete estimate
clear EULERMbias
%Here I create a plot and a table descrbining how the bias decreases as n increases
Pricemonitorn=DOCbls(S0,T,K,r,sigma,D);
Price=DOCbls(S0,T,K,r,sigma,D);

u=0;
N=1:1:250;
for j=N
    u=u+1;
    [Estimate(u),SE(u),HW(u)]=BarrieroptionEulerNAIV(S0,D,K,sigma,T,j,dim,r,seed);
end

%This is a table of the figure above
figure
plot(N,100*(Estimate-Price)/Price,'k-')

title('Figure 1. Bias of the Naiv estimator')
xlabel('n')
ylabel('E(X)^n-E(X) (%)')
xlim([1 250])
ylim([0 11])
hold on
set(gca,'color',[0.901960784313726 0.901960784313726 0.901960784313726])

%% Table 1
%Barrier option comparasion of ABAST and comntinuous value perhaps also
%Euler low monitoring frequency.
%Comparasion of Barrier options computed using a control variate and
%nothing else and using a brownian bridge. Table 1
n=300
Z=normrnd(0,1,n+1,dim);
B=[]; A=[];C=[];E=[];F=[];G=[];EE=[]; T=1;
S0=100; M=0; r=0.1; N=[]; SIGMA=[]; DD=[]; 

JJ=[];
sigma1=0.1:0.2:0.3; n1=[10 50 250]; D1=[90 99];
for sigma=sigma1
    for D=D1
    for n=n1
        M=M+1;
        tic
[EST2,SE2,HW2]=DOCBB(S0,D,K,sigma,T,1,dim,r,Z(1,1:dim));
toc
tic
[EST3,SE3,HW3]=BarrieroptionEulerNAIVrunswithout(S0,D,K,sigma,T,n,dim,r,Z(1:n,:));
toc
tic
[EST4,SE4,HW4]=BarrieroptionABAST(S0,D,K,sigma,T,n,dim,r,Z(1:n,:));
toc
B=[B;EST2,HW2];
A=[A;EST3,HW3];
G=[G;EST4,HW4];
JJ(M)=DOCbls(S0,T,K,r,sigma,D);
SIGMA=[SIGMA; sigma];
N=[N; n];
DD=[DD;D];
    end
    end
end

% Compute error of the estimated continuous option
JJ=JJ'; 
BBB=B(:,1); 
AAA=A(:,1); GGG=G(:,1);
for i=1:length(AAA)
BBB(i,1)=((BBB(i,1)-JJ(i))/JJ(i))*100;
AAA(i,1)=((AAA(i,1)-JJ(i))/JJ(i))*100;
GGG(i,1)=((GGG(i,1)-JJ(i))/JJ(i))*100;
end
Tabel=array2table([round(BBB,2),round(AAA,2),round(GGG,2),round(JJ,2)],'VariableNames',{'BB error (%)' 'Simple error (%)' 'ABAST error (%)' 'Cont. Price'});
Tabel.D=DD;
Tabel.sigma=SIGMA;
Tabel.n=N;
Tabel=[Tabel(:,end-2:end), Tabel(:,1:end-3)];
Tabel=sortrows(Tabel,[1 2 3]);
disp(Tabel)

%% Simulation of discrete option. How large VRF is there of using a control variate
%Table 2 and 3
dim=250000; 
n=254;
Z=normrnd(0,1,n+1,dim);
sigma=0.1;
CPUTIME=[]; sigma=0.1; Preciseest=[];
S0=100; K=90; T=1;
D=90;
    for sigma=0.1:0.2:0.5
        for n=[12 52 253]
        tic
        [EST2,SE1,HW1] = CVforDOCCONT(S0,D,K,sigma,T,n,dim,r,Z(1:n,1:dim));
        TimeCV=toc;
        tic
        [EST,SE3,HW] = CVEURUND(S0,D,K,sigma,T,n,dim,r,Z(1:n,1:dim));
        TimeCV2=toc;
        tic
        [Estimate,SE,~] = BarrieroptionEulerNAIVrunswithout(S0,D,K,sigma,T,n,dim,r,Z(1:n,1:dim));
        TimeSimple=toc;
        CPUTIME=[CPUTIME; [sigma,n,round(EST2,2),round(HW1,2),round((SE3.^2)/(SE1.^2),2),round((TimeCV2/TimeCV)^-1,2),round((SE.^2)/(SE1.^2),2),round((TimeSimple/TimeCV)^-1,2),round(DOCbls(S0,T,K,r,sigma,D),2)]];
   Preciseest=[EST2;Preciseest];
        end
    end
    Preciseest
    Preciseest=[];
CPUTIME=array2table(CPUTIME,'VariableNames',{'σ' 'n' 'Price' 'HW' 'VRF CV 2' 'TF CV 2' 'VRF Simple' 'TF Simple' 'Cont. Price'});
disp(CPUTIME);
sigma=0.1; T=1; CPUTIME=[];
    for D=85:5:95
        for n=[12 50 253]
        tic
        [EEE,SE3,HW] = CVEURUND(S0,D,K,sigma,T,n,dim,r,Z(1:n,:));
        M=toc;
        tic
        [EST,SE2,HW1] = CVforDOCCONT(S0,D,K,sigma,T,n,dim,r,Z(1:n,:));
        M2=toc;
        tic
        [Estimate,SE,~] = BarrieroptionEulerNAIVrunswithout(S0,D,K,sigma,T,n,dim,r,Z(1:n,:));
        TimeSimple=toc;
        Preciseest=[EST;Preciseest];
   
        CPUTIME=[CPUTIME; [D,n,round(EST,2),round(HW1,2),round((SE3.^2)/(SE2.^2),2),round((M/M2)^-1,2),round((SE.^2)/(SE2.^2),2),round((TimeSimple/M2)^-1,2),round(DOCbls(S0,T,K,r,sigma,D),2)]];
        end
    end
Preciseest
CPUTIME=array2table(round(CPUTIME,2),'VariableNames',{'D' 'n' 'Estimate' 'HW' 'VRF CV 2' 'TF CV 2' 'VRF Simple' 'TF Simple' 'Cont. Price'});
disp(CPUTIME);

%% Finite difference method central difference
% Figure 3 and 4
S0=100;
sigma=0.1; T=1; r=0.1; K=100; S0=100; D=95;

dim=50000;
Z=normrnd(0,1,n+1,dim);
DELTA=[]; VEGA=[]; DELTA1=[]; VEGA1=[]; jj=0;
sigma=0.1; T=1; r=0.1; K=100; S0=100; sigma=0.1
D1=80:1:99; K1=80:4:120;

    for K=K1
jj=jj+1;
ii=0;
        for D=D1
            ii=1+ii;
[DeltaC(jj,ii),DeltavarC]=CentraldifferencesBB(S0,0.0003,0,D,K,sigma,T,1,dim,r,Z);
[VegaC(jj,ii),VegavarC]=CentraldifferencesBB(S0,0,0.00004,D,K,sigma,T,1,dim,r,Z);
DELTA1=[DELTA;DeltaC];
VEGA1=[VEGA;VegaC];
        end
            DELTA=[DELTA;DELTA1];
            VEGA=[VEGA;VEGA1];
    end
    
[Sigma,DD]=meshgrid(K1,D1);
figure
surf(DD,Sigma,DeltaC')
title('Figur 3. Delta sensitivity surface')
xlabel('D')
ylabel('K')
zlabel('Delta')
figure
surf(DD,Sigma,VegaC')
title('Figur 4. Vega sensitivity surface')
xlabel('D')
ylabel('K')
zlabel('Vega')

%% Estimating the coefficients of control variate regression + underlying asset to show possible multicollinarity
%Figure 5
dim=250000; 
n=253; Z=normrnd(0,1,n+1,dim); 
A=[];
D=90; S=100; r=0.1; T=1; sigma=0.1;
for n=[12 121 253]
for D=85:10:95
    for K=[80 100]
[~,~,~,b] = CV3forDOCCONT(S0,D,K,sigma,T,n,dim,r,Z(1:n,1:dim));
A=[D,K,n,b';A];
    end
end
end

B=A(:,[1,2,3,5,6,7]);
B=sortrows([B],[3 2 1]);
%D
Legend1=num2str(B(:,1));
%K
Legend2=num2str(B(:,2));
%n
Legend3=num2str(B(:,3));

B=B(:,[4,5,6]);
figure
plot([1:length(Legend1)]',B(:,1),'K-o')
hold on
plot([1:length(Legend1)]',B(:,2),'K-+')
plot([1:length(Legend1)]',B(:,3),'K-*')
set(gca,'color',[0.901960784313726 0.901960784313726 0.901960784313726]);

Legend=[]; Legends=[];
for i=1:length(Legend1)
Legend=strcat(2,num2str(Legend1(i,:)),'/',num2str(Legend2(i,:)),'/',num2str(Legend3(i,:)));
Legends=[Legends;Legend];
end
legend('Continuously monitored DOC option','European call option','The underlying asset','Position',[0.726519065069579 0.114761418012897 0.169661461586753 0.102831198873683],'FontSize',16)
xlim([1,length(Legends)]);
clear xticks
title('Figur 5. regression coefficients','FontSize',20)
xticks([1:27]');
xticklabels({Legends});
xlabel('D / K / n','FontSize',18);
ylabel('b','FontSize',18)

%% Estimating the correlation of control variate regression of only European call option and the continuous DOC option. Plot of correlation
sigma=0.1; 
S0=100; K=95; T=1; r=0.1;
A=[];
A1=[];
for n=[52 253]
for D=[95 99]
    for K=[90 100 110]
[~,~,~,~,b,rho2] = CVforDOCCONTrho(S0,D,K,sigma,T,n,dim,r,Z(1:n,1:dim));
A1=[D,n,K,rho2;A1];
end
end
end

A=[];
for n=[52 253]
for D=[95 99]
    for K=[90 100 110]
[~,~,~,b,Rsquared] = CVEURUND(S0,D,K,sigma,T,n,dim,r,Z(1:n,1:dim));
A=[D,n,K,Rsquared;A];
end
end
end

B1=A1;
B1=sortrows([B1],[1 2 3]);
%D
Legend1=num2str(B1(:,1));
%n
Legend2=num2str(B1(:,2));
B1=B1(:,4);

B=A;
B=sortrows([B],[1 2 3]);
%D
Legend1=num2str(B(:,1));
%n
Legend2=num2str(B(:,2));

Legend3=num2str(B(:,3));

B=B(:,4);

figure
plot([1:length(Legend1)]',B(:,1),'K-+')
hold on
plot([1:length(Legend1)]',B1(:,1),'K-o')
set(gca,'color',[0.901960784313726 0.901960784313726 0.901960784313726])

Legend=[]; Legends=[];
for i=1:length(Legend1)
Legend=strcat(2,num2str(Legend1(i,:)),'/',num2str(Legend2(i,:)),'/',num2str(Legend3(i,:)));
Legends=[Legends;Legend];
end

legend({'2 CV','Continuously monitored DOC option'},'FontSize',12,'Location','southwest')
xlim([1,length(Legends)]);
clear xticks
xticks([1:27]');
xticklabels({Legends});
xlabel('D / n / K','FontSize',16);
ylabel('\rho^2 / Rsquared','FontSize',16)
title('Figur 2. Squared correlation','FontSize',24)

%% Plots of payoff of discretely and continuously monitored
K=95;
D=95;
sigma=0.1;
[Y,X] = DOC3CV(S0,D,K,sigma,T,n,dim,r,Z(1:n,1:dim));
figure 
scatter(X(2,:),Y,'MarkerEdgeColor','none','LineWidth',1e-05,...
    'MarkerFaceColor',[0 0 0])
xlabel('Continuously monitored DOC option','FontSize',14);
ylabel('Discretely monitored DOC option','FontSize',14)
set(gca,'color',[0.901960784313726 0.901960784313726 0.901960784313726])
title('Figure 6. Scatter plot of payoff','FontSize',14)
figure
scatter(X(3,:),Y,'MarkerEdgeColor','none','LineWidth',1e-05,...
    'MarkerFaceColor',[0 0 0])
xlabel('European call option','FontSize',14);
ylabel('Discretely monitored DOC option','FontSize',14)
set(gca,'color',[0.901960784313726 0.901960784313726 0.901960784313726])
title('Figure 7. Scatter plot of payoff','FontSize',14)
figure
scatter(X(4,:),Y,'MarkerEdgeColor','none','LineWidth',1e-05,...
    'MarkerFaceColor',[0 0 0])
xlabel('Underlying asset','FontSize',14);
ylabel('Discretely monitored DOC option','FontSize',14)
set(gca,'color',[0.901960784313726 0.901960784313726 0.901960784313726])
title('Figure 8. Scatter plot of payoff','FontSize',14)
