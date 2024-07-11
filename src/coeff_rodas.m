function [s,gamma,b,bd,a,alpha,gammatilde,g,c,d]=coeff_rodas
%-- Coefficients for RODASP with order 4 for linear parabolic problems
s=6; alpha=zeros(s,s); beta=alpha; a=zeros(s,1); g=zeros(s,1);
gamma=0.25;
alpha(2,1)=0.75;
alpha(3,1)= 8.6120400814152190E-2; alpha(3,2)=0.1238795991858478;
alpha(4,1)= 0.7749345355073236;    alpha(4,2) = 0.1492651549508680; alpha(4,3) = -0.2941996904581916;
alpha(5,1)= 5.308746682646142;     alpha(5,2) = 1.330892140037269;  alpha(5,3) = -5.374137811655562; 
alpha(5,4)= -0.2655010110278497;
alpha(6,1)= -1.764437648774483;    alpha(6,2) =-0.4747565572063027; alpha(6,3) = 2.369691846915802;
alpha(6,4)= 0.6195023590649829;    alpha(6,5) = 0.25;
beta(2,1)= 0.0; 
beta(3,1)= -0.049392;            beta(3,2)= -0.014112;
beta(4,1)= -0.4820494693877561;  beta(4,2)= -0.1008795555555556;    beta(4,3)=  0.9267290249433117;
beta(5,1)= -1.764437648774483;   beta(5,2)= -0.4747565572063027;    beta(5,3) =  2.369691846915802;
beta(5,4)=  0.6195023590649829; 
beta(6,1)= -8.0368370789113464E-2;beta(6,2)=-5.6490613592447572E-2; beta(6,3)=  0.4882856300427991;
beta(6,4)=  0.5057162114816189;   beta(6,5)= -0.1071428571428569;
b=[beta(6,1);beta(6,2);beta(6,3);beta(6,4);beta(6,5);gamma];
bd=[beta(5,1);beta(5,2);beta(5,3);beta(5,4);gamma;0];
c=[-40.98639964388325;-10.36668980524365;44.66751816647147;4.13001572709988;2.55555555555556;0];
d=[73.75018659483291;18.54063799119389;-81.66902074619779;-6.84402606205123;-3.77777777777778;0];
gammatilde=beta-alpha;
for i=1:s
 a(i)=sum(alpha(i,:));
 g(i)=sum(gammatilde(i,:))+gamma;
end
gammatilde=gammatilde/gamma;
alpha=alpha'; gammatilde=gammatilde';  % order by column