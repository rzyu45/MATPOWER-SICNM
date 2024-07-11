function [s,gamma,b,bd,alpha,gammatilde]=coeff_rodas3d
%-- Coefficients for rodas3d
s=4; 
gamma = 0.57281606;
alpha = ...
  [0 0 0 0;
   1.2451051999132263 0 0 0;
   1 0 0 0;
   0.32630307266483527 0.10088086733516474 0.57281606 0];
b = [0.69775271462407906; 0.056490613592447572; -0.32705938821652658; 0.57281606; ...
     ];
beta = ...
  [0 0 0 0;
   -3.1474142698552949 0 0 0;
   0.32630307266483527 0.10088086733516474 0 0;
   0.69775271462407906 0.056490613592447572 -0.32705938821652658 0];
%%
bd=[beta(s-1,1);beta(s-1,2);gamma;0];
gammatilde=beta-alpha;
gammatilde=gammatilde/gamma;
alpha=alpha'; gammatilde=gammatilde';  % order by column
