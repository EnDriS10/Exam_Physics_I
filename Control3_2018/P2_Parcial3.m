clear;clc;


%% PROBAR EL CODDIGO

%  k0= [1 ;1 ;2];
%  D = [1; 1 ;1];
%  T0 =[25; 25; 25];
%  T =20;
% 
% l= [3; 2; 5];
% m=[1 1];
% L=15;
% % Define the equilibrium function for the springs
% d = equilibriom_3_springs(k0, D, T0, l, m, L, T);

%%
function [d]= equilibriom_3_springs(k0,D,T0,l,m,L,T)

    g=9.8;
    k= get_k(k0,D,T0, T); ka=k(1); kb=k(2); kc=k(3);
    m1=m(1); m2=m(2); la= l(1); lb=l(2); lc=l(3);

    A = [ka -kb 0 ; 0 kb kc ; 1 1 1]
    b = [ g*m1 ; g*m2 ; L-la-lb-lc ] 

    d =inv(A)*b;

end

function [k]= get_k(k0,D,T0, T)
    k= k0 - (D)./( exp( T0./T ) - 1 );
end

