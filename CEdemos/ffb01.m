function [e,c]=ReplaceRes(Astar,Q,P,C,y,rPhi0,Phi1,phiVM,phiSP)

B=[rPhi0-Phi1./Astar;phiVM];
b=[feval(Q,Astar*y)*P;C];
c=B\b;
e=phiSP*c;