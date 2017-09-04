% Function file for Cournot oligopoly model.
function fval = exampcournot(q,alpha,beta)
P    = (q(1)+q(2))^(-alpha);
Pder = (-alpha)*(q(1)+q(2))^(-alpha-1);
fval = [P+Pder*q(1)-beta(1)*q(1); ...
        P+Pder*q(2)-beta(2)*q(2)];