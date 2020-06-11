function [E] = model_energy(Q, x, y, alpha)
eps = 1e-6;
E = zeros(1, length(alpha));
for k=1:length(alpha)
    a=alpha(k);
    i = 2;
    idx = helper(x,y,Q,a,eps,i);
    while (idx>length(x))
        eps=eps*10;
        i = 2;
        idx = helper(x,y,Q,a,eps,i);
    end
    E(k) = x(idx);
end
end
    
function result = helper(x, y, Q, a, eps, i)
while (true)
    if i > length(x)
        break;
    end
    sum_E = trapz(x(1:i), y(1:i))/Q;
    if ( abs(sum_E - a) < eps )
        break;
    end
    i=i+1;
end
result = i;
end
