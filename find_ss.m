function ID_ss = find_ss(x_pr,X_em,m,l_train,num_nb,tau)

d = zeros(l_train,1);
for t = m*tau:l_train
    d(t) = norm(x_pr-X_em(t,:));
end

[~,id] = sort(d,'ascend');
ID_ss = id(m*tau:(m*tau+num_nb-1));
