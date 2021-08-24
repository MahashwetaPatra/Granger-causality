function X_em = cstr_em(A, m, tau)

X_em = zeros(length(A),m);

for t=m*tau:length(A)
    X_em(t,:)=A((t-(m-1)*tau):tau:t);
end