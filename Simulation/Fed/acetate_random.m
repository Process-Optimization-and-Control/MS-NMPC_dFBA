nfe = 600;

ac_r_t = ones(10,nfe);

for i = 2:9
ac_r_t(i,1:(i-1)*(nfe/10)) = 0.0*ac_r_t(i,1:(i-1)*(nfe/10));
ac_r_t(i,:) = ac_r_t(i,randperm(length(ac_r_t(i,:))));
sum(ac_r_t(i,:))/nfe
end
ac_r_t(10,:) = 0*ac_r_t(10,:);

csvwrite('ac_id_600.csv',ac_r_t)