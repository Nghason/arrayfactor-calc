function coflows_sorted=  sortCoflowsByID(coflows)
[~, ind] = sort([coflows.id]);
coflows_sorted = coflows(ind);
end