parool('local', 2)

ncwrite('test.nc', 'test', ones(500.*labindex, 100), [1+500.*(labindex-1) 1], [1 1])