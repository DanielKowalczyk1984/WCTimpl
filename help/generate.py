from random import randint

for j in [3,5,8,10,12]:
    for n in [20,50,100,150]:
        for i in xrange(0,10):
            f = open('instance3_%d_%d_%d.txt'%(n,j,i),'w')
            f.write('p %d\n'%n)
            for it in xrange(0,n):
                f.write('n %d %d\n'%(randint(10,20),randint(10,20)))
