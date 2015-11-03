from random import randint

for i in xrange(1,10):
    f = open('test%d.txt'%i,'w')
    f.write('p 75\n')
    for it in xrange(1,76):
        f.write('n %d %d\n'%(randint(1,20),randint(1,100)))
