from random import seed
from random import randint
from random import choice
from time import time
import sys

# seed random number generator
seed(1)

n = 59468446
cigars = []
for _ in range(n):
    chrom = randint(1, 24)
    if chrom == 23:
        chrom = 'X'
    elif chrom == 24:
        chrom == 'Y'

    start = randint(1, 4294967295)
    stop = randint(start, 4294967295)
    pos = randint(0, 1000)

    nodes = []
    n_nodes = randint(1, 3)
    for i in range(n_nodes):
        ops = []
        n_ops = randint(1, 5)
        for _ in range(n_ops):
            num = randint(1, 150)
            op = choice(['M', 'X', 'I', 'D', 'S'])
            ops.append('{}{}'.format(num, op))            
        nodes.append('{}[{}]'.format(i, ''.join(ops)))

    cigars.append('chr{}_{}_{},{},{}'.format(chrom, start, stop, pos, ''.join(nodes)))

    # chr1_31555_31570,177,0[150M]

print('Done generating. Sorting...')
print('Size of list: {}'.format(sys.getsizeof(cigars)))

for i in range(0, 20):
    n = 59468446 // 20 * (i + 1)
    t0 = time()
    cigars[:n].sort()
    print("list.sort: %.3fs" % (time()-t0))