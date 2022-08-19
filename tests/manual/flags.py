import sys

FLAGS = (
    (1, 'paired'),
    (2, 'proper pair'),
    (4, 'read unmapped'),
    (8, 'mate unmapped'),
    (16, 'read reverse strand'),
    (32, 'mate reverse strand'),
    (64, 'is read 1'),
    (128, 'is read 2'),
    (256, 'secondary alignment'),
    (512, 'fail QC'),
    (1024, 'duplicate'),
    (2048, 'supplementary alignment')
)

for line in sys.stdin:
    flag = int(line.strip())
    annotations = ', '.join(label for fl, label in FLAGS if flag & fl)
    print(f'{flag:4}\t{annotations}')
