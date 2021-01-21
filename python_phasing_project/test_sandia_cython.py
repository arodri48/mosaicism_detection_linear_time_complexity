import sandia_stats
import numpy as np


def main():
    x = np.random.normal(10, 5, 1000)
    w1 = x[0:10]
    print(w1.mean())
    print(w1.var())
    w1_sand = sandia_stats.statistical_moment_generator(w1)
    print(w1_sand[0])
    print(w1_sand[1] / w1.size)
    w2 = x[1:11]
    print(w2.mean())
    print(w2.var())
    w2_sandia = sandia_stats.moment_updater(w1_sand, x[0], x[10], w1.size)
    print(w2_sandia[0])
    print(w2_sandia[1] / w1.size)


main()
