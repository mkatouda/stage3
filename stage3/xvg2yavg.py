#!/usr/bin/env python
import os
import re
import argparse
import matplotlib.pyplot as plt

def get_parser():
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description='convert Gromacs xvg file to csv file, draw graph, and get avarage'
    )
    parser.add_argument(
        'xvg', type=str,
        help = 'input Gromacs xvg flie'
    )
    parser.add_argument(
        'out', type=str,
        help = 'output avarage score flie'
    )
    args = parser.parse_args()

    return args

def plot2d(x, y, xlabel, ylabel, label, pngfile):
    fig = plt.figure(figsize=(12, 8))
    ax = fig.add_subplot(111)
    ax.plot(x, y, label=label)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    plt.savefig(pngfile)
    plt.clf()
    plt.close()

def xvg_yavg(xvgfile):
    with open(xvgfile) as fin:
        lines = fin.readlines()

        x = []; y = []
        p = '"([^"]*)"'
        for line in lines:
            if line[0] == '@':
                if 'xaxis' in line:
                    xlabel = re.search(p, line).group().strip('"')
                elif 'yaxis' in line:
                    ylabel = re.search(p, line).group().strip('"')
                elif 'title' in line:
                    title = re.search(p, line).group().strip('"')
            
            elif line[0] != '#' and line[0] != '@':
                list = line.split()
                x.append(float(list[0]))
                y.append(float(list[1]))

    #print(xlabel, ylabel)
    #for i in range(len(x)):
    #    print(x[i], y[i])

    s = '{},{}\n'.format(xlabel, ylabel)
    for i in range(len(x)):
        s += '{},{}\n'.format(x[i], y[i])

    basename = os.path.splitext(os.path.basename(xvgfile))[0]
    csvfile = basename + '.csv'
    with open(csvfile, 'w') as f:
        f.write(s)

    pngfile = basename + '.png'
    plot2d(x, y, xlabel, ylabel, title, pngfile)

    yavg = sum(y) / float(len(y))
    print('yavg: ', yavg)

    return yavg

def main():
    args = get_parser()
    print(args)

    xvgfile = args.xvg
    outfile = args.out
    with open(outfile, 'w') as f:
        f.write(str(xvg_yavg(xvgfile)))

if __name__ == '__main__':
    main()
