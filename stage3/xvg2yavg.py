#!/usr/bin/env python
import sys
import os
import re
import matplotlib.pyplot as plt

def plot2d(x, y, xlabel, ylabel, label, pngfile):
    fig = plt.figure(figsize=(12, 8))

    ax = fig.add_subplot(111)
    ax.plot(x, y, label=label)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    fig.savefig(pngfile)

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

    basename = os.path.splitext(os.path.basename(xvgfile))[0]
    pngfile = basename + '.png'
    plot2d(x, y, xlabel, ylabel, title, pngfile)

    yavg = sum(y) / float(len(y))
    print('yavg: ', yavg)

    return yavg

def main():
    xvgfile = sys.argv[1]
    outfile = sys.argv[2]
    yavg = xvg_yavg(xvgfile)
    if len(sys.argv) >= 2:
        with open(outfile, 'w') as f:
            f.write(str(yavg))

if __name__ == "__main__":
    main()
