# -*- coding: utf-8 -*-

__author__ = 'anna'
import colorsys

import matplotlib.pyplot as plt


def not_constant(ys):
    return len(set(ys)) > 1


def create_plot(xs, ys_list, legend, x_label, y_label, title, log_x_scale=False, log_y_scale=False, plot_num=121):
    initialise_fig((12, 6), 10)
    create_subplot(xs, ys_list, legend, x_label, y_label, title, log_x_scale, log_y_scale, plot_num)
    plt.show()


def create_subplot(xs, ys_list, legend, x_label, y_label, title, log_x_scale=False, log_y_scale=False, plot_num=131,
                   legend_loc='center left', bb_to_anchor=(1, 0.5), num_colors=0):
    plt.subplot(plot_num)
    n = num_colors if num_colors else len(legend)
    # colors = (colorsys.hsv_to_rgb(x * 1.0 / n, 0.6 + 0.4 * random(), 0.65 + 0.1 * random()) for x in range(n))
    colors = (colorsys.hsv_to_rgb(0.5, 0, 0.05 + 0.6 * x / n) for x in range(n))
    line_styles = ['-', '--', ':']
    line_widths = [1.0, 1.5, 2.0]
    markers = ["v", "o", "s", "*", "p", "8", "h", "^", "<", ">", "1", "2", "d", "3", "4", "H", "+", "D", "|"]
    i = 0
    different_xs = False
    if xs and isinstance(xs[0], list):
        different_xs = True
        xs = iter(xs)
    for ys in ys_list:
        plt.plot(next(xs) if different_xs else xs, ys, color=next(colors), linewidth=line_widths[i % len(line_widths)],
                 ls=line_styles[i % len(line_styles)], marker=markers[i % len(markers)], markersize=8)
        i += 1
    if legend:
        # plt.rc('font', **{'sans-serif': 'Arial', 'family': 'sans-serif'})
        plt.legend(legend, fontsize="small", loc=legend_loc, bbox_to_anchor=bb_to_anchor, fancybox=True,
                   framealpha=True, frameon=True, numpoints=1)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    if log_x_scale:
        plt.xscale('log', basex=2)
    if log_y_scale:
        plt.yscale('log', basey=2)
    plt.title(title)


def initialise_fig(fig_size, font_size):
    plt.rcParams.update({'font.size': font_size})
    plt.figure(tight_layout=False, figsize=fig_size)


def save_fig(path):
    plt.savefig(path)
    plt.close()


def create_figure(xs, legend2ys, x_label, path, legend2color=None):
    if not legend2color:
        legend2color = {}
    plt.rcParams.update({'font.size': 10})
    plt.figure(tight_layout=False, figsize=(9, 20))
    i = 1
    n = len(legend2ys)
    m = 1
    for legend, ys in legend2ys.iteritems():
        plt.subplot(n, m, i)
        if legend2color and legend in legend2color:
            plt.plot(xs, ys, linewidth=1.0, ls='-', color=legend2color[legend])
        else:
            plt.plot(xs, ys, linewidth=1.0, ls='-')
        plt.title(legend)
        i += 1
    plt.xlabel(x_label)
    plt.subplots_adjust(hspace=1.0)
    save_fig(path)
