#!/usr/bin/env python

import sys, os, subprocess
import random
import re, string, copy
import time
import math

use_message = '''
Usage:
'''

def draw_bwt_search(match_filename):
    random.seed(0)

    show_steps = [[10 * i for i in range(10)], [10 * i for i in range(10)]]
    # show_steps = [[0, 20, 40, 60, 80], [0, 20, 40, 60, 80]]
    # show_steps = [[0], [0]]
    
    steps = [[], []]
    search_ranges = [[int(math.pow(2, 32)), 0], [int(math.pow(2, 32)), 0]]

    read_begin, read_end = 100000, 0 # inclusive, e.g., [0, 99]

    match_file = open(match_filename, "r")
    for line in match_file:
        dir, start, dep, top, bot, topp, botp = line[:-1].split()
        start, dep, top, bot, topp, botp = \
            int(start), int(dep), int(top), int(bot), int(topp), int(botp)

        if dir == "fw":
            fwi = 0
        else:
            fwi = 1

        assert fwi in range(2)
        fwi_steps = steps[fwi]

        assert start <= len(fwi_steps)
        if start == len(fwi_steps):
            fwi_steps.append([])

        last_steps = fwi_steps[-1]

        step = [top, bot, topp, botp]
        last_steps.append(step)

        search_ranges[fwi][0] = min(top, search_ranges[fwi][0])
        search_ranges[fwi][1] = max(bot, search_ranges[fwi][1])

        if dep < read_begin:
            read_begin = start
        if dep > read_end:
            read_end = dep

    read_len = read_end - read_begin + 1

    graph_width, graph_height = 1600 * 4, 800
    bwt_graph_height = 100
    left_extra, right_extra = 30, 30
    canvas_width, canvas_height = graph_width + left_extra + right_extra, graph_height * 2 + bwt_graph_height + 100

    html_file = open(match_filename + ".html", "w")
    html_prev, html_body, java_script = [], [], []

    html_prev.append(r'<HTML>')
    html_prev.append(r'<HEAD>')
    html_prev.append(r'<TITLE>result</TITLE>')
    html_prev.append(r'<META NAME="description" CONTENT="result">')
    html_prev.append(r'<META NAME="keywords" CONTENT="result">')
    html_prev.append(r'<META NAME="resource-type" CONTENT="document">')
    html_prev.append(r'<META NAME="distribution" CONTENT="global">')
    html_prev.append(r'<META HTTP-EQUIV="Content-Style-Type" CONTENT="text/css">')
    html_prev.append(r'<!--[if IE]><script src="excanvas.js"></script><![endif]-->')
    html_prev.append(r'<style type="text/css">')
    # html_prev.append(r'canvas { border: 1px solid black; }')
    html_prev.append(r'</style>')
    html_body.append(r'</HEAD>')

    html_body.append(r'<BODY text="#000000" bgcolor="#FFFFFF" onload="draw_depth();">')
    html_body.append(r'<canvas id="draw_depth" width="%d" height="%d" name="draw_depth"></canvas>' % \
                         (canvas_width, canvas_height))

    # HTML Canvas Reference - http://www.w3schools.com/tags/ref_canvas.asp
    java_script.append(r'<script type="text/javascript">')
    java_script.append(r"function draw_depth(){")
    # java_script.append(r"alert('message')")
    java_script.append(r"var canvas = document.getElementById('draw_depth');")
    java_script.append(r"if (canvas.getContext){")
    java_script.append(r"var ctx = canvas.getContext('2d');")

    java_script.append(r"ctx.font = '10pt Helvetica';")
    java_script.append(r"ctx.strokeStyle = 'rgba(0, 0, 0, 255)';")

    # http://html-color-codes.info/color-names/
    fill_colors = [[255, 0, 0],\
                       [0, 255, 0], \
                       [0, 0, 255], \
                       [0x8B, 0, 0], \
                       [0xFF, 0x14, 0x93], \
                       [0xFF, 0xA5, 0], \
                       [0xBC, 0x8F, 0x8F], \
                       [0x80, 0, 0x80], \
                       [0, 0x8B, 0x8B], \
                       [0x2F, 0x4F, 0x4F]]
    assert len(steps[0]) == read_len and len(steps[1]) == read_len
    graph_width_diff = graph_width / read_len
    for fwi in range(2):
        search_range = search_ranges[fwi]
        assert search_range[1] > search_range[0]
        search_height = search_range[1] - search_range[0]
        assert search_height > 0
        color_idx = 0
        add_height = 0
        if fwi == 1:
            add_height = graph_height + bwt_graph_height

        fwi_steps = steps[fwi]
        for ebwti in range(2):
            if ebwti == 1:
                break
            
            height_2darray = [[[0, 0] for j in range(read_len)] \
                                  for i in range(len(show_steps[fwi]))]
            for start_idx in range(len(show_steps[fwi])):
                start = show_steps[fwi][start_idx]
                if start >= len(fwi_steps):
                    break

                curr_steps = fwi_steps[start]
                assert len(curr_steps) > 0

                if start == 0:
                    java_script.append(r"ctx.fillStyle = 'rgb(0, 0, 0)';")
                    for i in range(read_len):
                        text_x = (start + i) * graph_width_diff + left_extra - 7
                        text_y = graph_height + add_height + 50
                        java_script.append(r"ctx.fillText('%d', %d, %d);" % (i + 1, text_x, text_y))

                curr_x = start * graph_width_diff + left_extra

                fill_color = fill_colors[color_idx % len(fill_colors)]
                color_idx += 1
                java_script.append(r"ctx.strokeStyle = 'rgba(%d, %d, %d, 255)';" % \
                                       (fill_color[0], fill_color[1], fill_color[2]))
                java_script.append(r"ctx.fillStyle = 'rgba(%d, %d, %d, %f)';" % \
                                       (fill_color[0], fill_color[1], fill_color[2], 0.5))

                for dep in range(len(curr_steps) - 1):
                    top, bot = curr_steps[dep][2 * ebwti:2 + 2 * ebwti]
                    next_top, next_bot = curr_steps[dep+1][2 * ebwti:2 + 2 * ebwti]

                    top -= search_range[0]
                    bot -= search_range[0]
                    next_top -= search_range[0]
                    next_bot -= search_range[0]

                    curr_y = int(top / float(search_height) * graph_height) + add_height
                    next_x = curr_x + graph_width_diff
                    next_y = int(next_top / float(search_height) * graph_height) + add_height

                    height = min(bot - top, 50)
                    next_height = min(next_bot - next_top, 50)

                    curr_y_bot = curr_y + height
                    next_y_bot = next_y + next_height

                    read_pos = start + dep
                    assert start_idx < len(height_2darray) and read_pos + 1 < len(height_2darray[start_idx])
                    if start_idx > 0:
                        p_start = show_steps[fwi][start_idx - 1]
                        last_steps = fwi_steps[p_start]

                        if read_pos - p_start + 1 < len(last_steps):
                            p_top, p_bot = last_steps[read_pos - p_start][2 * ebwti:2 + 2 * ebwti]
                            p_next_top, p_next_bot = last_steps[read_pos - p_start + 1][2 * ebwti:2 + 2 * ebwti]
                            p_curr_y, p_curr_y_bot = height_2darray[start_idx - 1][read_pos]
                            p_next_y, p_next_y_bot = height_2darray[start_idx - 1][read_pos + 1]

                            if top <= p_top and bot >= p_bot:
                                curr_y = min(curr_y, p_curr_y)
                                curr_y_bot = max(curr_y_bot, p_curr_y_bot)

                            if next_top <= p_next_top and next_bot >= p_next_bot:
                                next_y = min(next_y, p_next_y)
                                next_y_bot = max(next_y_bot, p_next_y_bot)                    

                    height_2darray[start_idx][read_pos] = [curr_y, curr_y_bot]
                    height_2darray[start_idx][read_pos + 1] = [next_y, next_y_bot]

                    java_script.append(r"ctx.beginPath ();")
                    java_script.append(r"ctx.moveTo (%d, %d);" % (curr_x, curr_y))
                    java_script.append(r"ctx.lineTo (%d, %d);" % (next_x, next_y))

                    if height > 1:
                        java_script.append(r"ctx.lineTo (%d, %d);" % (next_x, next_y_bot))
                        java_script.append(r"ctx.lineTo (%d, %d);" % (curr_x, curr_y_bot))

                    java_script.append(r"ctx.closePath ();")
                    java_script.append(r"ctx.stroke ();")
                    java_script.append(r"ctx.fill ();")

                    if dep + 2 == len(curr_steps):
                        # java_script.append(r"ctx.fillStyle = 'rgb(0, 0, 0)';")
                        java_script.append(r"ctx.fillText('%dX', %d, %d);" % (start + 1, next_x, next_y + start * 1.5))

                    curr_x = next_x


    java_script.append(r"}")
    java_script.append(r"}")
    java_script.append(r'</script>')

    html_body.append(r'</BODY>')
    html_body.append(r'</HTML>')
    
    match_file.close()

    for line in html_prev + java_script + html_body:
        print >> html_file, line

    html_file.close()
    

if __name__ == "__main__":
    if len(sys.argv) == 2:
        draw_bwt_search(sys.argv[1])
    else:
        print use_message
