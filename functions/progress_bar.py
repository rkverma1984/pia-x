import sys
from collections import OrderedDict
import numpy as np


def cli_progress(dict_average, total_frames, metrics, bar_length=20):
    """
    :param total_frames: no of frames considered in the analysis
    :param metrics:
    :param bar_length: length of the progress bar
    :type dict_average: object
    """
    temp = OrderedDict()

    for matrix in metrics:
        temp.update({matrix: []})
    
    for key in dict_average:
        for matrix in dict_average[key].keys():
            temp[matrix].append(len(dict_average[key][matrix]))
       
    progress = []
    for a in temp.keys():
        percent = float(np.sum(temp.get(a))) / total_frames
        hashes = '#' * int(round(percent * bar_length))
        spaces = ' ' * (bar_length - len(hashes))
        progress.append([hashes + spaces + str(int(percent * 100))+'%'])
        
    command = "\rProgress : "
    for i in range(0, len(metrics)):
        command += metrics[i] + ' ' + ' '.join(progress[i]) + ' [|||] '

    sys.stdout.write(command)
    sys.stdout.flush()
