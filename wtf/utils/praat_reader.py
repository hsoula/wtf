"""
    : praat_reader.py
    Created: 18/03/16
    Description: simple reader for praat TxtGrid output
    Should output event in a standard form : type start end label duration
"""
import sys


def read_tag_line(lines, ix):
    """
        Read a whole tag line (line that starts with !) using the idx
        and return
    :param lines: the file in lines
    :param ix: line ix
    :return: type start end label duration
    """

    # first recheck its a tag line
    line = lines[ix]
#    print line
    if len(line) == 0 or line[0] != '!':
        return None
    tag_line = line.split()
    tag_type = tag_line[1][:-1]
    time_line = lines[ix+1].split()
    #print time_line
    if len(time_line) == 3:
        start = float(time_line[1])
        end = float(time_line[2])
    if len(time_line) == 2:
        start = float(time_line[1])
        end = start

    duration = end - start
    label_line = lines[ix+2].split()
    label = label_line[0].strip("\"")
    if len(label) == 0:
        return None

    return tag_type, start, end, label, duration


def read_txtgrid(fname):
    res = []
    with open(fname) as iOF:
        llines = iOF.readlines()
        lines = [line.strip() for line in llines]
        nl = len(lines)
        for ix in range(nl):
            ret = read_tag_line(lines, ix)
            if ret is not None:
             res.append(ret)
    return res

if __name__=='__main__':
    fname = sys.argv[1]
    res = read_txtgrid(fname)
    for l in res:
        print l[0], l[1], l[2], l[3], l[4]


