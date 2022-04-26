import math
import os


def get_filename(name):
    return os.path.basename(name).split('.')[0]

def to_int(data):
    data = data.apply(lambda x: x.astype(int) if x is not None else math.nan)
    return data

def unique_path(f):
    fnew = f
    root, ext = os.path.splitext(f)
    i = 0
    while os.path.exists(fnew):
        i += 1
        fnew = '%s_%i%s' % (root, i, ext)
    return fnew

'''def get_reverse_complement_count(act_row, data):
    try:
        count = data.loc[act_row['rev_complement']]['seq_count']
        return float(count)
    except KeyError:
        return math.nan'''


'''def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""

    return which(name) is not None'''