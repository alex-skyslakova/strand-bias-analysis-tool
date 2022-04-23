import math
import os


def get_filename(name):
    return os.path.basename(name).split('.')[0]

def to_int(data):
    data = data.apply(lambda x: x.astype(int) if x is not None else math.nan)
    return data


'''ef to_string(my_list):
    return ''.join(my_list)'''


'''def convert(seq, forward):
    if seq in forward:
        return get_reverse_complement(seq)
    return seq'''


'''def get_reverse_complement_count(act_row, data):
    try:
        count = data.loc[act_row['rev_complement']]['seq_count']
        return float(count)
    except KeyError:
        return math.nan'''


'''def is_tool(name):
    """Check whether `name` is on PATH and marked as executable."""

    return which(name) is not None'''