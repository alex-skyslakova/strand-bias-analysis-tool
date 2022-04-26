import math
import os
import dateutil.rrule as rrule
from dateutil.parser import parse as dparse


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


def is_or_create_dir(dir):
    if not os.path.isdir(dir):
        os.makedirs(dir)


def get_n_percent(df, n, tail=False):
    if tail:
        return df.tail(int(len(df) * (n / 100)))
    else:
        return df.head(int(len(df) * (n / 100)))


def filter_time(self, descr, time_from, time_to):
    # extract start time from description
    time = dparse([i for i in descr.split() if i.startswith('start_time')][0].split('=')[1])
    if dparse(time_from) <= time <= dparse(time_to):
        return True
    else:
        return False


def hours_aligned(self, start, end, interval=1):
    chunks = []
    rule = rrule.rrule(rrule.HOURLY, byminute=0, bysecond=0, dtstart=start, interval=interval)
    for x in rule.between(start, end, inc=False):
        chunks.append(x)
    chunks.append(end)
    return chunks
