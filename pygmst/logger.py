import logging

logger = logging.getLogger("log")
logger.setLevel(logging.DEBUG)

# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.DEBUG)

# create formatter
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)

import functools


class Log(object):
    def __init__(self, verbose):
        self.logger = logging.getLogger("decorator-log")
        self.verbose = verbose
        # self.logfile = logfile

    def __call__(self, fn):
        @functools.wraps(fn)
        def decorated(*args, **kwargs):
            if self.verbose:
                try:
                    self.logger.info(
                        "{0} - {1} - {2}".format(fn.__name__, args, kwargs)
                    )
                    result = fn(*args, **kwargs)
                    self.logger.info(result)
                    return result
                except Exception as ex:
                    self.logger.info("Exception {0}".format(ex))
                    raise ex
                return result

        return decorated
