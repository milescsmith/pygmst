import logging
import coloredlogs
from typing import Optional


def setup_logging(name: Optional[str] = None):
    coloredlogs.DEFAULT_FIELD_STYLES = {
        'asctime': {'color': 'green'},
        'levelname': {'bold': True, 'color': 'red'},
        'module': {'color': 73},
        'funcName': {'color': 74},
        'lineno': {'bold': True, 'color': 75},
        'message': {'color': 'yellow'}
        }
    coloredlogs.install(level="DEBUG", fmt='[%(asctime)s] {%(module)s:%(funcName)s():%(lineno)d} %(levelname)s - %(message)s')
    logger = logging.getLogger("sqanti3_qc")
    logger.setLevel(logging.DEBUG)
    logger.propagate = True

    if name:
        fh = logging.FileHandler(filename=name)
    else:
        fh = logging.FileHandler(filename=f"{__name__}.log")
    formatter = logging.Formatter('[%(asctime)s] {%(module)s:%(funcName)s:%(lineno)d} %(levelname)s - %(message)s')
    fh.setLevel(logging.DEBUG)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
