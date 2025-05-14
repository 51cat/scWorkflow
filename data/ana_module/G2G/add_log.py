import logging
from functools import wraps
import time
from datetime import timedelta
import sys

def add_log(func):
    logFormatter = logging.Formatter(
        '[%(levelname)s][%(asctime)s] %(name)s %(message)s'
    )

    # 创建一个日志器
    logger = logging.getLogger(
        f'{func.__module__}.{func.__name__}'
    )
    
    # 设置日志级别
    logger.setLevel(logging.INFO)

    # 创建一个控制台处理器
    consoleHandler = logging.StreamHandler(sys.stderr)
    consoleHandler.setFormatter(logFormatter)
    logger.addHandler(consoleHandler)

    @wraps(func)
    def wrapper(*args, **kwargs):
        logger.info('start...')
        start = time.time()

        result = func(*args, **kwargs)
    
        end = time.time()
        used = timedelta(seconds=end - start)
        logger.info('done. time used: %s', used)
        
        return result

    wrapper.logger = logger
    return wrapper