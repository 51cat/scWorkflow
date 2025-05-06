import logging
from functools import wraps
import time
from datetime import timedelta
import sys
from rich.logging import RichHandler

def add_log(func):
    log_formatter = logging.Formatter("%(message)s", datefmt="[%X]")

    logger = logging.getLogger(f'{func.__module__}')
    logger.setLevel(logging.INFO)

    # 避免重复添加 handler
    if not logger.hasHandlers():
        rich_handler = RichHandler(
            rich_tracebacks=True,
            markup=True,
            show_time=True,
            show_path=False
        )
        rich_handler.setFormatter(log_formatter)
        logger.addHandler(rich_handler)

    @wraps(func)
    def wrapper(*args, **kwargs):
        logger.info("[bold green]Start[/bold green] 🚀")
        start = time.time()

        result = func(*args, **kwargs)

        end = time.time()
        used = timedelta(seconds=end - start)
        logger.info(f"[bold blue]Done[/bold blue] ✅ Time used: {used}")

        return result

    wrapper.logger = logger
    return wrapper

