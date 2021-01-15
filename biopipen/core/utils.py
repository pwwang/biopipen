import logging

from rich.logging import RichHandler as _RichHandler
from rich.console import Console
from rich.text import Text

from .defaults import PKG_NAME

class RichHandler(_RichHandler):
    """Subclass of rich.logging.RichHandler, showing log levels as a single
    character"""
    def get_level_text(self, record: logging.LogRecord) -> Text:
        """Get the level name from the record.
        Args:
            record (LogRecord): LogRecord instance.
        Returns:
            Text: A tuple of the style and level name.
        """
        level_name = record.levelname
        level_text = Text.styled(
            level_name[:5].upper(), f"logging.level.{level_name.lower()}"
        )
        return level_text


_logger_handler = RichHandler(show_path=False,
                              show_level=True,
                              console=Console(),
                              rich_tracebacks=True,
                              markup=True)

logger = logging.getLogger(PKG_NAME)
logger.addHandler(_logger_handler)
logger.setLevel(logging.DEBUG)
