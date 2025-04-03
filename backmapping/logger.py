import logging
import os


# TODO: move to utils
def check_if_file_exists(path: str):
    if not os.path.isfile(path):
        logger.error(f"Couldn't find file: {path}")
        return False
    return True


LOG_FILE = "backmapping.log"

logging.basicConfig(
    level=logging.DEBUG,
    format="%(asctime)s - %(levelname)s - %(filename)s - %(message)s",
    handlers=[logging.FileHandler(LOG_FILE, mode="w"), logging.StreamHandler()],
)

logger = logging.getLogger(__name__)
