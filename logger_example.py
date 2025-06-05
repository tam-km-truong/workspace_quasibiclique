import logging
import sys

# Custom verbosity levels
SILENT = 50      # 0 - Critical only
MINIMAL = 30     # 1 - Warning+
VERBOSE = 20     # 2 - Info+ 
DEBUG = 10       # 3 - Debug+

def setup_logger(name: str, verbosity_level=1, log_file=None):
    """
    Setup logger with verbosity levels.
    verbosity_level: 0=silent, 1=minimal, 2=verbose, 3=debug
    """
    level_map = {
        0: SILENT, 1: MINIMAL, 2: VERBOSE, 3: DEBUG
    }
    
    level = level_map.get(verbosity_level, VERBOSE)
    logger = logging.getLogger(name)
    logger.setLevel(level)
    
    # Console handler (uses stdout for redirection with >)
    console_handler = logging.StreamHandler(sys.stdout)
    # Simple format
    formatter = logging.Formatter('%(levelname)s - %(message)s')
    # Full format option:
    # formatter = logging.Formatter('[%(asctime)s] %(levelname)s - %(name)s: %(message)s')
    console_handler.setFormatter(formatter)
    
    if not logger.hasHandlers():
        logger.addHandler(console_handler)
        
        # File handler (optional)
        if log_file:
            file_handler = logging.FileHandler(log_file)
            file_handler.setFormatter(formatter)
            logger.addHandler(file_handler)
    
    return logger

# Basic usage
logger = setup_logger("strainMinerGraph", verbosity_level=3, log_file="app.log")

logger.debug("This is a debug message.")
logger.info("This is an info message.")
logger.warning("This is a warning message.")
logger.error("This is an error message.")
logger.critical("This is a critical message.")

# Suppress third-party library logs
logging.getLogger("requests").setLevel(logging.INFO)
logging.getLogger("urllib3").setLevel(logging.INFO)

# Example: a request that will only log warnings or errors from 'requests'
import requests
try:
    logger.info("Attempting HTTP request...")
    requests.get("http://example.com/404")
except Exception as e:
    logger.error(f"Error during request: {e}")

# Verbosity examples
print("=== Level 0 (SILENT) ===")
logger0 = setup_logger("test0", verbosity_level=0)
logger0.debug("Debug invisible")
logger0.info("Info invisible") 
logger0.warning("Warning invisible")
logger0.error("Error invisible")
logger0.critical("Critical visible")

print("=== Level 1 (MINIMAL) ===")
logger1 = setup_logger("test1", verbosity_level=1)
logger1.debug("Debug invisible")
logger1.info("Info invisible") 
logger1.warning("Warning visible")
logger1.error("Error visible")
logger1.critical("Critical visible")

print("=== Level 2 (VERBOSE) ===")
logger2 = setup_logger("test2", verbosity_level=2)
logger2.debug("Debug invisible")
logger2.info("Info visible")
logger2.warning("Warning visible") 
logger2.error("Error visible")
logger2.critical("Critical visible")

print("=== Level 3 (DEBUG) ===")
logger3 = setup_logger("test3", verbosity_level=3)
logger3.debug("Debug visible")
logger3.info("Info visible")
logger3.warning("Warning visible")
logger3.error("Error visible")
logger3.critical("Critical visible")



# Usage: python logger_example.py > output.txt
# Or to capture stderr too: python logger_example.py > output.txt 2>&1