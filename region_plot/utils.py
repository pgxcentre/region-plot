
# This file is part of region-plot.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import logging
from subprocess import Popen, PIPE


def execute_command(name, command):
    """Executes a single command."""
    # Logging
    logging.info("Starting task '{}'".format(name))
    logging.debug("Command: {}".format(" ".join(command)))

    # Launching the command
    proc = Popen(command, stdout=PIPE, stderr=PIPE)

    # Waiting for the process to terminate
    outs, errs = proc.communicate()

    # Checking the return code
    rc = proc.returncode
    if rc != 0:
        logging.error("Task '{}' did not complete".format(name))
        return False

    # Everything when well
    logging.info("Task '{}' successful".format(name))

    return True
