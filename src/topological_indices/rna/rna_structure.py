import subprocess
import logging

logger = logging.getLogger('rna')

def generate_structures(seq, number_of_structures=500):
    command = f'RNAsubopt -p {number_of_structures}'

    with subprocess.Popen(command.split(),
        stdin =subprocess.PIPE,
        stdout = subprocess.PIPE) as RNAsubopt:

        RNAsubopt.stdin.write((seq + "\n").encode())
        output, error = RNAsubopt.communicate()

    if error is not None:
        logger.error('RNAsubopt error! %s', error)
        return False

    structures = output.decode().split()[1:]

    return structures