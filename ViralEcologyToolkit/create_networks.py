import os


def test():
    return ('This is a test')


def checkFasta(fastaFile):
    """
    Checks to see if a fasta file is a valid protein file

    Args:
        fastaFile: The fasta file to check

    Raises:
        IncorrectFileError: if the file is not a Fasta file
        IOError: if the file can't be located.
    """
    if not os.path.isfile(fastaFile):
        raise IOError()

    return True


def allVsAllBLAST(fastafile):
    pass
